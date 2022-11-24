from flask import render_template, flash, request, send_from_directory, redirect, url_for, Response, abort, Flask
from time import time
import requests
from src.molecule import Molecule
import os
import shutil
from multiprocessing import Pool
from src.SQEqp import calculate_charges
import numpy as np
import zipfile
from datetime import datetime
from numba import jit
from numba.core import types
from numba.typed import Dict
import json


@jit(nopython=True, cache=True) # paralelizovat to? Zkusit na velké struktuře
def precalculate_parameters_SQEqps(atomic_types, bonds_types, surfaces, parameters, bond_hardnesses):
    n_atoms = len(atomic_types)
    n_bonds = len(bonds_types)

    precalc_params = np.empty((n_atoms, 4), dtype=np.float32)
    precalc_bond_hardnesses = np.empty(n_bonds, dtype=np.float32)

    for i in range(n_atoms): # pokud nebudeme paralelizovat, tak nazipovat!


        symbol_i = atomic_types[i]
        surface = surfaces[i]

        electronegativity, hardness, width, q0, q0_cor, hardness_cor, electronegativity_cor, width_cor = parameters[symbol_i]
        precalc_params[i] = (-electronegativity + electronegativity_cor * surface,
                         hardness + hardness_cor * surface,
                         2 * (width + width_cor * surfaces[i]) ** 2,
                         q0 + q0_cor * surface)

    for i in range(n_bonds): # upravit podle toho, zda budeme paralelizovat
        precalc_bond_hardnesses[i] = bond_hardnesses[bonds_types[i]]

    return precalc_params, precalc_bond_hardnesses



@jit(nopython=True, cache=True) # paralelizovat to? Zkusit na velké struktuře
def precalculate_parameters_SQEqp(atomic_types, bonds_types, parameters, bond_hardnesses):
    n_atoms = len(atomic_types)
    n_bonds = len(bonds_types)
    precalc_params = np.empty((n_atoms, 4), dtype=np.float64)
    precalc_bond_hardnesses = np.empty(n_bonds, dtype=np.float64)
    for i in range(n_atoms): # pokud nebudeme paralelizovat, tak nazipovat!
        symbol_i = atomic_types[i]
        electronegativity, hardness, width, q0, = parameters[symbol_i]
        precalc_params[i] = np.array((-electronegativity,
                             hardness,
                             2 * width ** 2,
                             q0), dtype=np.float64)
    for i in range(n_bonds): # upravit podle toho, zda budeme paralelizovat
        precalc_bond_hardnesses[i] = bond_hardnesses[bonds_types[i]]
    return precalc_params, precalc_bond_hardnesses









root_dir = os.path.dirname(os.path.abspath(__file__))

try:
    os.mkdir(f"{root_dir}/calculated_structures")
except FileExistsError:
    pass

pdb_to_pqr_path = shutil.which("pdb2pqr30")




application = Flask(__name__)
application.jinja_env.trim_blocks = True
application.jinja_env.lstrip_blocks = True
application.config['SECRET_KEY'] = "asdfasdf"



def load_parameters():
    params_SQEqps = json.load(open(f"{root_dir}/parameters/parameters_SQEqps.json"))
    parameters_SQEqps = Dict.empty(key_type=types.unicode_type,
                                   value_type=types.float32[:])
    for key, value in params_SQEqps["atom"]["data"].items():
        parameters_SQEqps[key] = np.array(value, dtype=np.float32)
    bond_hardnesses_SQEqps = Dict.empty(key_type=types.unicode_type,
                                        value_type=types.float64)
    for key, value in params_SQEqps["bond"]["data"].items():
        bond_hardnesses_SQEqps[key] = value

    params_SQEqp = json.load(open(f"{root_dir}/parameters/parameters_SQEqp.json"))
    parameters_SQEqp = Dict.empty(key_type=types.unicode_type,
                                   value_type=types.float32[:])
    for key, value in params_SQEqp["atom"]["data"].items():
        parameters_SQEqp[key] = np.array(value, dtype=np.float32)
    bond_hardnesses_SQEqp = Dict.empty(key_type=types.unicode_type,
                                        value_type=types.float64)
    for key, value in params_SQEqp["bond"]["data"].items():
        bond_hardnesses_SQEqp[key] = value

    return parameters_SQEqp, bond_hardnesses_SQEqp, parameters_SQEqps, bond_hardnesses_SQEqps

parameters_SQEqp, bond_hardnesses_SQEqp, parameters_SQEqps, bond_hardnesses_SQEqps = load_parameters()


n_cpu = 1
# calculation time do flash?

# upgradovat věci od tomáše

# kontrola přímo v javascriptu

# zkontrolovat rychlost

# kontrola atomových typů

# dodělat organizmus a název proteinu

# upravit cachování ať nežere tolik paměti

 # updatovat python

 # change molecule.code na molecule.name
 
 # zkontrolovat, zda funguje rdkit exception!!

 # dát v results data o molekule a viev vedle sebe

def page_log(data_dir,
             step,
             log,
             delete_last_line=False):
    # předělat? kontrolovat pokud jsou na konci posledního řádku ..., tak smazat poslední řádek
    log = f"<p><span style='font-weight:bold'> Step {step}/6:</span> {log}</p>\n"
    logs = open(f"{data_dir}/page_log.txt").readlines()
    if delete_last_line:
        logs = "".join(logs[:-1]) + log
    else:
        logs = "".join(logs) + log
    with open(f"{data_dir}/page_log.txt", "w") as page_log_file:
        page_log_file.write(logs)



@application.route('/', methods=['GET', 'POST'])
def main_site():
    if request.method == 'POST':
        code = request.form['code']
        action = request.form["action"]
        if action == "settings":
            return render_template('settings.html',
                                   code=code)
        elif action == "calculate charges":
            ph = request.form['ph']
            basis_set = request.form['basis_set']
            structure_optimization = request.form['structure_optimization']
            prediction_version = request.form['prediction_version']


            try:
                ph = float(ph)
            except ValueError:
                flash('pH must be a float value!')
                return render_template('settings.html',
                                        code=code)
            if not 0 <= ph <= 14:
                flash('pH value must be between 1-14!')
                return render_template('settings.html',
                                        code=code)


            ID = f"{code}_{ph}"
            data_dir = f"{root_dir}/calculated_structures/{ID}"

            # check whether the structure with the given setting has already been calculated
            if os.path.isdir(f"{root_dir}/calculated_structures/{ID}"):
                return redirect(url_for('results',
                                        ID=ID,
                                        from_cache=True))
            os.mkdir(data_dir)
            os.mknod(f"{data_dir}/page_log.txt")

            s = time()
            response = requests.get(f"https://alphafold.ebi.ac.uk/files/AF-{code}-F1-model_v{prediction_version}.pdb")
            if response.status_code != 200:
                # shutil.rmtree(data_dir)
                os.system(f"rm -r {data_dir}")
                flash(f'No structure with UniProt code "{code}" in prediction version "{prediction_version}" found in the AlphaFold database.')
                return render_template('index.html')
            with open(f"{data_dir}/{code}.pdb", "w") as pdb_file:
                pdb_file.write(response.text)
            page_log(data_dir, 1, f"Structure downloaded. ({round(time() - s, 2)}s)")

            # start calculation
            return render_template('computation_progress.html',
                                   ID=ID,
                                   code=code,
                                   ph=ph)

    else:
        return render_template('index.html')



@application.route("/calculation")
def calculation():
    # download and save PDB file (and cif todo)
    ID = request.args.get("ID")
    code, ph = ID.split("_")
    data_dir = f"{root_dir}/calculated_structures/{ID}"
    with open(f"{root_dir}/logs.txt", "a") as log_file:
        log_file.write(f"{request.remote_addr} {code} {ph} {datetime.now().strftime('%d/%m/%Y %H:%M:%S')}\n")

    step_counter = 2
    page_log(data_dir,step_counter, "Protonation of structure...")
    s = time()
    pdb_file = f"{data_dir}/{code}.pdb"
    pdb_file_with_hydrogens = f"{pdb_file[:-4]}_added_H.pdb"
    os.system(f"{pdb_to_pqr_path} --log-level DEBUG --noopt --titration-state-method propka --with-ph {ph} "
              f"--pdb-output {pdb_file_with_hydrogens} {pdb_file} {pdb_file[:-4]}_added_H.pqr  > {data_dir}/propka.log 2>&1 ")
    # add root for obabel!!!!
    os.system(f"obabel -ipdb {pdb_file_with_hydrogens} -ommcif -O {data_dir}/{code}_added_H.cif")
    page_log(data_dir,step_counter, f"Structure protonated. ({round(time() - s, 2)}s)", delete_last_line=True)
    step_counter += 1


    page_log(data_dir,step_counter, "Loading of molecule...")
    s = time()
    try:
        molecule = Molecule(code, pdb_file_with_hydrogens)
    except ValueError as e:
        return redirect(url_for('wrong_structure',
                                ID=ID,
                                code=ID.split("_")[0],
                                message=str(e)))
    page_log(data_dir,step_counter, f"Molecule loaded. ({round(time() - s, 2)}s)", delete_last_line=True)
    step_counter += 1


    empirical_method = "SQEqp"
    if empirical_method == "SQEqp":
        page_log(data_dir,step_counter, "Precalculate parameters...")
        molecule.precalc_params, molecule.precalc_bond_hardnesses = precalculate_parameters_SQEqp(molecule.ats_srepr,
                                                                                            molecule.bonds_srepr,
                                                                                            parameters_SQEqp,
                                                                                            bond_hardnesses_SQEqp)
        page_log(data_dir, step_counter, f"Parameters precalculated. ({round(time() - s, 2)}s)",
                 delete_last_line=True)
        step_counter += 1


    elif empirical_method == "SQEqps":
        page_log(data_dir,step_counter, "Calculation of solvatable surface...")
        s = time()
        molecule.calculate_surfaces(cpu=n_cpu)
        page_log(data_dir, step_counter, f"Solvatable surface calculated. ({round(time() - s, 2)}s)",
                 delete_last_line=True)
        step_counter += 1

        page_log(data_dir,step_counter, "Precalculate parameters...")
        molecule.precalc_params, molecule.precalc_bond_hardnesses = precalculate_parameters_SQEqps(molecule.ats_srepr,
                                                                                            molecule.bonds_srepr,
                                                                                            molecule.surfaces,
                                                                                            parameters_SQEqps,
                                                                                            bond_hardnesses_SQEqps)
        page_log(data_dir, step_counter, f"Parameters precalculated. ({round(time() - s, 2)}s)",
                 delete_last_line=True)
        step_counter += 1



    page_log(data_dir,step_counter, "Creation of submolecules...")
    s = time()
    molecule.create_submolecules()
    page_log(data_dir,step_counter, f"Submolecules created. ({round(time() - s, 2)}s)", delete_last_line=True)
    step_counter += 1


    page_log(data_dir,step_counter, "Calculation of partial atomic charges...")
    s = time()
    with Pool(n_cpu) as p:
        all_charges = p.map(calculate_charges, [substructure for substructure in molecule.substructures])
    all_charges = [chg for chgs in all_charges for chg in chgs]
    all_charges -= (np.sum(all_charges) - molecule.total_chg) / len(all_charges)
    charges = all_charges
    with open(f"{data_dir}/charges.txt", "w") as chg_file:
        chg_file.write(f"{molecule.code}\n" + ' '.join([str(charge) for charge in charges]) + " \n")
    pqr_file_lines = open(f"{pdb_file[:-4]}_added_H.pqr").readlines()
    c = 0
    new_lines = []
    for line in pqr_file_lines:
        if line[:4] == "ATOM":
            new_lines.append(line[:54] + '{:>8.4f}'.format(charges[c])  + line[62:])
            c += 1
        else:
            new_lines.append(line)
    with open(f"{pdb_file[:-4]}_added_H.pqr", "w") as pqr_file:
        pqr_file.write("".join(new_lines))

    page_log(data_dir,step_counter, f"Partial atomic charges calculated. ({round(time() - s, 2)}s)", delete_last_line=True)
    return redirect(url_for('results',
                            ID=ID))




@application.route('/wrong_structure')
def wrong_structure():
    ID = request.args.get('ID')
    code, ph = ID.split("_")
    message = request.args.get("message")
    return render_template('wrong_structure.html',
                           code=code,
                           ID=ID,
                           message=message)


@application.route('/progress')
def progress():
    ID = request.args.get('ID')
    data_dir = f"{root_dir}/calculated_structures/{ID}"
    return open(f"{data_dir}/page_log.txt", "r").read()



@application.route('/results')
def results():
    ID = request.args.get('ID')
    data_dir = f"{root_dir}/calculated_structures/{ID}"
    code, ph = ID.split("_")
    try:
        absolute_charges = [abs(float(x)) for x in open(f"{data_dir}/charges.txt", "r").readlines()[1].split()]
    except FileNotFoundError:
        # shutil.rmtree(data_dir)
        os.system(f"rm {data_dir}")
        flash(f'Alphafold or pdb2pqr30 error with structure "{code}"!')
        return redirect(url_for('main_site'))

    chg_range = round(max(absolute_charges), 4)
    n_ats = len(absolute_charges)
    total_time = f"{round(sum([float(line.split('(')[1].split(')')[0][:-1]) for line in open(f'{data_dir}/page_log.txt').readlines()]), 2)} seconds"
    if request.args.get("from_cache"):
        total_time = f"{total_time} (The charges have already been calculated earlier and the results are taken from memory.)"


    return render_template('results.html',
                           ID=ID,
                           chg_range=chg_range,
                           code=code,
                           n_ats=n_ats,
                           total_time=total_time,
                           ph=ph)


@application.route('/download_wrong_structure')
def download_wrong_structure():
    ID = request.args.get("ID")
    code, _ = ID.split("_")
    data_dir = f"{root_dir}/calculated_structures/{ID}"
    with zipfile.ZipFile(f'{data_dir}/{ID}.zip', 'w') as zip:
        zip.write(f"{data_dir}/{code}.pdb", arcname=f"{code}.pdb")
        zip.write(f"{data_dir}/{code}_added_H.pdb", arcname=f"{code}_added_H.pdb")
    return send_from_directory(data_dir, f'{ID}.zip', as_attachment=True)


@application.route('/download')
def download_charges():
    ID = request.args.get("ID")
    code, _ = ID.split("_")
    data_dir = f"{root_dir}/calculated_structures/{ID}"
    with zipfile.ZipFile(f'{data_dir}/{ID}.zip', 'w') as zip:
        zip.write(f"{data_dir}/charges.txt", arcname=f"{code}_charges.txt")
        zip.write(f"{data_dir}/{code}_added_H.pdb", arcname=f"{code}_added_H.pdb")
        zip.write(f"{data_dir}/{code}_added_H.pqr", arcname=f"{code}_added_H.pqr")
        zip.write(f"{data_dir}/{code}_added_H.cif", arcname=f"{code}_added_H.cif")
    return send_from_directory(data_dir, f'{ID}.zip', as_attachment=True)


@application.route('/structure')
def get_structure():
    ID = request.args.get("ID")
    return Response(open(f"{root_dir}/calculated_structures/{ID}/{ID.split('_')[0]}_added_H.pdb", "r").read(),
                    mimetype='text/plain')


@application.route('/format')
def get_format():
    # todo remove it
    return Response("PDB",
                    mimetype='text/plain')


@application.route('/charges')
def get_charges():
    ID = request.args.get("ID")
    return Response(open(f"{root_dir}/calculated_structures/{ID}/charges.txt", "r").read(),
                    mimetype='text/plain')

