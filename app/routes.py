from flask import render_template, flash, request, send_from_directory, redirect, url_for, Response, abort
from time import time
from . import application
import requests
from .src.molecule import Molecule
from .src.SQEqp_h import SQEqp_h
import shutil

import os

request_data = {}

@application.route('/', methods=['GET', 'POST'])
def main_site():
    # todo jak čistit request_data? Jak dlouho by měly být výsledky dostupné?
    # todo za jak dlouho smazat tmp_disr?
    if request.method == 'POST':
        code = request.form['code']
        ph = request.form['ph']

        if not code:
            flash('Alphafold code is required!')
            return render_template('index.html')
        if ph:
            try:
                ph = float(ph)
            except ValueError:
                flash('pH must be a float value!')
                return render_template('index.html')
            if not 0 < ph < 14:
                flash('pH value must be between 1-14!')
                return render_template('index.html')
        if not ph:
            ph = 7.2

        s = time()
        response = requests.get(f"https://alphafold.ebi.ac.uk/files/AF-{code}-F1-model_v3.pdb")
        if response.status_code != 200:
            flash(f'Structure is not downloaded. Code "{code}" si probably wrong!')
            return render_template('index.html')
        print(f"Structure downloaded. ({time() - s})")
        n_heavy_atoms = response.text.count("ATOM")
        ID = f"{code}_{ph}"
        tmp_dir = f"calculate_{ID}"
        try:
            os.mkdir(tmp_dir)
        except FileExistsError:
            shutil.rmtree(tmp_dir)
            os.mkdir(tmp_dir)
        request_data[ID] = {}
        request_data[ID]["code"] = code
        request_data[ID]["ph"] = ph
        request_data[ID]["pdb_text"] = response.text
        request_data[ID]["tmp_dir"] = tmp_dir
        request_data[ID]["progress"] = f"<p><span style='font-weight:bold'> Step 1/6:</span> Structure downloaded. ({round(time() - s, 2)}s)</p>"

        if n_heavy_atoms < 500:
            return redirect(url_for('calculation',
                                    ID=ID))
        else:
            return render_template('computation_progress.html',
                                   ID=ID,
                                   code=code,
                                   ph=ph)

    else:
        return render_template('index.html')



@application.route("/calculation")
def calculation():
    ID = request.args.get("ID")
    calculation_data = request_data[ID]
    code = calculation_data["code"]
    ph = calculation_data["ph"]
    tmp_dir = calculation_data["tmp_dir"]
    pdb_text = calculation_data["pdb_text"]


    s = time()
    pdb_file = os.path.join(tmp_dir, f"structure_{code}.pdb")

    open(pdb_file, "w").write(pdb_text)
    pdb_file_with_hydrogens = f"{pdb_file[:-4]}_added_H.pdb"
    os.system(f"pdb2pqr30 --log-level DEBUG --noopt --with-ph {ph} "
              f"--pdb-output {pdb_file_with_hydrogens} {pdb_file} {pdb_file[:-4]}_added_H.pqr  > {tmp_dir}/propka.log 2>&1 ")
    print(f"Structure protonated. ({time() - s})")
    request_data[ID]["progress"] += f"<p><span style='font-weight:bold'> Step 2/6:</span> Structure protonated. ({round(time() - s, 2)}s) </p>"


    s = time()
    molecule = Molecule(code, pdb_file_with_hydrogens)
    print(f"Molecule loaded. ({time() - s})")
    request_data[ID]["progress"] += f"<p><span style='font-weight:bold'> Step 3/6:</span> Molecule loaded. ({round(time() - s, 2)}s) </p>"

    s = time()
    molecule.calculate_distace_matrix()
    print(f"Distance matrix calculated. ({time() - s})")
    request_data[ID]["progress"] += f"<p><span style='font-weight:bold'> Step 4/6:</span> Distance matrix calculated. ({round(time() - s, 2)}s) </p>"

    s = time()
    molecule.calculate_surfaces()
    print(f"Surface calculated. ({time() - s})")
    request_data[ID]["progress"] += f"<p><span style='font-weight:bold'> Step 5/6:</span> Solvatable surface caluclated. ({round(time() - s, 2)}s) </p>"

    s = time()
    empirical_method = SQEqp_h()
    request_data[ID]["pdb_file"] = open(pdb_file_with_hydrogens, "r").read()
    charges = empirical_method.calculate_charges(molecule)
    print(f"Charges calculated. ({time() - s})")
    request_data[ID]["progress"] += f"<p><span style='font-weight:bold'> Step 6/6:</span> Charges calculated. ({round(time() - s, 2)}s) </p>"

    request_data[ID]["charges"] = f"{molecule.code}\n" + ' '.join([str(charge) for charge in charges]) + " \n"
    request_data[ID]["n_ats"] = molecule.n_ats
    request_data[ID]["chg_range"] = round(max(abs(charges)), 4)
    return render_template('results.html',
                                   ID=ID)


@application.route('/progress')
def progress():
    ID = request.args.get('ID')
    return request_data[ID]["progress"]




@application.route('/results')
def results():
    ID = request.args.get('ID')
    results = request_data[ID]
    return render_template('results.html',
                           ID=ID,
                           chg_range=results['chg_range'],
                           code=results['code'],
                           n_ats=results['n_ats'],
                           ph=results['ph'])




@application.route('/download')
def download_charges():
    ID = request.args.get("ID")
    tmp_dir = request_data[ID]['tmpdir']
    with open(os.path.join(tmp_dir, "charges.txt"), "w") as chg_file:
        chg_file.write(request_data[ID]['charges'])

    return send_from_directory(tmp_dir, 'charges.txt', as_attachment=True)


@application.route('/structure')
def get_structure():
    ID = request.args.get("ID")
    return Response(request_data[ID]["pdb_file"], mimetype='text/plain')


@application.route('/format')
def get_format():
    # todo remove it
    return Response("PDB", mimetype='text/plain')


@application.route('/charges')
def get_charges():
    ID = request.args.get("ID")
    return Response(request_data[ID]["charges"], mimetype='text/plain')

