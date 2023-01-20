import os
import requests
import zipfile
from time import time, sleep
from datetime import datetime
from flask import render_template, flash, request, send_from_directory, redirect, url_for, Response, Flask, Markup, jsonify
from src.SQEqp import calculate_charges, load_parameters, precalculate_parameters_SQEqp, precalculate_parameters_SQEqps
from src.molecule import Molecule
from random import random

application = Flask(__name__)
application.jinja_env.trim_blocks = True
application.jinja_env.lstrip_blocks = True
application.config['SECRET_KEY'] = str(random())

root_dir = os.path.dirname(os.path.abspath(__file__))

parameters_SQEqp, bond_hardnesses_SQEqp, parameters_SQEqps, bond_hardnesses_SQEqps = load_parameters(root_dir)

currently_running = set()

class Logs:
    def __init__(self,
                 data_dir: str,
                 empirical_method: str):
        self.data_dir = data_dir
        self.step = 1
        self.max_steps = 6 if empirical_method == "SQEqp" else 7

    def add_log(self,
                log: str):
        html_log = f"<p><span style='font-weight:bold'> Step {self.step}/{self.max_steps}:</span> {log}</p>\n"
        previous_logs = open(f"{self.data_dir}/page_log.txt").readlines()
        with open(f"{self.data_dir}/page_log.txt", "w") as page_log_file:
            if len(previous_logs) and "..." in previous_logs[-1]:
                previous_logs = previous_logs[:-1]
                self.step += 1
            page_log_file.write("".join(previous_logs) + html_log)


def valid_pH(ph):
    try:
        ph = float(ph)
    except ValueError:
        return False
    if not 0 <= ph <= 14:
        return False
    return True

def is_calculated(ID):
    # check whether the structure with the given setting has already been calculated
    if os.path.isdir(f"{root_dir}/calculated_structures/{ID}"):
        if os.path.isfile(f"{root_dir}/calculated_structures/{ID}/charges.txt"):
            return True
        else:  # for case that results directory exists without results (e.g. charges.txt)
            os.system(f"rm -r {root_dir}/calculated_structures/{ID}")
    return False

def is_valid_alphafold_request(code, alphafold_prediction_version):
    # check whether UniProt code is valid, ping AlphaFold website
    response = requests.head(f"https://alphafold.ebi.ac.uk/files/AF-{code}-F1-model_v{alphafold_prediction_version}.pdb")
    if response.status_code != 200:
        return False
    return True



@application.route('/', methods=['GET', 'POST'])
def main_site():

    if request.method == 'POST':

        code = request.form['code'] # UniProt code
        action = request.form["action"]

        if action == "settings":
            return render_template('settings.html',
                                   code=code)

        elif action == "calculate charges":

            ph = request.form['ph']
            alphafold_prediction_version = request.form['prediction_version']

            if not valid_pH(ph):
                flash("Error! pH must be a float value from 0 to 14!")
                return render_template('index.html',
                                       code=code)

            ID = f"{code}_{ph}_{alphafold_prediction_version}"

            # check whether the structure is currently calculated
            if ID in currently_running:
                flash(Markup(f'The partial atomic charges for your input are just calculated. '
                             f'For results visit  <a href="https://alphacharges.ncbr.muni.cz/results?ID={ID}" '
                             f'target="_blank" rel="noreferrer">https://alphacharges.ncbr.muni.cz/results?ID={ID}</a>'
                             f' after a while.'))
                return render_template('index.html')

            if is_calculated(ID):
                return redirect(url_for('results',
                                        ID=ID))

            if not is_valid_alphafold_request(code, alphafold_prediction_version):
                flash(f'The structure with UniProt code {code} in prediction version {alphafold_prediction_version} '
                      f'is either not found in AlphaFoldDB or the UniProt code is entered in the wrong format. '
                      f'UniProt code is allowed only in its short form (e.g., A0A1P8BEE7, B7ZW16). '
                      f'Other notations (e.g., A0A159JYF7_9DIPT, Q8WZ42-F2) are not supported.')
                return render_template('index.html')

            # start calculation
            return render_template('computation_progress.html',
                                   ID=ID,
                                   alphafold_prediction_version=alphafold_prediction_version,
                                   code=code,
                                   ph=ph)

    else:
        return render_template('index.html')


class Calculation:
    def __init__(self,
                 ID,
                 remote_addr,
                 empirical_method: str):
        self.ID = ID
        self.empirical_method = empirical_method
        self.code, self.ph, self.alphafold_prediction_version = self.ID.split("_")
        self.data_dir = f"{root_dir}/calculated_structures/{self.ID}"
        self.pdb_file = f"{self.data_dir}/{self.code}.pdb" # original pdb from alphafold, without hydrogens
        self.pdb_file_with_hydrogens = f"{self.data_dir}/{self.code}_added_H.pdb"
        self.pqr_file = f"{self.data_dir}/{self.code}.pqr"
        self.logs = Logs(data_dir=self.data_dir,
                         empirical_method=self.empirical_method)
        currently_running.update([self.ID])
        os.mkdir(self.data_dir)
        os.mknod(f"{self.data_dir}/page_log.txt")
        with open(f"{root_dir}/logs.txt", "a") as log_file:
            log_file.write(f"{remote_addr} {self.code} {self.ph} {self.alphafold_prediction_version} {datetime.now().strftime('%d/%m/%Y %H:%M:%S')}\n")

    def download_PDB(self):
        self.logs.add_log("Structure download...")
        s = time()
        response = requests.get(f"https://alphafold.ebi.ac.uk/files/AF-{self.code}-F1-model_v{self.alphafold_prediction_version}.pdb")
        with open(f"{self.pdb_file}", "w") as pdb_file:
            pdb_file.write(response.text)
        self.logs.add_log(f"Structure downloaded. ({round(time() - s, 2)}s)")

    def protonate_structure(self):
        self.logs.add_log("Protonation of structure...")
        s = time()
        os.system(f"/opt/venv/bin/pdb2pqr30 --log-level DEBUG --noopt --titration-state-method propka "
                  f"--with-ph {self.ph} --pdb-output {self.pdb_file_with_hydrogens} {self.pdb_file} "
                  f"{self.pqr_file} > {self.data_dir}/propka.log 2>&1 ")
        self.logs.add_log(f"Structure protonated. ({round(time() - s, 2)}s)")

    def load_molecule(self):
        self.logs.add_log("Loading of molecule...")
        s = time()
        try:
            self.molecule = Molecule(self.pdb_file_with_hydrogens,
                                     self.pqr_file)
        except ValueError as e:
            return False, str(e)
        self.logs.add_log(f"Molecule loaded. ({round(time() - s, 2)}s)")
        return True, None

    def precalculate_parameters(self):
        if self.empirical_method == "SQEqp":
            self.logs.add_log("Assigning parameters...")
            s = time()
            self.molecule.precalc_params, \
                self.molecule.precalc_bond_hardnesses = precalculate_parameters_SQEqp(self.molecule.ats_srepr,
                                                                                      self.molecule.bonds_srepr,
                                                                                      parameters_SQEqp,
                                                                                      bond_hardnesses_SQEqp)
            self.logs.add_log(f"Parameters assigned. ({round(time() - s, 2)}s)")

        elif self.empirical_method == "SQEqps":
            self.logs.add_log("Calculation of solvatable surface...")
            s = time()
            self.molecule.calculate_surfaces(cpu=1)
            self.logs.add_log(f"Solvatable surface calculated. ({round(time() - s, 2)}s)")

            self.logs.add_log("Precalculate parameters...")
            s = time()
            self.molecule.precalc_params, \
                self.molecule.precalc_bond_hardnesses = precalculate_parameters_SQEqps(self.molecule.ats_srepr,
                                                                                       self.molecule.bonds_srepr,
                                                                                       self.molecule.surfaces,
                                                                                       parameters_SQEqps,
                                                                                       bond_hardnesses_SQEqps)
            self.logs.add_log(f"Parameters precalculated. ({round(time() - s, 2)}s)")

    def create_submolecules(self):
        self.logs.add_log("Creation of submolecules...")
        s = time()
        self.molecule.create_submolecules()
        self.logs.add_log(f"Submolecules created. ({round(time() - s, 2)}s)")

    def calculate_charges(self):
        self.logs.add_log("Calculation of partial atomic charges...")
        s = time()

        # calculation of charges
        # with Pool(n_cpu) as p:
        #    all_charges = p.map(calculate_charges, [substructure for substructure in molecule.substructures])
        # all_charges = [chg for chgs in all_charges for chg in chgs]
        all_charges = []
        for substructure in self.molecule.substructures:
            all_charges.extend(calculate_charges(substructure))
        all_charges -= (sum(all_charges) - self.molecule.total_chg) / len(all_charges)
        charges = all_charges

        # writing charges to txt
        with open(f"{self.data_dir}/charges.txt", "w") as chg_file:
            chg_file.write(f"{self.code}\n" + ' '.join([str(round(charge, 4)) for charge in charges]) + " \n")

        # writing charges to pqr
        pqr_file_lines = open(self.pqr_file).readlines()
        c = 0
        new_lines = []
        for line in pqr_file_lines:
            if line[:4] == "ATOM":
                new_lines.append(line[:54] + '{:>8.4f}'.format(charges[c]) + line[62:])
                c += 1
            else:
                new_lines.append(line)
        with open(self.pqr_file, "w") as pqr_file:
            pqr_file.write("".join(new_lines))

        # writing charges to mmcif
        os.system(f"gemmi convert {self.pdb_file_with_hydrogens} {self.data_dir}/{self.code}_added_H.cif")
        # mmcif_lines_chgs = []
        # c = 0
        # for line in open(f"{data_dir}/{code}_added_H.cif", "r").readlines():
        #     sl = line.split()
        #     if len(sl) > 2 and sl[0].isdigit() and sl[1] in "HCNOS":
        #         sl[-4] = str(round(charges[c], 4))
        #         c += 1
        #         mmcif_lines_chgs.append(" ".join(sl) + "\n")
        #     else:
        #         mmcif_lines_chgs.append(line)
        # with open(f"{data_dir}/{code}_added_H.cif", "w") as mmcif_file:
        #     mmcif_file.write("".join(mmcif_lines_chgs))
        self.logs.add_log(f"Partial atomic charges calculated. ({round(time() - s, 2)}s)")
        currently_running.remove(self.ID)


@application.route("/calculation", methods=['POST'])
def calculation():
    empirical_method = "SQEqp"
    calculation = Calculation(request.args.get("ID"),
                              request.remote_addr,
                              empirical_method)
    calculation.download_PDB()
    calculation.protonate_structure()
    loaded, problematic_atom = calculation.load_molecule()
    if not loaded:
        return redirect(url_for('wrong_structure',
                                ID=calculation.ID,
                                code=calculation.code,
                                message=problematic_atom))
    calculation.precalculate_parameters()
    calculation.create_submolecules()
    calculation.calculate_charges()
    return redirect(url_for('results',
                            ID=calculation.ID))

@application.route('/wrong_structure')
def wrong_structure():
    ID = request.args.get('ID')
    message = request.args.get("message")
    return render_template('wrong_structure.html',
                           code=ID.split("_")[0],
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

    if ID is None:
        flash(f'No ID parameter was entered in your request! The ID should be of the form <UniProt code>_<ph>_<AlphaFold2 prediction version>.')
        return redirect(url_for('main_site'))

    try:
        code, ph, alphafold_prediction_version = ID.split("_")
    except:
        flash(f'The ID was entered in the wrong format. The ID should be of the form <UniProt code>_<ph>_<AlphaFold2 prediction version>.')
        return redirect(url_for('main_site'))

    data_dir = f"{root_dir}/calculated_structures/{ID}"

    try:
        absolute_charges = [abs(float(x)) for x in open(f"{data_dir}/charges.txt", "r").readlines()[1].split()]
    except FileNotFoundError:
        os.system(f"rm -r {data_dir}")
        flash(f'There are no results for structure with UnitProt {code} in AlphaFold2 prediction version {alphafold_prediction_version} and pH {ph}.')
        return redirect(url_for('main_site'))

    chg_range = round(max(absolute_charges), 4)
    n_ats = len(absolute_charges)

    return render_template('results.html',
                           ID=ID,
                           chg_range=chg_range,
                           code=code,
                           n_ats=n_ats,
                           ph=ph,
                           alphafold_prediction_version=alphafold_prediction_version)


@application.route('/download_files_for_wrong_structure')
def download_files_for_wrong_structure():
    ID = request.args.get("ID")
    code = ID.split("_")[0]
    data_dir = f"{root_dir}/calculated_structures/{ID}"
    with zipfile.ZipFile(f'{data_dir}/{ID}.zip', 'w') as zip:
        zip.write(f"{data_dir}/{code}.pdb", arcname=f"{code}.pdb")
        zip.write(f"{data_dir}/{code}_added_H.pdb", arcname=f"{code}_added_H.pdb")
    return send_from_directory(data_dir, f'{ID}.zip', as_attachment=True)


@application.route('/download_files')
def download_files():
    ID = request.args.get("ID")
    code, _, _ = ID.split("_")
    data_dir = f"{root_dir}/calculated_structures/{ID}"
    with zipfile.ZipFile(f'{data_dir}/{ID}.zip', 'w') as zip:
        zip.write(f"{data_dir}/charges.txt", arcname=f"{code}_charges.txt")
        zip.write(f"{data_dir}/{code}_added_H.pdb", arcname=f"{code}_added_H.pdb")
        zip.write(f"{data_dir}/{code}.pqr", arcname=f"{code}.pqr")
        zip.write(f"{data_dir}/{code}_added_H.cif", arcname=f"{code}_added_H.cif")
    return send_from_directory(data_dir, f'{ID}.zip', as_attachment=True)


@application.route('/structure')
def get_structure():
    ID = request.args.get("ID")
    return Response(open(f"{root_dir}/calculated_structures/{ID}/{ID.split('_')[0]}_added_H.pdb", "r").read(),
                    mimetype='text/plain')


@application.route('/format')
def get_format():
    return Response("PDB",
                    mimetype='text/plain')


@application.route('/charges')
def get_charges():
    ID = request.args.get("ID")
    return Response(open(f"{root_dir}/calculated_structures/{ID}/charges.txt", "r").read(),
                    mimetype='text/plain')


@application.route("/get_sqeqp_charges")
def get_sqeqp_charges():
    # API
    empirical_method = "SQEqp"
    code = request.args.get('code')
    ph = request.args.get('ph')
    alphafold_prediction_version = request.args.get('alphafold_prediction_version')
    if any([argument is None for argument in [code, ph, alphafold_prediction_version]]):
        return Response("Url arguments <code>, <ph> and <alphafold_prediction_version> are required.",
                        status=400)
    if not valid_pH(ph):
        return Response("Error! pH must be a float value from 0 to 14!",
                        status=400)

    ID = f"{code}_{ph}_{alphafold_prediction_version}"
    if ID in currently_running:
        while ID in currently_running:
            sleep(1)
    else:
        if not is_calculated(ID):
            if not is_valid_alphafold_request(code, alphafold_prediction_version):
                return Response(f'The structure with UniProt code {code} in prediction version {alphafold_prediction_version} '
                                f'is either not found in AlphaFoldDB or the UniProt code is entered in the wrong format. '
                                f'UniProt code is allowed only in its short form (e.g., A0A1P8BEE7, B7ZW16). '
                                f'Other notations (e.g., A0A159JYF7_9DIPT, Q8WZ42-F2) are not supported.',
                                status=400)
            calculation = Calculation(ID,
                                      request.remote_addr,
                                      empirical_method)
            calculation.download_PDB()
            calculation.protonate_structure()
            loaded, problematic_atom = calculation.load_molecule()
            if loaded is False:
                return Response(f"There is an error with atom with index {problematic_atom}!"
                                f" The structure is probably incorrectly predicted by AlphaFold2, or incorrectly protonated by PROPKA3.",
                                status=400)
            calculation.precalculate_parameters()
            calculation.create_submolecules()
            calculation.calculate_charges()
    return jsonify({"UniProt code": code,
                    "pH": ph,
                    "AlphaFold2 prediction version": alphafold_prediction_version,
                    "ID": ID,
                    "empirical method": empirical_method})



@application.route('/download_file')
def download_file():
    ID = request.args.get("ID")
    code = ID.split("_")[0]
    format = request.args.get("format")
    if any([argument is None for argument in [ID, format]]):
        return Response("Url arguments <ID> and <format> are required. Supproted formats are txt, pdb, pqr and mmcif.",
                        status=400)
    data_dir = f"{root_dir}/calculated_structures/{ID}"
    if format == "txt":
        file = "charges.txt"
    elif format == "pdb":
        file = f"{code}_added_H.pdb"
    elif format == "pqr":
        file = f"{code}.pqr"
    elif format == "mmcif":
        file = f"{code}_added_H.cif"
    else:
        return Response("Wrong format. Supproted formats are txt, pdb, pqr and mmcif.",
                        status=400)
    return send_from_directory(data_dir, file, as_attachment=True)
