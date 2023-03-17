import json
import os
import zipfile
from random import random
from time import sleep, time

from flask import render_template, flash, request, send_from_directory, redirect, url_for, Response, Flask, Markup, jsonify

from src.calculation import Calculation
from src.input_validators import valid_pH, valid_prediction_version, valid_alphafold_request

application = Flask(__name__)
application.jinja_env.trim_blocks = True
application.jinja_env.lstrip_blocks = True
application.config['SECRET_KEY'] = str(random())
root_dir = os.path.dirname(os.path.abspath(__file__))


def already_calculated(ID):
    path = f'{root_dir}/calculated_structures/{ID}'
    if os.path.isdir(path):
        if os.path.isfile(f'{path}/charges.txt') or os.path.isfile(f'{path}/problematic_atoms.json'):
            return True
        elif time() - os.stat(path).st_mtime > 180:
            # for case that results directory exists without results (charges.txt or problematic_atoms.json)
            # it means that something unexpected happen during calculation
            os.system(f'rm -r {root_dir}/calculated_structures/{ID}')
    return False


def is_running(ID):
    path = f'{root_dir}/calculated_structures/{ID}'
    if os.path.isdir(path):
        if os.path.isfile(f'{path}/charges.txt') or os.path.isfile(f'{path}/problematic_atoms.json'):
            return False
        elif time() - os.stat(path).st_mtime > 180:
            return False
        return True
    return False

@application.route('/', methods=['GET', 'POST'])
def main_site():
    if request.method == 'POST':
        code = request.form['code'].strip().upper() # UniProt code, not case-sensitive
        code = code.replace("AF-","").replace("-F1", "") # Also AlphaFold DB identifiers are supproted (e.g., AF-A8H2R3-F1)

        if request.form['action'] == 'settings':
            return render_template('settings.html',
                                   code=code)

        elif request.form['action'] == 'calculate charges':
            ph, alphafold_prediction_version = request.form['ph'], request.form['prediction_version']

            ph, is_ph_valid = valid_pH(ph)
            if not is_ph_valid:
                message = 'pH must be a float value from 0 to 14!'
                flash(message, 'warning')
                return render_template('index.html',
                                       code=code)

            ID = f'{code}_{ph}_{alphafold_prediction_version}'

            # check whether the structure is currently calculated
            if is_running(ID):
                message = Markup(f'The partial atomic charges for your input are just calculated. '
                             f'For results visit  <a href="https://alphacharges.ncbr.muni.cz/results?ID={ID}" class="alert-link"'
                             f'target="_blank" rel="noreferrer">https://alphacharges.ncbr.muni.cz/results?ID={ID}</a>'
                             f' after a while.')
                flash(message, 'info')
                return render_template('index.html')

            if already_calculated(ID):
                return redirect(url_for('results',
                                        ID=ID))

            if not valid_alphafold_request(code, alphafold_prediction_version):
                message = Markup(f'The structure with code <strong>{code}</strong> in prediction version <strong>{alphafold_prediction_version}</strong> '
                      f'is either not found in AlphaFoldDB or the code is entered in the wrong format. '
                      f'UniProt code is allowed only in its short form (e.g., A0A1P8BEE7, B7ZW16). '
                      f'Other notations (e.g., A0A159JYF7_9DIPT, Q8WZ42-F2) are not supported. '
                      f'An alternative option is AlpfaFold DB Identifier (e.g., AF-L8BU87-F1).')
                flash(message, 'warning')
                return render_template('index.html')

            # start calculation
            return render_template('computation_progress.html',
                                   ID=ID,
                                   alphafold_prediction_version=alphafold_prediction_version,
                                   code=code,
                                   ph=ph)
    return render_template('index.html')


@application.route('/calculation', methods=['POST'])
def calculation():
    empirical_method = 'SQEqp'
    calculation = Calculation(request.args.get('ID'),
                              request.remote_addr,
                              empirical_method,
                              root_dir)
    calculation.download_PDB()
    calculation.protonate_structure()
    loaded, _ = calculation.load_molecule()
    if not loaded:
        return redirect(url_for('wrong_structure', ID=calculation.ID))
    calculation.precalculate_parameters()
    calculation.create_submolecules()
    calculation.calculate_charges()
    return redirect(url_for('results', ID=calculation.ID))


@application.route('/wrong_structure')
def wrong_structure():
    ID = request.args.get('ID')
    code, ph, alphafold_prediction_version = ID.split('_')
    problematic_atoms_file = open(f'{root_dir}/calculated_structures/{ID}/problematic_atoms.json', 'r')
    problematic_atoms = json.load(problematic_atoms_file)
    message = Markup('There is a structural error with atoms <span id="problematic_atoms"></span>! Calculation of partial atomic charges is not possible.')
    flash(message, 'danger')
    return render_template('wrong_structure.html',
                           ID=ID,
                           code=code,
                           ph=ph,
                           problematic_atoms=problematic_atoms,
                           alphafold_prediction_version=alphafold_prediction_version)


@application.route('/progress')
def progress():
    data_dir = f'{root_dir}/calculated_structures/{request.args.get("ID")}'
    return open(f'{data_dir}/page_log.txt', 'r').read()


@application.route('/results')
def results():
    ID = request.args.get('ID')

    try:
        code, ph, alphafold_prediction_version = ID.split('_')
    except:
        message = Markup('The ID was entered in the wrong format. The ID should be of the form <strong>&ltUniProt code&gt_&ltph&gt_&ltAlphaFold2 prediction version&gt</strong>.')
        flash(message, 'danger')
        return redirect(url_for('main_site'))

    data_dir = f'{root_dir}/calculated_structures/{ID}'

    if not already_calculated(ID) and not is_running(ID):
        message = Markup(f'There are no results for structure with UniProt <strong>{code}</strong> in AlphaFold2 prediction version <strong>{alphafold_prediction_version}</strong> and pH <strong>{ph}</strong>.')
        flash(message, 'danger')
        return redirect(url_for('main_site'))

    if os.path.isfile(f'{data_dir}/problematic_atoms.json'):
        return redirect(url_for('wrong_structure',
                                ID=ID))

    absolute_charges = [abs(float(x)) for x in open(f'{data_dir}/charges.txt', 'r').readlines()[1].split()]
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
    ID = request.args.get('ID')
    code = ID.split('_')[0]
    data_dir = f'{root_dir}/calculated_structures/{ID}'
    with zipfile.ZipFile(f'{data_dir}/{ID}.zip', 'w') as zip:
        zip.write(f'{data_dir}/{code}.pdb', arcname=f'{code}.pdb')
        zip.write(f'{data_dir}/{code}_added_H.pdb', arcname=f'{code}_added_H.pdb')
    return send_from_directory(data_dir, f'{ID}.zip', as_attachment=True)


@application.route('/download_files')
def download_files():
    ID = request.args.get('ID')
    code, _, _ = ID.split('_')
    data_dir = f'{root_dir}/calculated_structures/{ID}'
    with zipfile.ZipFile(f'{data_dir}/{ID}.zip', 'w') as zip:
        zip.write(f'{data_dir}/charges.txt', arcname=f'{code}_charges.txt')
        zip.write(f'{data_dir}/{code}_added_H.pdb', arcname=f'{code}_added_H.pdb')
        zip.write(f'{data_dir}/{code}.pqr', arcname=f'{code}.pqr')
        zip.write(f'{data_dir}/{code}_added_H.cif', arcname=f'{code}_added_H.cif')
    return send_from_directory(data_dir, f'{ID}.zip', as_attachment=True)


@application.route('/structure/<ID>/<FORMAT>')
def get_structure(ID: str,
                  FORMAT: str):
    filepath = f'{root_dir}/calculated_structures/{ID}/{ID.split("_")[0]}_added_H.{FORMAT}'
    return Response(open(filepath, 'r').read(), mimetype='text/plain')


@application.route('/calculate_charges/<string:code>')
def calculate_charges(code: str):
    # API
    code = code.upper()  # alphacharges should be not case sensitive
    empirical_method = 'SQEqp'
    message_dict = {'UniProt code': code,
                    'empirical method': empirical_method}
    allowed_url_arguments = set(['ph', 'alphafold_prediction_version'])
    if not set(request.args.keys()).issubset(allowed_url_arguments):
        message_dict.update({'status': 'failed',
                             'error message': 'Only URL arguments "ph" and "alphafold_prediction_version" are allowed.'})
        return jsonify(message_dict), 400

    ph, is_ph_valid = valid_pH(request.args.get('ph'))
    message_dict['pH'] = ph
    if not is_ph_valid:
        message_dict.update({'status': 'failed',
                             'error message': 'pH must be a float value from 0 to 14!'})
        return jsonify(message_dict), 400

    alphafold_prediction_version, is_version_valid = valid_prediction_version(request.args.get('alphafold_prediction_version'))
    message_dict['AlphaFold2 prediction version'] = alphafold_prediction_version
    if not is_version_valid:
        message_dict.update({'status': 'failed',
                             'error message': 'AlphaFold2 prediction version can be integer from 1 to 4'})
        return jsonify(message_dict), 400

    ID = f'{code}_{ph}_{alphafold_prediction_version}'
    if is_running(ID):
        while is_running(ID):
            sleep(5)
    else:
        if not already_calculated(ID):
            if not valid_alphafold_request(code, alphafold_prediction_version):
                message_dict.update({'status': 'failed',
                                     'error message': f'The structure with UniProt code {code} in prediction version {alphafold_prediction_version} '
                                                      f'is either not found in AlphaFoldDB or the UniProt code is entered in the wrong format. '
                                                      f'UniProt code is allowed only in its short form (e.g., A0A1P8BEE7, B7ZW16). '
                                                      f'Other notations (e.g., A0A159JYF7_9DIPT, Q8WZ42-F2) are not supported.'})
                return jsonify(message_dict), 400
            calculation = Calculation(ID,
                                      request.remote_addr,
                                      empirical_method,
                                      root_dir)
            calculation.download_PDB()
            calculation.protonate_structure()
            loaded, problematic_atoms = calculation.load_molecule()
            if not loaded:
                message_dict.update({'status': 'failed',
                                     'error message': f'There is an error with atoms <strong>{problematic_atoms}</strong>! '
                                                      f'The structure is probably incorrectly predicted by AlphaFold2, or incorrectly protonated by PROPKA3.'})
                return jsonify(message_dict), 501
            calculation.precalculate_parameters()
            calculation.create_submolecules()
            calculation.calculate_charges()
    message_dict.update({'status': 'partial atomic charges successfully calculated',
                         'ID': ID})
    return jsonify(message_dict)


@application.route('/download_file/<string:ID>/<string:format>')
def download_file(ID: str,
                  format: str):
    if not already_calculated(ID):
        return Response(f'No results calculated for this ID.',
                        status=400)
    code = ID.split('_')[0]
    data_dir = f'{root_dir}/calculated_structures/{ID}'
    if format == 'txt':
        file = 'charges.txt'
    elif format == 'pdb':
        file = f'{code}_added_H.pdb'
    elif format == 'pqr':
        file = f'{code}.pqr'
    elif format == 'mmcif':
        file = f'{code}_added_H.cif'
    else:
        return Response('Wrong format. Supported formats are txt, pdb, pqr and mmcif.',
                        status=400)
    return send_from_directory(data_dir, file, as_attachment=True)


@application.errorhandler(404)
def page_not_found(error):
    return render_template('404.html'), 404
