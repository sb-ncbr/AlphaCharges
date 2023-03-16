import gemmi
from . import SQEqp
from .molecule import Molecule
from .logs import Logs
import os
from datetime import datetime
from time import time
import requests
import json

class Calculation:
    def __init__(self,
                 ID: str,
                 remote_addr: str,
                 empirical_method: str,
                 root_dir: str):
        self.ID = ID
        self.empirical_method = empirical_method
        self.root_dir = root_dir
        self.code, self.ph, self.alphafold_prediction_version = self.ID.split('_')
        self.data_dir = f'{self.root_dir}/calculated_structures/{self.ID}'
        self.pdb_file = f'{self.data_dir}/{self.code}.pdb' # original pdb from alphafold, without hydrogens
        self.mmcif_file = f'{self.data_dir}/{self.code}.cif' # original mmCIF from alphafold, without hydrogens
        self.pdb_file_with_hydrogens = f'{self.data_dir}/{self.code}_added_H.pdb'
        self.pqr_file = f'{self.data_dir}/{self.code}.pqr'
        self.logs = Logs(data_dir=self.data_dir,
                         empirical_method=self.empirical_method)
        os.mkdir(self.data_dir)
        os.mknod(f'{self.data_dir}/page_log.txt')
        with open(f'{self.root_dir}/calculated_structures/logs.txt', 'a') as log_file:
            log_file.write(f'{remote_addr} {self.code} {self.ph} {self.alphafold_prediction_version} {datetime.now().strftime("%d/%m/%Y %H:%M:%S")}\n')

    def download_PDB(self):
        self.logs.add_log('Structure download...')
        s = time()
        response = requests.get(f'https://alphafold.ebi.ac.uk/files/AF-{self.code}-F1-model_v{self.alphafold_prediction_version}.pdb')
        with open(f'{self.pdb_file}', 'w') as pdb_file:
            pdb_file.write(response.text)
        response = requests.get(f'https://alphafold.ebi.ac.uk/files/AF-{self.code}-F1-model_v{self.alphafold_prediction_version}.cif')
        with open(f'{self.mmcif_file}', 'w') as mmcif_file:
            mmcif_file.write(response.text)
        self.logs.add_log(f'Structure downloaded. ({round(time() - s, 2)}s)')

    def protonate_structure(self):
        self.logs.add_log('Protonation of structure...')
        s = time()
        os.system(f'/opt/venv/bin/pdb2pqr30 --log-level DEBUG --noopt --titration-state-method propka '
                  f'--with-ph {self.ph} --pdb-output {self.pdb_file_with_hydrogens} {self.pdb_file} '
                  f'{self.pqr_file} > {self.data_dir}/propka.log 2>&1 ')
        self.logs.add_log(f'Structure protonated. ({round(time() - s, 2)}s)')

    def load_molecule(self):
        self.logs.add_log('Loading of molecule...')
        s = time()
        try:
            self.molecule = Molecule(self.pdb_file_with_hydrogens,
                                     self.pqr_file)
        except ValueError as error:
            with open(f"{self.data_dir}/problematic_atoms.json", "w") as problematic_atoms_file:
                json.dump(error.args[0], problematic_atoms_file)
            return False, str(error)
        self.logs.add_log(f'Molecule loaded. ({round(time() - s, 2)}s)')
        return True, None

    def precalculate_parameters(self):
        parameters, bond_hardnesses = SQEqp.load_parameters(self.root_dir,
                                                            self.empirical_method)
        if self.empirical_method == 'SQEqp':
            self.logs.add_log('Assigning parameters...')
            s = time()
            self.molecule.precalc_params,\
                self.molecule.precalc_bond_hardnesses = SQEqp.precalculate_parameters_SQEqp(self.molecule.ats_srepr,
                                                                                            self.molecule.bonds_srepr,
                                                                                            parameters,
                                                                                            bond_hardnesses)
            self.logs.add_log(f'Parameters assigned. ({round(time() - s, 2)}s)')

        elif self.empirical_method == 'SQEqps':
            self.logs.add_log('Calculation of solvatable surface...')
            s = time()
            self.molecule.calculate_surfaces(cpu=1)
            self.logs.add_log(f'Solvatable surface calculated. ({round(time() - s, 2)}s)')

            self.logs.add_log('Precalculate parameters...')
            s = time()
            self.molecule.precalc_params, \
                self.molecule.precalc_bond_hardnesses = SQEqp.precalculate_parameters_SQEqps(self.molecule.ats_srepr,
                                                                                             self.molecule.bonds_srepr,
                                                                                             self.molecule.surfaces,
                                                                                             parameters,
                                                                                             bond_hardnesses)
            self.logs.add_log(f'Parameters precalculated. ({round(time() - s, 2)}s)')

    def create_submolecules(self):
        self.logs.add_log('Creation of submolecules...')
        s = time()
        self.molecule.create_submolecules()
        self.logs.add_log(f'Submolecules created. ({round(time() - s, 2)}s)')

    def calculate_charges(self):
        self.logs.add_log('Calculation of partial atomic charges...')
        s = time()

        # calculation of charges
        # with Pool(n_cpu) as p:
        #    all_charges = p.map(calculate_charges, [substructure for substructure in molecule.substructures])
        # all_charges = [chg for chgs in all_charges for chg in chgs]
        all_charges = []
        for substructure in self.molecule.substructures:
            all_charges.extend(SQEqp.calculate_charges(substructure))
        all_charges -= (sum(all_charges) - self.molecule.total_chg) / len(all_charges)
        self.charges = all_charges

        # write charges to files
        self._write_txt()
        self._write_pqr()
        self._write_mmcif()

        self.logs.add_log(f'Partial atomic charges calculated. ({round(time() - s, 2)}s)')

    def _write_txt(self):
        with open(f'{self.data_dir}/charges.txt', 'w') as chg_file:
            chg_file.write(f'{self.code}\n' + ' '.join([str(round(charge, 4)) for charge in self.charges]) + ' \n')

    def _write_pqr(self):
        pqr_file_lines = open(self.pqr_file).readlines()
        c = 0
        new_lines = []
        for line in pqr_file_lines:
            if line[:4] == 'ATOM':
                new_lines.append(line[:54] + '{:>8.4f}'.format(self.charges[c]) + line[62:])
                c += 1
            else:
                new_lines.append(line)
        with open(self.pqr_file, 'w') as pqr_file:
            pqr_file.write(''.join(new_lines))

    def _add_AF_confidence_score(self, write_block):
        document = gemmi.cif.read(self.mmcif_file)
        block = document.sole_block()

        ma_qa_metric_prefix = '_ma_qa_metric'
        ma_qa_metric_local_prefix = '_ma_qa_metric_local'
        ma_qa_metric_global_prefix = '_ma_qa_metric_global'

        categories = {
            ma_qa_metric_prefix: block.get_mmcif_category(ma_qa_metric_prefix),
            ma_qa_metric_local_prefix: block.get_mmcif_category(ma_qa_metric_local_prefix),
            ma_qa_metric_global_prefix: block.get_mmcif_category(ma_qa_metric_global_prefix)
        }

        asym_id = write_block.get_mmcif_category('_struct_asym').get('id')[0]

        length = len(categories[ma_qa_metric_local_prefix]['label_asym_id'])
        categories[ma_qa_metric_local_prefix]['label_asym_id'] = [asym_id] * length

        for name, data in categories.items():
            write_block.set_mmcif_category(name, data)

    def _write_mmcif(self):
        input_file = self.pdb_file_with_hydrogens
        filename, _ = os.path.splitext(input_file)
        output_file = f"{filename}.cif"
        structure = gemmi.read_pdb(input_file)
        structure.setup_entities()
        structure.assign_label_seq_id()
        block = structure.make_mmcif_block()
        block.find_mmcif_category('_chem_comp.').erase() # remove pesky _chem_comp category >:(
        partial_atomic_charges_meta_prefix = "_partial_atomic_charges_meta."
        partial_atomic_charges_meta_attributes = ["id",
                                                  "type",
                                                  "method"]
        metadata_loop = block.init_loop(partial_atomic_charges_meta_prefix,
                                        partial_atomic_charges_meta_attributes)
        metadata_loop.add_row(['1',
                               "'empirical'",
                               "'SQE+qp/Schindler 2021 (PUB_pept)'"])
        partial_atomic_charges_prefix = "_partial_atomic_charges."
        partial_atomic_charges_attributes = ["type_id",
                                             "atom_id",
                                             "charge"]
        charges_loop = block.init_loop(partial_atomic_charges_prefix,
                                       partial_atomic_charges_attributes)
        for atomId, charge in enumerate(self.charges):
            charges_loop.add_row(["1",
                                  f"{atomId + 1}",
                                  f"{charge: .4f}"])

        self._add_AF_confidence_score(block)
        block.write_file(output_file)
