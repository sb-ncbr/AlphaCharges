import requests

from src.molecule import Molecule
from src.SQEqp_h import SQEqp_h


import sys
import os

from time import time

def main_site():
    code = sys.argv[1]
    ph = 7.0


    s = time()
    response = requests.get(f"https://alphafold.ebi.ac.uk/files/AF-{code}-F1-model_v3.pdb")
    if response.status_code != 200:
        exit("wrong code")
    print(f"Structure downloaded. ({time() - s})")


    s = time()
    tmp_dir = f"calculation_{code}_{str(time())}"
    os.system(f"mkdir {tmp_dir}")
    pdb_file = os.path.join(tmp_dir, f"structure_{code}.pdb")
    open(pdb_file, "w").write(response.text)
    pdb_file_with_hydrogens = f"{pdb_file[:-4]}_added_H.pdb"
    os.system(f"pdb2pqr30 --log-level DEBUG --noopt --with-ph {ph} "
              f"--pdb-output {pdb_file_with_hydrogens} {pdb_file} {pdb_file[:-4]}_added_H.pqr > {tmp_dir}/propka.log 2>&1 ")
    print(f"Structure protonated. ({time() - s})")


    s = time()
    molecule = Molecule(code, pdb_file_with_hydrogens)
    print(f"Molecule loaded. ({time() - s})")


    s = time()
    molecule.calculate_distace_matrix()
    print(f"Distance matrix calculated. ({time() - s})")


    s = time()
    molecule.calculate_surfaces()
    print(f"Surface calculated. ({time() - s})")
    exit()

    s = time()

    empirical_method = SQEqp_h()
    print(f"Empirical method loaded. ({time() - s})")

    s = time()
    charges = empirical_method.calculate_charges(molecule)
    print(f"Charges calculated. ({time() - s})")





main_site()







