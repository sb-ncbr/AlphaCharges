from src.molecule import Molecule
from src.SQEqp_h import SQEqp_h, calculate_charges
import time
import sys
import numpy as np

s = time.time()
sss = time.time()
code = sys.argv[1]
molecule = Molecule(code, f"examples/structure_{code}_added_H.pdb")
print(f"Loading of molecule: {time.time() - s}")

molecule.calculate_surfaces()

empirical_method = SQEqp_h(f"parameters/parameters.json")


from numba import jit
@jit(nopython=True, cache=True) # paralelizovat to? Zkusit na velké struktuře
def precalculate_parameters(atomic_types, bonds_types, surfaces, parameters, bond_hardnesses):
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


s = time.time()
molecule.precalc_params, molecule.precalc_bond_hardnesses = precalculate_parameters(molecule.ats_srepr,
                                                                                    molecule.bonds_srepr,
                                                                                    molecule.surfaces,
                                                                                    empirical_method.parameters,
                                                                                    empirical_method.bond_hardnesses)
print(f"Parameters precalculated: {time.time() - s}")



s = time.time()
molecule.create_submolecules()
print(f"Create submolecules: {time.time() - s}")



s = time.time()
all_charges = []
from multiprocessing import Pool

with Pool(2) as p:
    charges = p.map(calculate_charges, [substructure for substructure in molecule.substructures])

all_charges = [chg for chgs in charges for chg in chgs]

# for substructure in molecule.substructures:
#     all_charges.extend(empirical_method.calculate_charges(substructure))


all_charges -= (np.sum(all_charges) - molecule.total_chg) / len(all_charges)
print(f"Calculate submolecules: {time.time() - s}")


#
# s = time.time()
# charges = empirical_method.calculate_charges(molecule)
# print(f"Calculation of charges: {time.time() - s}")
#
#
# #
# print(np.max(np.abs(all_charges-charges)))
# print(np.mean(np.abs(all_charges-charges)))
# print(np.sqrt(np.mean(np.abs(all_charges - charges) ** 2)))

#

