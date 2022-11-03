import sys
from io import StringIO
from multiprocessing import Pool
import numpy as np
from numba import jit
from numba.typed import List
from rdkit import Chem
from sklearn.neighbors import KDTree as kdtreen


class Substructure:
    def __init__(self,
                 residues: list):

        self.residues = residues
        self.calculated_atoms = len(residues[0].coordinates)
        self.surfaces = np.concatenate([res.surfaces for res in residues])
        self.coordinates = np.concatenate([res.coordinates for res in residues], axis=0)
        self.precalc_params = np.concatenate([res.precalc_params for res in residues], axis=0)
        self.total_chg = sum([chg for res in residues for chg in res.propka_charges])
        precalc_bond_hardnesses = [hardness for res in residues for hardness in res.precalc_bond_hardnesses]
        bonds = [bond for res in residues for bond in res.bonds]

        # add bonds between residues (peptide and disulfide bonds)
        residues_numbers = set([res.number for res in residues])
        for res in residues:
            for connected_res, connected_by, bond_hardness in zip(res.connected_with,
                                                                  res.connected_bonds,
                                                                  res.connected_precalc_bond_hardnesses):
                if connected_res in residues_numbers:
                    bonds.append(connected_by)
                    precalc_bond_hardnesses.append(bond_hardness)

        # reindexing
        reindexes = {}
        c = 0
        for res in residues:
            for idx in res.indices:
                reindexes[idx] = c
                c += 1
        reindexed_bonds = []
        for a1, a2, bond_type in bonds:
            reindexed_bonds.append([reindexes[a1], reindexes[a2], bond_type])
        self.bonds = np.array(reindexed_bonds, dtype=np.int32)
        self.precalc_bond_hardnesses = np.array(precalc_bond_hardnesses)


class Residue:
    def __init__(self,
                 name: str,
                 number: int,
                 coordinates: np.array,
                 surfaces: np.array,
                 propka_charges: np.array,
                 indices: np.array,
                 precalc_params: np.array):

        self.name = name
        self.number = number
        self.coordinates = coordinates
        self.coordinates_mean = np.mean(self.coordinates, axis=0)
        self.surfaces = surfaces
        self.propka_charges = propka_charges
        self.indices = indices
        self.precalc_params = precalc_params


class Molecule:
    def __init__(self,
                 code: str,
                 pdb_file: str):

        self.code = code

        # load charges from propka
        self.propka_charges = [float(line.split()[8]) for line in
                               open(pdb_file.replace(".pdb", ".pqr"), "r").readlines()[:-2]]
        self.total_chg = round(sum(self.propka_charges))

        # load molecule by rdkit
        Chem.WrapLogs()
        terminal_stdout = sys.stderr

        sio = sys.stderr = StringIO()
        self.rdkit_mol = Chem.MolFromPDBFile(pdb_file,
                                             removeHs=False,
                                             sanitize=True)
        if self.rdkit_mol is None:
            raise ValueError(f"{sio.getvalue().split()[6]}")
        sys.stderr = terminal_stdout
        self.symbols = [atom.GetSymbol() for atom in self.rdkit_mol.GetAtoms()]
        self.n_ats = len(self.symbols)
        self.calculated_atoms = self.n_ats
        bonds = []
        bond_types = {"SINGLE": 1,
                      "DOUBLE": 2,
                      "TRIPLE": 3,
                      "AROMATIC": 4}
        for bond in self.rdkit_mol.GetBonds():
            a1 = bond.GetBeginAtom().GetIdx()
            a2 = bond.GetEndAtom().GetIdx()
            bond_type = bond_types[str(bond.GetBondType())]
            if a1 < a2:
                bonds.append((a1, a2, bond_type))
            else:
                bonds.append((a2, a1, bond_type))
        self.bonds = np.array(bonds, dtype=np.int32)
        self.num_of_bonds = len(self.bonds)
        coordinates = []
        for i in range(0, self.rdkit_mol.GetNumAtoms()):
            pos = self.rdkit_mol.GetConformer().GetAtomPosition(i)
            coordinates.append((pos.x, pos.y, pos.z))
        self.coordinates = np.array(coordinates, dtype=np.float32)



        ats_sreprba = self.create_ba()
        ats_sreprba2 = self.create_ba2()
        bonds_srepr = [f"{'-'.join(sorted([ats_sreprba[ba1], ats_sreprba[ba2]]))}-{bond_type}"
                       for ba1, ba2, bond_type in bonds]

        # aromaticity of HIP
        bonds_srepr_modified = []
        for bond in bonds_srepr:
            if bond not in ["C/HNN-N/CCH-1", "C/CCN-N/CCH-1", "C/CHN-N/CCH-1", "C/CCN-C/CHN-2"]:
                bonds_srepr_modified.append(bond)
            else:
                bonds_srepr_modified.append(bond[:-1] + "4")
        bonds_srepr = bonds_srepr_modified
        # aromaticity of HIP

        real_ats_types = {'C/CCC', 'C/CCCH', 'C/CCH', 'C/CCHH', 'C/CCHN', 'C/CCHO', 'C/CCN', 'C/CCO', 'C/CHHH',
                          'C/CHHN', 'C/CHHO', 'C/CHHS',
                          'C/CHN', 'C/CNO', 'C/COO', 'C/HHHS', 'C/HNN', 'C/NNN', 'H/C', 'H/N', 'H/O', 'H/S', 'N/CC',
                          'N/CCC', 'N/CCH',
                          'N/CCHH', 'N/CHH', 'N/CHHH', 'O/C', 'O/CH', 'S/CC', 'S/CH', 'S/CS'}

        # baex
        ats_srepr = []
        for i, (atba, atba2) in enumerate(zip(ats_sreprba, ats_sreprba2)):
            if atba not in real_ats_types:
                raise ValueError(f"{i+1} {atba}")

            if atba in ['H/N', "H/C", "O/C", "N/CHH", "N/CHHH", "C/CHHS", "C/CCC"]:
                ats_srepr.append(atba2)
            elif atba == "O/CH":
                if atba2 == "O/CH/CC":
                    ats_srepr.append("O/CH1")
                elif atba2 in ["O/CH/CHH", "O/CH/CCH"]:
                    ats_srepr.append("O/CH2")
                else:
                    raise ValueError(f"{i+1} {atba}")
            elif atba == "N/CCH":
                if len(atba2.split("/")[2]) == 4:
                    ats_srepr.append("N/CCH1")
                elif len(atba2.split("/")[2]) == 5:
                    ats_srepr.append("N/CCH2")
                elif len(atba2.split("/")[2]) == 3:
                    ats_srepr.append("N/CCH3")
                else:
                    raise ValueError(f"{i+1} {atba}")
            elif atba == "C/CHN":
                if atba2 == "C/CHN/CCCH":
                    ats_srepr.append("C/CHN1")
                elif atba2 in ["C/CHN/CCN", "C/CHN/CCHN"]:
                    ats_srepr.append("C/CHN2")
                else:
                    raise ValueError(f"{i+1} {atba}")
            elif atba == "C/CHHN":
                if atba2 == "C/CHHN/CHOO":
                    ats_srepr.append("C/CHHN1")
                elif atba2 in ["C/CHHN/CHNO", "C/CHHN/HHHNO"]:
                    ats_srepr.append("C/CHHN2")
                elif atba2 in ["C/CHHN/CCHHHH", "C/CHHN/CCCHH"]:
                    ats_srepr.append("C/CHHN3")
                elif atba2 in ["C/CHHN/CHHHHH", "C/CHHN/CCHHH"]:
                    ats_srepr.append("C/CHHN4")
                else:
                    raise ValueError(f"{i+1} {atba}")
            elif atba == "C/CHHH":
                if atba2 == "C/CHHH/CHO":
                    ats_srepr.append("C/CHHH1")
                elif atba2 in ["C/CHHH/CHN", "C/CHHH/CHH", "C/CHHH/CCH"]:
                    ats_srepr.append("C/CHHH2")
                else:
                    raise ValueError(f"{i+1} {atba}")
            elif atba == "C/CCN":
                if atba2 == "C/CCN/CCCCHH":
                    ats_srepr.append("C/CCN1")
                elif atba2 in ["C/CCN/CCHHHHN", "C/CCN/CCHHHN"]:
                    ats_srepr.append("C/CCN2")
                else:
                    raise ValueError(f"{i+1} {atba}")
            elif atba == "C/CCHN":
                if atba2.split("/")[2].count("O") == 1:
                    ats_srepr.append("C/CCHN1")
                elif atba2.split("/")[2].count("O") == 2:
                    ats_srepr.append("C/CCHN2")
                elif atba2.split("/")[2].count("O") == 3:
                    ats_srepr.append("C/CCHN3")
                else:
                    raise ValueError(f"{i+1} {atba}")
            elif atba == "C/CCHH":
                if atba2.split("/")[2].count("O") == 2:
                    ats_srepr.append("C/CCHH1")
                elif atba2.split("/")[2].count("O") == 1:
                    ats_srepr.append("C/CCHH2")
                elif atba2.split("/")[2].count("N") == 2:
                    ats_srepr.append("C/CCHH3")
                elif atba2 in ["C/CCHH/CHHHNS", "C/CCHH/CHHHHN", "C/CCHH/CCHHHN", "C/CCHH/CCHHHH", "C/CCHH/CCCHN",
                               "C/CCHH/CCCHHN"]:
                    ats_srepr.append("C/CCHH4")
                else:
                    raise ValueError(f"{i+1} {atba}")
            elif atba == "C/CCH":
                if atba2.split("/")[2].count("O") == 1:
                    ats_srepr.append("C/CCH1")
                elif atba2.split("/")[2].count("N") == 1:
                    ats_srepr.append("C/CCH2")
                elif atba2 == "C/CCH/CCCH":
                    ats_srepr.append("C/CCH3")
                elif atba2 == "C/CCH/CCHH":
                    bonded_atoms = []
                    for at1, at2, _ in bonds:
                        if at1 == i:
                            bonded_atoms.append(at2)
                        if at2 == i:
                            bonded_atoms.append(at1)
                    bonded_atoms = set(bonded_atoms)
                    if bonded_atoms == {"C/CCH/CCCH", "H/C/CC", "C/CCH/CCCHH"}:
                        ats_srepr.append("C/CCH4")
                    else:
                        ats_srepr.append("C/CCH3")
                else:
                    raise ValueError(f"{i+1} {atba}")
            else:
                ats_srepr.append(atba)
        self.ats_srepr = List(ats_srepr)
        self.bonds_srepr = List(bonds_srepr)


    def create_ba(self) -> list:
        bonded_ats = [[] for _ in range(self.n_ats)]
        for bonded_at1, bonded_at2, _ in self.bonds:
            bonded_ats[bonded_at1].append(self.symbols[bonded_at2])
            bonded_ats[bonded_at2].append(self.symbols[bonded_at1])
        return [f"{symbol}/{''.join(sorted(bonded_ats))}"
                for symbol, bonded_ats in zip(self.symbols, bonded_ats)]

    def create_ba2(self) -> list:
        bonded_ats = [[] for _ in range(self.n_ats)]
        for bonded_at1, bonded_at2, _ in self.bonds:
            bonded_ats[bonded_at1].append(self.symbols[bonded_at2])
            bonded_ats[bonded_at2].append(self.symbols[bonded_at1])

        bonded_bonded_ats = [[] for _ in range(self.n_ats)]
        for bonded_at1, bonded_at2, _ in self.bonds:
            bonded_bonded_ats[bonded_at1].extend(bonded_ats[bonded_at2])
            bonded_bonded_ats[bonded_at1].remove(self.symbols[bonded_at1])
            bonded_bonded_ats[bonded_at2].extend(bonded_ats[bonded_at1])
            bonded_bonded_ats[bonded_at2].remove(self.symbols[bonded_at2])

        return [f"{symbol}/{''.join(sorted(bonded_ats))}/{''.join(sorted(bonded_bonded_ats))}"
                for symbol, bonded_ats, bonded_bonded_ats in zip(self.symbols, bonded_ats, bonded_bonded_ats)]


    def calculate_surfaces(self, cpu) -> np.array:
        num_pts = 1000
        kdtree = kdtreen(self.coordinates, leaf_size=40)
        with Pool(cpu) as p:
            surfaces = p.starmap(f2, [[index, at, self.coordinates, self.symbols, kdtree, num_pts] for index, at in enumerate(self.symbols)])
        self.surfaces = np.array(surfaces)


    def create_submolecules(self):
        self.residues = []
        number = 1
        start_index = 0
        residues_numbers = []
        residues_names = []
        for i, at in enumerate(self.rdkit_mol.GetAtoms()):
            if at.GetPDBResidueInfo().GetResidueNumber() == number + 1:
                self.residues.append(Residue(residues_names[start_index],
                                             number-1,
                                             self.coordinates[start_index: i],
                                             self.surfaces[start_index: i],
                                             self.propka_charges[start_index: i],
                                             [x for x in range(start_index, i)],
                                             self.precalc_params[start_index: i]))
                start_index = i
                number += 1
            residues_names.append(at.GetPDBResidueInfo().GetResidueName())
            residues_numbers.append(number)
        self.residues.append(Residue(residues_names[start_index],
                                     number-1,
                                     self.coordinates[start_index: ],
                                     self.surfaces[start_index: ],
                                     self.propka_charges[start_index: ],
                                     [x for x in range(start_index, self.n_ats)],
                                     self.precalc_params[start_index: ]))


        from collections import defaultdict

        residue_bonds = defaultdict(list)
        residue_bonds_srepr = defaultdict(list)
        inter_residues_bonds = defaultdict(list)
        inter_residues_bonds_indices = defaultdict(list)
        inter_residues_bonds_srepr = defaultdict(list)

        for bond, bond_s_repr in zip(self.bonds, self.precalc_bond_hardnesses):
            re0 = residues_numbers[bond[0]]-1
            re1 = residues_numbers[bond[1]]-1
            if re0  == re1:
                residue_bonds[re0].append(bond)
                residue_bonds_srepr[re0].append(bond_s_repr)

            elif re0 < re1:
                inter_residues_bonds_indices[re0].append(re1)
                inter_residues_bonds[re0].append(bond)
                inter_residues_bonds_srepr[re0].append(bond_s_repr)

            elif re0 > re1:
                inter_residues_bonds_indices[re1].append(re0)
                inter_residues_bonds[re1].append(bond)
                inter_residues_bonds_srepr[re1].append(bond_s_repr)

        for residuum in self.residues:
            residue_number = residuum.number
            residuum.bonds = residue_bonds[residue_number]
            residuum.precalc_bond_hardnesses = residue_bonds_srepr[residue_number]
            residuum.connected_with = inter_residues_bonds_indices[residue_number]
            residuum.connected_bonds =  inter_residues_bonds[residue_number]
            residuum.connected_precalc_bond_hardnesses = inter_residues_bonds_srepr[residue_number]



        # residues created


        residues_averages = [res.coordinates_mean for res in self.residues]
        res_kdtree = kdtreen(residues_averages, leaf_size=50)

        self.substructures = []
        for res in self.residues:
            residues = [res]

            distances, indices = res_kdtree.query([res.coordinates_mean], k=len(residues_averages))
            distances = distances[0][1:]
            indices = indices[0][1:]
            for d,i in zip(distances, indices):


                amk_radius = {'ALA': 2.48013472197102,
                             'ARG': 4.861893836930157,
                             'ASN': 3.223781749594369,
                             'ASP': 2.803611164950305,
                             'CYS': 2.5439900881094437,
                             'GLN': 3.845641228833085,
                             'GLU': 3.396388805414707,
                             'GLY': 2.145581026362788,
                             'HIS': 3.837607343643752,
                             'ILE': 3.4050022674834866,
                             'LEU': 3.5357084005904222,
                             'LYS': 4.452109446576905,
                             'MET': 4.18214798399969,
                             'PHE': 4.117010781374703,
                             'PRO': 2.8418414713774762,
                             'SER': 2.499710830364658,
                             'THR': 2.74875502488962,
                             'TRP': 4.683613811874498,
                             'TYR': 4.514843482397014,
                             'VAL': 2.951599144669813}

                if d < amk_radius[res.name] + amk_radius[self.residues[i].name] + 5:
                    residues.append(self.residues[i])
            self.substructures.append(Substructure(residues))



@jit(nopython=True, cache=True)
def dist(grid, c, d):
    indices_to_remove = []
    e, f, g = c
    for gi, (a, b, c) in enumerate(grid):
        if np.sqrt((a - e) ** 2 + (b - f) ** 2 + (c - g) ** 2) < d:
            indices_to_remove.append(gi)
    return indices_to_remove


@jit(nopython=True, cache=True)
def fibonacci_sphere(xc, yc, zc, radius, num_pts):
    indices = np.arange(0, num_pts, dtype=np.float32) + 0.5
    phi = np.arccos(1 - 2 * indices / num_pts)
    theta = np.pi * (1 + 5 ** 0.5) * indices
    x, y, z = np.cos(theta) * np.sin(phi) * radius, np.sin(theta) * np.sin(phi) * radius, np.cos(phi) * radius
    return np.column_stack((x + xc, y + yc, z + zc))



def f2(index, at, coordinates, symbols, kdtree, num_pts):
    vdw = {"H": 1.17183574,
           "C": 1.74471729,
           "N": 1.59013069,
           "O": 1.46105651,
           "S": 1.84999765}
    grid = fibonacci_sphere(*coordinates[index], vdw[at[0]], num_pts)
    distances, indices = kdtree.query([coordinates[index]], k=30)
    distances = distances[0][1:]
    indices = indices[0][1:]
    atom_radius = vdw[symbols[index]]
    for i, index_near in enumerate(indices):
        near_atom_radius = vdw[symbols[index_near]]
        if distances[i] < atom_radius + near_atom_radius:
            grid = np.delete(grid, dist(grid, coordinates[index_near], near_atom_radius), 0)
    return len(grid)/num_pts
