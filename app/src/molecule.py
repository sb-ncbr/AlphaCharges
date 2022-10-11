from rdkit import Chem
import numpy as np
from scipy.spatial.distance import cdist
from scipy.spatial import KDTree
from numba import jit
from numba.typed import List

from multiprocessing import Pool


class Molecule:
    def __init__(self,
                 code: str,
                 pdb_file: str):

        self.cpu = 4
        self.code = code

        pqr_file_lines = open(pdb_file.replace(".pdb", ".pqr"), "r").readlines()[:-2]
        self.total_chg = round(sum([float(line.split()[8]) for line in pqr_file_lines]))

        rdkit_mol = Chem.MolFromPDBFile(pdb_file, removeHs=False, sanitize = True)
        self.symbols = [atom.GetSymbol() for atom in rdkit_mol.GetAtoms()]
        print(f"Number of atoms: {len(self.symbols)}")
        self.n_ats = len(self.symbols)

        bonds = []
        bond_types = {"SINGLE": 1,
                      "DOUBLE": 2,
                      "TRIPLE": 3,
                      "AROMATIC": 4}
        for bond in rdkit_mol.GetBonds():
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
        for i in range(0, rdkit_mol.GetNumAtoms()):
            pos = rdkit_mol.GetConformer().GetAtomPosition(i)
            coordinates.append((pos.x, pos.y, pos.z))
        self.coordinates = np.array(coordinates)

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

        # baex
        ats_srepr = []
        for i, (atba, atba2) in enumerate(zip(ats_sreprba, ats_sreprba2)):
            if atba in ['H/N', "H/C", "O/C", "N/CHH", "N/CHHH", "C/CHHS", "C/CCC"]:
                ats_srepr.append(atba2)
            elif atba == "O/CH":
                if atba2 == "O/CH/CC":
                    ats_srepr.append("O/CH1")
                elif atba2 in ["O/CH/CHH", "O/CH/CCH"]:
                    ats_srepr.append("O/CH2")
                else:
                    exit("FATAL ERROR1")
            elif atba == "N/CCH":
                if len(atba2.split("/")[2]) == 4:
                    ats_srepr.append("N/CCH1")
                elif len(atba2.split("/")[2]) == 5:
                    ats_srepr.append("N/CCH2")
                elif len(atba2.split("/")[2]) == 3:
                    ats_srepr.append("N/CCH3")
                else:
                    exit("FATAL ERROR2")
            elif atba == "C/CHN":
                if atba2 == "C/CHN/CCCH":
                    ats_srepr.append("C/CHN1")
                elif atba2 in ["C/CHN/CCN", "C/CHN/CCHN"]:
                    ats_srepr.append("C/CHN2")
                else:
                    exit("FATAL ERROR3")
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
                    exit("FATAL ERROR4")
            elif atba == "C/CHHH":
                if atba2 == "C/CHHH/CHO":
                    ats_srepr.append("C/CHHH1")
                elif atba2 in ["C/CHHH/CHN", "C/CHHH/CHH", "C/CHHH/CCH"]:
                    ats_srepr.append("C/CHHH2")
                else:
                    exit("FATAL ERROR10")
            elif atba == "C/CCN":
                if atba2 == "C/CCN/CCCCHH":
                    ats_srepr.append("C/CCN1")
                elif atba2 in ["C/CCN/CCHHHHN", "C/CCN/CCHHHN"]:
                    ats_srepr.append("C/CCN2")
                else:
                    exit("FATAL ERROR5")
            elif atba == "C/CCHN":
                if atba2.split("/")[2].count("O") == 1:
                    ats_srepr.append("C/CCHN1")
                elif atba2.split("/")[2].count("O") == 2:
                    ats_srepr.append("C/CCHN2")
                elif atba2.split("/")[2].count("O") == 3:
                    ats_srepr.append("C/CCHN3")
                else:
                    exit("FATAL ERROR6")
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
                    exit("FATAL ERROR7")
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
                    exit("FATAL ERROR8")
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

    def calculate_distace_matrix(self):
        self.distance_matrix = cdist(self.coordinates, self.coordinates)

        # from time import time
        #
        #
        #
        # s = time()
        # distance_matrix = np.empty((self.n_ats, self.n_ats), dtype=np.float32)
        # with Pool(self.cpu) as p:
        #     p.starmap(calculate_distance_matrix_partial, [[distance_matrix, self.coordinates, index] for index in range(self.n_ats)])
        # print(time()-s)
        # exit()
        #




    def calculate_surfaces(self) -> np.array:
        num_pts = 1000
        kdtree = KDTree(self.coordinates, leafsize=50)
        with Pool(self.cpu) as p:
            surfaces = p.starmap(f, [[index, at, self.coordinates, self.symbols, kdtree, num_pts] for index, at in enumerate(self.symbols)])
        self.surfaces = np.array(surfaces)


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
    indices = np.arange(0, num_pts, dtype=np.float64) + 0.5
    phi = np.arccos(1 - 2 * indices / num_pts)
    theta = np.pi * (1 + 5 ** 0.5) * indices
    x, y, z = np.cos(theta) * np.sin(phi) * radius, np.sin(theta) * np.sin(phi) * radius, np.cos(phi) * radius
    return np.column_stack((x + xc, y + yc, z + zc))

def f(index, at, coordinates, symbols, kdtree, num_pts):
    vdw = {"H": 1.17183574,
           "C": 1.74471729,
           "N": 1.59013069,
           "O": 1.46105651,
           "S": 1.84999765}
    grid = fibonacci_sphere(*coordinates[index], vdw[at[0]], num_pts)
    near_atoms = kdtree.query_ball_point(coordinates[index], vdw["S"] + vdw[at])
    near_atoms.remove(index)
    for index_near in near_atoms:
        grid = np.delete(grid, dist(grid, coordinates[index_near], vdw[symbols[index_near][0]]), 0)
    return len(grid)/num_pts


def calculate_distance_matrix_partial(distance_matrix, coordinates, index):
    r = cdist([coordinates[index]], coordinates[index:])
    distance_matrix[index, index:] = r
    distance_matrix[index:, index] = r
