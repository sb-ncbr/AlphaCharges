import numpy as np
from numba.typed import List
from rdkit import Chem
from sklearn.neighbors import KDTree as kdtreen

from .amino_acids_atomic_types import real_ats_types
from .amino_acids_radii import amk_radius
from .mean_qm_charges import mean_qm_charges
from .problematic_atom_info import problematic_atom_info


class Substructure:
    def __init__(self,
                 residues: list):
        self.residues = residues
        self.calculated_atoms = len(residues[0].coordinates)
        self.coordinates = np.concatenate([res.coordinates for res in residues], axis=0)
        self.precalc_params = np.concatenate([res.precalc_params for res in residues], axis=0)
        self.total_chg = sum(chg for res in residues for chg in res.mean_qm_chgs)
        precalc_bond_hardnesses = [hardness for res in residues for hardness in res.precalc_bond_hardnesses]
        bonds = [bond for res in residues for bond in res.bonds]

        # add bonds between residues (peptide and disulfide bonds)
        residues_numbers = set(res.number for res in residues)
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
                 mean_qm_chgs: np.array,
                 indices: np.array,
                 precalc_params: np.array):

        self.name = name
        self.number = number
        self.coordinates = coordinates
        self.coordinates_mean = np.mean(self.coordinates, axis=0)
        self.mean_qm_chgs = mean_qm_chgs
        self.indices = indices
        self.precalc_params = precalc_params


class Molecule:
    def __init__(self,
                 pdb_file: str,
                 pqr_file: str):

        # load charges from propka, sum of charges are used as total charge of molecule
        self.total_chg = round(sum(float(line[55:62]) for line in open(pqr_file, "r").readlines()[:-2]))

        # load molecule by rdkit
        self.rdkit_mol = Chem.MolFromPDBFile(pdb_file,
                                             removeHs=False,
                                             sanitize=False)

        # load atoms and bonds
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
        bonds_srepr = [f"{'-'.join(sorted([ats_sreprba[ba1], ats_sreprba[ba2]]))}-{bond_type}"
                       for ba1, ba2, bond_type in bonds]

        # control, whether molecule consist of standart aminoacids
        problematic_atoms = {}
        for i, (atba, rdkit_at) in enumerate(zip(ats_sreprba,
                                                 self.rdkit_mol.GetAtoms())):
            if atba not in real_ats_types:
                res_info = rdkit_at.GetPDBResidueInfo()
                label_comp_id = f"{res_info.GetResidueName()}"
                label_seq_id = int(res_info.GetResidueNumber())
                label_atom_id = f"{res_info.GetName().strip()}"
                id = f"{label_comp_id} {label_seq_id} {label_atom_id}"
                message = problematic_atom_info(rdkit_at, atba, self.rdkit_mol)
                problematic_atoms[id] = {
                    "key": {
                        "labelCompId": label_comp_id,
                        "labelSeqId": label_seq_id,
                        "labelAtomId": label_atom_id,
                    },
                    "message": message,
                }
        if problematic_atoms:
            raise ValueError(problematic_atoms)

        self.mean_qm_chgs = [mean_qm_charges[ats_srepr] for ats_srepr in ats_sreprba]

        # convert to numba data structure
        self.ats_srepr = List(ats_sreprba)
        self.bonds_srepr = List(bonds_srepr)

    def create_ba(self) -> list:
        bonded_ats = [[] for _ in range(self.n_ats)]
        for bonded_at1, bonded_at2, _ in self.bonds:
            bonded_ats[bonded_at1].append(self.symbols[bonded_at2])
            bonded_ats[bonded_at2].append(self.symbols[bonded_at1])
        return [f"{symbol}/{''.join(sorted(bonded_ats))}"
                for symbol, bonded_ats in zip(self.symbols, bonded_ats)]

    def create_submolecules(self):
        # create residues
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
                                             self.mean_qm_chgs[start_index: i],
                                             [x for x in range(start_index, i)],
                                             self.precalc_params[start_index: i]))
                start_index = i
                number += 1
            residues_names.append(at.GetPDBResidueInfo().GetResidueName())
            residues_numbers.append(number)
        self.residues.append(Residue(residues_names[start_index],
                                     number-1,
                                     self.coordinates[start_index: ],
                                     self.mean_qm_chgs[start_index: ],
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

        # create submolecules
        residues_averages = [res.coordinates_mean for res in self.residues]
        res_kdtree = kdtreen(residues_averages, leaf_size=50)
        self.substructures = []
        for res in self.residues:
            residues = [res]
            distances, indices = res_kdtree.query([res.coordinates_mean], k=len(residues_averages))
            distances = distances[0][1:]
            indices = indices[0][1:]
            for d,i in zip(distances, indices):
                if d < amk_radius[res.name] + amk_radius[self.residues[i].name] + 5:
                    residues.append(self.residues[i])
            self.substructures.append(Substructure(residues))
