import json
import numba
from math import erf
from numba.core import types
from numba.typed import Dict
import numpy as np

class SQEqp_h:
    def __init__(self, parameters):
        params = json.load(open(parameters))
        for i, parameter_name in enumerate(["electronegativity",
                                            "hardness",
                                            "width",
                                            "q0",
                                            "q0_cor",
                                            "hardness_cor",
                                            "electronegativity_cor",
                                            "width_cor"]):
             parameters = Dict.empty(key_type=types.unicode_type,
                                     value_type=types.float64)
             for key, value in params["atom"]["data"].items():
                 parameters[key] = value[i]
             setattr(self, parameter_name, parameters)
        bond_hardness = Dict.empty(key_type=types.unicode_type,
                                   value_type=types.float64)
        for key, value in params["bond"]["data"].items():
            bond_hardness[key] = value
        self.bond_hardness = bond_hardness

    def calculate_charges(self, molecule):
        charges = sqeqp_calculate(molecule.ats_srepr,
                        molecule.bonds,
                        molecule.bonds_srepr,
                        molecule.distance_matrix,
                        molecule.total_chg,
                        molecule.surfaces,

                    self.electronegativity,
                                            self.hardness,
                                            self.width,
                                            self.q0,
                                            self.q0_cor,
                                            self.hardness_cor,
                                            self.electronegativity_cor,
                                            self.width_cor, self.bond_hardness)

        return charges


@numba.jit(nopython=True)
def sqeqp_calculate(ats_srepr, bonds, bonds_srepr, distance_matrix, total_chg, surfaces,

                    pelectronegativity,
                                            phardness,
                                            pwidth,
                                            pq0,
                                            pq0_cor,
                                            phardness_cor,
                                            pelectronegativity_cor,
                                            pwidth_cor, bonds_hardness):
    num_of_at = len(ats_srepr)
    num_of_bonds = len(bonds)
    T = np.zeros((num_of_bonds, num_of_at))
    for i in range(len(bonds)):
        at1, at2, _ = bonds[i]
        T[i, at1] += 1
        T[i, at2] -= 1
    matrix = np.zeros((num_of_at, num_of_at))
    vector = np.zeros(num_of_at)
    list_of_q0 = np.empty(num_of_at, dtype=np.float64)
    list_of_hardness = np.empty(num_of_at, dtype=np.float64)
    for i, symbol_i in enumerate(ats_srepr):
        corected_hardness = phardness[symbol_i] + phardness_cor[symbol_i] * surfaces[i]
        matrix[i, i] = corected_hardness
        list_of_hardness[i] = corected_hardness

        vector[i] = -pelectronegativity[symbol_i] + pelectronegativity_cor[symbol_i] * surfaces[i]
        list_of_q0[i] = pq0[symbol_i] + pq0_cor[symbol_i] * surfaces[i]

        for j, (symbol_j, distance) in enumerate(zip(ats_srepr[i + 1:],
                                                        distance_matrix[i, i + 1:]),
                                                    i + 1):
            d0 = np.sqrt(2 * (pwidth[symbol_i] + pwidth_cor[symbol_i] * surfaces[i]) ** 2 + 2 * (pwidth[symbol_j]+ pwidth_cor[symbol_j] * surfaces[j]) ** 2)
            matrix[i, j] = matrix[j, i] = erf(distance / d0) / distance
    list_of_q0 -= (np.sum(list_of_q0) - total_chg) / len(list_of_q0)
    vector -= np.dot(matrix, list_of_q0)
    vector += list_of_hardness * list_of_q0
    A_sqe = np.dot(T, np.dot(matrix, T.T))
    B_sqe = np.dot(T, vector)
    for i, bond_srepr in enumerate(bonds_srepr):
        A_sqe[i, i] += bonds_hardness[bond_srepr]
    return np.dot(np.linalg.solve(A_sqe, B_sqe), T) + list_of_q0
