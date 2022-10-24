import json
from math import erf

import numpy as np
from numba import jit
from numba.core import types
from numba.typed import Dict


class SQEqp_h:
    def __init__(self, parameters):
        params = json.load(open(parameters))

        self.parameters = Dict.empty(key_type=types.unicode_type,
                                     value_type=types.float32[:])
        for key,value in params["atom"]["data"].items():
            self.parameters[key] = np.array(value, dtype=np.float32)

        self.widths = Dict.empty(key_type=types.unicode_type,
                                 value_type=types.float32[:])
        width_index = params["atom"]["names"].index("width")
        width_cor_index = params["atom"]["names"].index("width_cor")
        for key, value in params["atom"]["data"].items():
            self.widths[key] = np.array((value[width_index], value[width_cor_index]), dtype=np.float32)

        bond_hardnesses = Dict.empty(key_type=types.unicode_type,
                                   value_type=types.float64)
        for key, value in params["bond"]["data"].items():
            bond_hardnesses[key] = value
        self.bond_hardnesses = bond_hardnesses


def calculate_charges(molecule):
    charges = sqeqp_calculate(molecule.bonds,
                              molecule.precalc_bond_hardnesses,
                              molecule.coordinates,
                              molecule.total_chg,
                              molecule.precalc_params)
    return charges[:molecule.valid_atoms]



@jit(nopython=True, cache=True)
def sqeqp_calculate(bonds,
                    precalc_bond_hardnesses,
                    coordinates,
                    total_chg,
                    precalc_params):

    electronegativities = precalc_params[:, 0]
    hardnesses = precalc_params[:, 1]
    radiuses = precalc_params[:, 2]
    initial_charges = precalc_params[:, 3]


    num_of_ats = len(coordinates)
    num_of_bonds = len(bonds)

    T = np.zeros((num_of_bonds, num_of_ats), dtype=np.float32)
    matrix = np.empty((num_of_ats, num_of_ats), dtype=np.float32)

    for i,(a1, a2, _) in enumerate(bonds):
        T[i, a1] += 1
        T[i, a2] -= 1

    for i in range(num_of_ats):
        matrix[i, i] = hardnesses[i]
        i_radius = radiuses[i]
        ix, iy, iz = coordinates[i]

        for j, (j_radius, (jx,jy,jz)) in enumerate(zip(radiuses[i+1:],
                                                       coordinates[i+1:]),
                                                       i + 1):
            d0 = np.sqrt(i_radius + j_radius)
            distance = np.sqrt((ix - jx) ** 2 + (iy - jy) ** 2 + (iz - jz) ** 2)
            matrix[i, j] = matrix[j, i] = erf(distance / d0) / distance

    initial_charges -= (np.sum(initial_charges) - total_chg) / len(initial_charges)


    A_sqe = np.dot(T, np.dot(matrix, T.T))
    for i, hardness in enumerate(precalc_bond_hardnesses):
        A_sqe[i, i] += hardness


    electronegativities -= np.dot(matrix, initial_charges)
    electronegativities += hardnesses * initial_charges
    B_sqe = np.dot(T, electronegativities)



    r = np.dot(np.linalg.solve(A_sqe, B_sqe), T) + initial_charges
    return r

