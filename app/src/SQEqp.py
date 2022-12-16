from math import erf
import os
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

import numpy as np
from numba import jit





def calculate_charges(molecule):
    charges = sqeqp_calculate(molecule.bonds,
                              molecule.precalc_bond_hardnesses,
                              molecule.coordinates,
                              molecule.total_chg,
                              molecule.precalc_params)
    return charges[:molecule.calculated_atoms]



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

    T = np.zeros((num_of_bonds, num_of_ats), dtype=np.float64)
    matrix = np.empty((num_of_ats, num_of_ats), dtype=np.float64)

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

