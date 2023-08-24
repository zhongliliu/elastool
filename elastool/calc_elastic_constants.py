"""
  Elastool -- Elastic toolkit for zero and finite-temperature elastic constants and mechanical properties calculations

  Copyright (C) 2019-2024 by Zhong-Li Liu and Chinedu Ekuma

  This program is free software; you can redistribute it and/or modify it under the
  terms of the GNU General Public License as published by the Free Software Foundation
  version 3 of the License.

  This program is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
  PARTICULAR PURPOSE.  See the GNU General Public License for more details.

  E-mail: zl.liu@163.com, cekuma1@gmail.com
"""

from copy import deepcopy
from numpy import array, vstack, linalg, pi
from strain_matrix import strain_matrix
from read_input import indict


def Cubic(latt_system, elastic_constants_dict, stress_set_dict, conversion_factor):
    for n, up in enumerate(stress_set_dict.keys()):
        s = strain_matrix(latt_system, up)[0]

        #print(s)
        stress_list = array(stress_set_dict[up][0])

        eplisons_now = array([[s[0][0], s[1][1]+s[2][2], 0.],
                              [s[1][1], s[0][0]+s[2][2], 0.],
                              [s[2][2], s[0][0]+s[1][1], 0.],
                              [0.,      0.,        2*s[1][2]],
                              [0.,      0.,        2*s[0][2]],
                              [0.,      0.,        2*s[0][1]]])

        stresses_now = array([[stress_list[0]],
                              [stress_list[1]],
                              [stress_list[2]],
                              [stress_list[3]],
                              [stress_list[4]],
                              [stress_list[5]]])
        if n == 0:
            eplisons = deepcopy(eplisons_now)
            stresses = deepcopy(stresses_now)
        else:
            eplisons = vstack((eplisons, eplisons_now))
            stresses = vstack((stresses, stresses_now))

    cij = linalg.lstsq(eplisons, stresses, rcond=None)[0] * conversion_factor

    c11 = cij[0][0]
    c12 = cij[1][0]
    c44 = cij[2][0]

    B_v = (c11+2*c12)/3.
    B_r = B_v

    G_v = (c11-c12+3*c44)/5.
    G_r = 5*(c11-c12)*c44/(4*c44+3*(c11-c12))

    B_vrh = (B_v+B_r)/2.
    G_vrh = (G_v+G_r)/2.
    E = 9*B_vrh*G_vrh/(3*B_vrh+G_vrh)
    v = (3*B_vrh-2*G_vrh)/(2*(3*B_vrh+G_vrh))

    elastic_constants_dict['c11'] = c11
    elastic_constants_dict['c12'] = c12
    elastic_constants_dict['c44'] = c44

    elastic_constants_dict['B_v'] = B_v
    elastic_constants_dict['B_r'] = B_r
    elastic_constants_dict['G_v'] = G_v
    elastic_constants_dict['G_r'] = G_r
    elastic_constants_dict['B_vrh'] = B_vrh
    elastic_constants_dict['G_vrh'] = G_vrh
    elastic_constants_dict['E'] = E
    elastic_constants_dict['v'] = v

    return elastic_constants_dict


def Hexagonal(latt_system, elastic_constants_dict, stress_set_dict, conversion_factor):
    for n, up in enumerate(stress_set_dict.keys()):
        if up == 0:
            s1 = strain_matrix(latt_system, up)[0]
            stress_list1 = array(stress_set_dict[up][0])

            eplisons_now = array([[s1[0][0],  s1[1][1],  s1[2][2],  0., 0.],
                                  [s1[1][1],  s1[0][0],  s1[2][2],  0., 0.],
                                  [0., 0., s1[0][0]+s1[1][1], s1[2][2], 0.],
                                  [0.,    0.,     0.,     0.,   2*s1[1][2]],
                                  [0.,    0.,     0.,     0.,   2*s1[0][2]],
                                  [s1[0][1],   -s1[0][1], 0.,       0., 0.]])

            stresses_now = array([[stress_list1[0]],
                                  [stress_list1[1]],
                                  [stress_list1[2]],
                                  [stress_list1[3]],
                                  [stress_list1[4]],
                                  [stress_list1[5]]])
        else:
            s1 = strain_matrix(latt_system, up)[0]
            s2 = strain_matrix(latt_system, up)[1]
            stress_list1 = array(stress_set_dict[up][0])
            stress_list2 = array(stress_set_dict[up][1])

            eplisons_now = array([[s1[0][0],  s1[1][1],  s1[2][2],  0., 0.],
                                  [s1[1][1],  s1[0][0],  s1[2][2],  0., 0.],
                                  [0., 0., s1[0][0]+s1[1][1], s1[2][2], 0.],
                                  [0.,    0.,     0.,     0.,   2*s1[1][2]],
                                  [0.,    0.,     0.,     0.,   2*s1[0][2]],
                                  [s1[0][1],   -s1[0][1], 0.,       0., 0.],
                                  [s2[0][0],  s2[1][1],  s2[2][2],  0., 0.],
                                  [s2[1][1],  s2[0][0],  s2[2][2],  0., 0.],
                                  [0., 0., s2[0][0]+s2[1][1], s2[2][2], 0.],
                                  [0.,    0.,     0.,     0.,   2*s2[1][2]],
                                  [0.,    0.,     0.,     0.,   2*s2[0][2]],
                                  [s2[0][1],   -s2[0][1], 0.,       0., 0.]])

            stresses_now = array([[stress_list1[0]],
                                  [stress_list1[1]],
                                  [stress_list1[2]],
                                  [stress_list1[3]],
                                  [stress_list1[4]],
                                  [stress_list1[5]],
                                  [stress_list2[0]],
                                  [stress_list2[1]],
                                  [stress_list2[2]],
                                  [stress_list2[3]],
                                  [stress_list2[4]],
                                  [stress_list2[5]]])
        if n == 0:
            eplisons = deepcopy(eplisons_now)
            stresses = deepcopy(stresses_now)
        else:
            eplisons = vstack((eplisons, eplisons_now))
            stresses = vstack((stresses, stresses_now))

    cij = linalg.lstsq(eplisons, stresses, rcond=None)[0] * conversion_factor

    c11 = cij[0][0]
    c12 = cij[1][0]
    c13 = cij[2][0]
    c33 = cij[3][0]
    c44 = cij[4][0]

    M = c11+c12+2*c33-4*c13
    C2 = (c11+c12)*c33-2*c13*c13
    c66 = (c11-c12)/2.

    B_v = (2*(c11+c12)+4*c13+c33)/9.
    G_v = (M+12*c44+12*c66)/30.
    B_r = C2/M
    G_r = 2.5*(C2*c44*c66)/(3*B_v*c44*c66+C2*(c44+c66))

    B_vrh = (B_v+B_r)/2.
    G_vrh = (G_v+G_r)/2.
    E = 9*B_vrh*G_vrh/(3*B_vrh+G_vrh)
    v = (3*B_vrh-2*G_vrh)/(2*(3*B_vrh+G_vrh))

    elastic_constants_dict['c11'] = c11
    elastic_constants_dict['c12'] = c12
    elastic_constants_dict['c13'] = c13
    elastic_constants_dict['c33'] = c33
    elastic_constants_dict['c44'] = c44

    elastic_constants_dict['B_v'] = B_v
    elastic_constants_dict['B_r'] = B_r
    elastic_constants_dict['G_v'] = G_v
    elastic_constants_dict['G_r'] = G_r
    elastic_constants_dict['B_vrh'] = B_vrh
    elastic_constants_dict['G_vrh'] = G_vrh
    elastic_constants_dict['E'] = E
    elastic_constants_dict['v'] = v

    return elastic_constants_dict


def Trigonal1(latt_system, elastic_constants_dict, stress_set_dict, conversion_factor):
    for n, up in enumerate(stress_set_dict.keys()):
        if up == 0:
            s1 = strain_matrix(latt_system, up)[0]
            stress_list1 = array(stress_set_dict[up][0])

            eplisons_now = array([[s1[0][0],  s1[1][1],  s1[2][2], 2*s1[1][2],  0., 0.],
                                  [s1[1][1],  s1[0][0],  s1[2][2],-2*s1[1][2],  0., 0.],
                                  [0.,    0., s1[0][0]+s1[1][1],    0.,   s1[2][2], 0.],
                                  [0.,    0.,  0.,  s1[0][0]-s1[1][1],  0., 2*s1[1][2]],
                                  [0.,    0.,     0., 2*s1[0][1],       0., 2*s1[0][2]],
                                  [s1[0][1],  -s1[0][1],  0.,   2*s1[0][2],   0.,   0.]])

            stresses_now = array([[stress_list1[0]],
                                  [stress_list1[1]],
                                  [stress_list1[2]],
                                  [stress_list1[3]],
                                  [stress_list1[4]],
                                  [stress_list1[5]]])
        else:
            s1 = strain_matrix(latt_system, up)[0]
            s2 = strain_matrix(latt_system, up)[1]
            stress_list1 = array(stress_set_dict[up][0])
            stress_list2 = array(stress_set_dict[up][1])

            eplisons_now = array([[s1[0][0],  s1[1][1],  s1[2][2], 2*s1[1][2],  0., 0.],
                                  [s1[1][1],  s1[0][0],  s1[2][2],-2*s1[1][2],  0., 0.],
                                  [0.,    0., s1[0][0]+s1[1][1],    0.,   s1[2][2], 0.],
                                  [0.,    0.,  0.,  s1[0][0]-s1[1][1],  0., 2*s1[1][2]],
                                  [0.,    0.,     0., 2*s1[0][1],       0., 2*s1[0][2]],
                                  [s1[0][1],  -s1[0][1],  0.,   2*s1[0][2],   0.,   0.],
                                  [s2[0][0],  s2[1][1],  s2[2][2], 2*s2[1][2],  0., 0.],
                                  [s2[1][1],  s2[0][0],  s2[2][2],-2*s2[1][2],  0., 0.],
                                  [0.,    0., s2[0][0]+s2[1][1],    0.,   s2[2][2], 0.],
                                  [0.,    0.,  0.,  s2[0][0]-s2[1][1],  0., 2*s2[1][2]],
                                  [0.,    0.,     0., 2*s2[0][1],       0., 2*s2[0][2]],
                                  [s2[0][1],  -s2[0][1],  0.,   2*s2[0][2],   0.,   0.]])

            stresses_now = array([[stress_list1[0]],
                                  [stress_list1[1]],
                                  [stress_list1[2]],
                                  [stress_list1[3]],
                                  [stress_list1[4]],
                                  [stress_list1[5]],
                                  [stress_list2[0]],
                                  [stress_list2[1]],
                                  [stress_list2[2]],
                                  [stress_list2[3]],
                                  [stress_list2[4]],
                                  [stress_list2[5]]])
        if n == 0:
            eplisons = deepcopy(eplisons_now)
            stresses = deepcopy(stresses_now)
        else:
            eplisons = vstack((eplisons, eplisons_now))
            stresses = vstack((stresses, stresses_now))

    cij = linalg.lstsq(eplisons, stresses, rcond=None)[0] * conversion_factor

    c11 = cij[0][0]
    c12 = cij[1][0]
    c13 = cij[2][0]
    c14 = cij[3][0]
    c33 = cij[4][0]
    c44 = cij[5][0]
    elastic_constants_dict['c11'] = c11
    elastic_constants_dict['c12'] = c12
    elastic_constants_dict['c13'] = c13
    elastic_constants_dict['c14'] = c14
    elastic_constants_dict['c33'] = c33
    elastic_constants_dict['c44'] = c44

    return elastic_constants_dict


def Trigonal2(latt_system, elastic_constants_dict, stress_set_dict, conversion_factor):
    for n, up in enumerate(stress_set_dict.keys()):
        if up == 0:
            s1 = strain_matrix(latt_system, up)[0]
            stress_list1 = array(stress_set_dict[up][0])

            eplisons_now = array([[s1[0][0],  s1[1][1],  s1[2][2],  2*s1[1][2],  2*s1[0][2],  0., 0.],
                                 [s1[1][1],  s1[0][0],  s1[2][2], -2*s1[1][2], -2*s1[0][2],  0., 0.],
                                 [0.,    0., s1[0][0]+s1[1][1],      0.,      0.,      s1[2][2], 0.],
                                 [0.,    0.,     0., s1[0][0]-s1[1][1], -2*s1[0][1], 0., 2*s1[1][2]],
                                 [0.,    0.,     0., 2*s1[0][1], s1[0][0]-s1[1][1],  0., 2*s1[0][2]],
                                 [s1[0][1],  -s1[0][1],  0.,   2*s1[0][2],   -2*s1[1][2],   0.,  0.]])

            stresses_now = array([[stress_list1[0]],
                                  [stress_list1[1]],
                                  [stress_list1[2]],
                                  [stress_list1[3]],
                                  [stress_list1[4]],
                                  [stress_list1[5]]])
        else:
            s1 = strain_matrix(latt_system, up)[0]
            s2 = strain_matrix(latt_system, up)[1]
            stress_list1 = array(stress_set_dict[up][0])
            stress_list2 = array(stress_set_dict[up][1])

            eplisons_now = array([[s1[0][0],  s1[1][1],  s1[2][2],  2*s1[1][2],  2*s1[0][2],  0., 0.],
                                  [s1[1][1],  s1[0][0],  s1[2][2], -2*s1[1][2], -2*s1[0][2],  0., 0.],
                                  [0.,    0., s1[0][0]+s1[1][1],      0.,      0.,      s1[2][2], 0.],
                                  [0.,    0.,     0., s1[0][0]-s1[1][1], -2*s1[0][1], 0., 2*s1[1][2]],
                                  [0.,    0.,     0., 2*s1[0][1], s1[0][0]-s1[1][1],  0., 2*s1[0][2]],
                                  [s1[0][1],  -s1[0][1],  0.,   2*s1[0][2],   -2*s1[1][2],   0.,  0.],
                                  [s2[0][0],  s2[1][1],  s2[2][2],  2*s2[1][2],  2*s2[0][2],  0., 0.],
                                  [s2[1][1],  s2[0][0],  s2[2][2], -2*s2[1][2], -2*s2[0][2],  0., 0.],
                                  [0.,    0., s2[0][0]+s2[1][1],      0.,      0.,      s2[2][2], 0.],
                                  [0.,    0.,     0., s2[0][0]-s2[1][1], -2*s2[0][1], 0., 2*s2[1][2]],
                                  [0.,    0.,     0., 2*s2[0][1], s2[0][0]-s2[1][1],  0., 2*s2[0][2]],
                                  [s2[0][1],  -s2[0][1],  0.,   2*s2[0][2],   -2*s2[1][2],   0.,  0.]])

            stresses_now = array([[stress_list1[0]],
                                  [stress_list1[1]],
                                  [stress_list1[2]],
                                  [stress_list1[3]],
                                  [stress_list1[4]],
                                  [stress_list1[5]],
                                  [stress_list2[0]],
                                  [stress_list2[1]],
                                  [stress_list2[2]],
                                  [stress_list2[3]],
                                  [stress_list2[4]],
                                  [stress_list2[5]]])
        if n == 0:
            eplisons = deepcopy(eplisons_now)
            stresses = deepcopy(stresses_now)
        else:
            eplisons = vstack((eplisons, eplisons_now))
            stresses = vstack((stresses, stresses_now))

    cij = linalg.lstsq(eplisons, stresses, rcond=None)[0] * conversion_factor

    c11 = cij[0][0]
    c12 = cij[1][0]
    c13 = cij[2][0]
    c14 = cij[3][0]
    c15 = cij[4][0]
    c33 = cij[5][0]
    c44 = cij[6][0]
    elastic_constants_dict['c11'] = c11
    elastic_constants_dict['c12'] = c12
    elastic_constants_dict['c13'] = c13
    elastic_constants_dict['c14'] = c14
    elastic_constants_dict['c15'] = c15
    elastic_constants_dict['c33'] = c33
    elastic_constants_dict['c44'] = c44

    return elastic_constants_dict
    

def Tetragonal1(latt_system, elastic_constants_dict, stress_set_dict, conversion_factor):
    for n, up in enumerate(stress_set_dict.keys()):
        if up == 0:
            s1 = strain_matrix(latt_system, up)[0]
            stress_list1 = array(stress_set_dict[up][0])

            eplisons_now = array([[s1[0][0],  s1[1][1],  s1[2][2],   0.,   0.,   0.],
                                  [s1[1][1],  s1[0][0],  s1[2][2],   0.,   0.,   0.],
                                  [0.,    0., s1[0][0]+s1[1][1], s1[2][2], 0.,   0.],
                                  [0.,    0.,     0.,     0.,      2*s1[1][2],   0.],
                                  [0.,    0.,     0.,     0.,      2*s1[0][2],   0.],
                                  [0.,    0.,     0.,     0.,       0.,  2*s1[0][1]]])

            stresses_now = array([[stress_list1[0]],
                                  [stress_list1[1]],
                                  [stress_list1[2]],
                                  [stress_list1[3]],
                                  [stress_list1[4]],
                                  [stress_list1[5]]])
        else:
            s1 = strain_matrix(latt_system, up)[0]
            s2 = strain_matrix(latt_system, up)[1]
            stress_list1 = array(stress_set_dict[up][0])
            stress_list2 = array(stress_set_dict[up][1])

            eplisons_now = array([[s1[0][0],  s1[1][1],  s1[2][2],   0.,   0.,   0.],
                                  [s1[1][1],  s1[0][0],  s1[2][2],   0.,   0.,   0.],
                                  [0.,    0., s1[0][0]+s1[1][1], s1[2][2], 0.,   0.],
                                  [0.,    0.,     0.,     0.,      2*s1[1][2],   0.],
                                  [0.,    0.,     0.,     0.,      2*s1[0][2],   0.],
                                  [0.,    0.,     0.,     0.,       0.,  2*s1[0][1]],
                                  [s2[0][0],  s2[1][1],  s2[2][2],   0.,   0.,   0.],
                                  [s2[1][1],  s2[0][0],  s2[2][2],   0.,   0.,   0.],
                                  [0.,    0., s2[0][0]+s2[1][1], s2[2][2], 0.,   0.],
                                  [0.,    0.,     0.,     0.,      2*s2[1][2],   0.],
                                  [0.,    0.,     0.,     0.,      2*s2[0][2],   0.],
                                  [0.,    0.,     0.,     0.,       0.,  2*s2[0][1]]])

            stresses_now = array([[stress_list1[0]],
                                  [stress_list1[1]],
                                  [stress_list1[2]],
                                  [stress_list1[3]],
                                  [stress_list1[4]],
                                  [stress_list1[5]],
                                  [stress_list2[0]],
                                  [stress_list2[1]],
                                  [stress_list2[2]],
                                  [stress_list2[3]],
                                  [stress_list2[4]],
                                  [stress_list2[5]]])
        if n == 0:
            eplisons = deepcopy(eplisons_now)
            stresses = deepcopy(stresses_now)
        else:
            eplisons = vstack((eplisons, eplisons_now))
            stresses = vstack((stresses, stresses_now))

    cij = linalg.lstsq(eplisons, stresses, rcond=None)[0] * conversion_factor

    c11 = cij[0][0]
    c12 = cij[1][0]
    c13 = cij[2][0]
    c33 = cij[3][0]
    c44 = cij[4][0]
    c66 = cij[5][0]

    M = c11+c12+2*c33-4*c13
    C2 = (c11+c12)*c33-2*c13*c13

    B_v = (2*(c11+c12)+c33+4*c13)/9.
    G_v = (M+3*c11-3*c12+12*c44+6*c66)/30.
    B_r = C2/M
    G_r = 15./(18*B_v/C2+6/(c11-c12)+6/c44+3/c66)

    B_vrh = (B_v+B_r)/2.
    G_vrh = (G_v+G_r)/2.
    E = 9*B_vrh*G_vrh/(3*B_vrh+G_vrh)
    v = (3*B_vrh-2*G_vrh)/(2*(3*B_vrh+G_vrh))

    elastic_constants_dict['c11'] = c11
    elastic_constants_dict['c12'] = c12
    elastic_constants_dict['c13'] = c13
    elastic_constants_dict['c33'] = c33
    elastic_constants_dict['c44'] = c44
    elastic_constants_dict['c66'] = c66

    elastic_constants_dict['B_v'] = B_v
    elastic_constants_dict['B_r'] = B_r
    elastic_constants_dict['G_v'] = G_v
    elastic_constants_dict['G_r'] = G_r

    elastic_constants_dict['B_vrh'] = B_vrh
    elastic_constants_dict['G_vrh'] = G_vrh
    elastic_constants_dict['E'] = E
    elastic_constants_dict['v'] = v

    return elastic_constants_dict


def Tetragonal2(latt_system, elastic_constants_dict, stress_set_dict, conversion_factor):
    for n, up in enumerate(stress_set_dict.keys()):
        if up == 0:
            s1 = strain_matrix(latt_system, up)[0]
            stress_list1 = array(stress_set_dict[up][0])

            eplisons_now = array([[s1[0][0],  s1[1][1],  s1[2][2], 2*s1[0][1],  0.,   0.,   0.],
                                  [s1[1][1],  s1[0][0],  s1[2][2],-2*s1[0][1],  0.,   0.,   0.],
                                  [0.,    0., s1[0][0]+s1[1][1],      0.,   s1[2][2], 0.,   0.],
                                  [0.,    0.,     0.,     0.,     0.,     2*s1[1][2],       0.],
                                  [0.,    0.,     0.,     0.,     0.,     2*s1[0][2],       0.],
                                  [0.,    0.,     0., s1[0][0]-s1[1][1],   0.,   0.,2*s1[0][1]]])

            stresses_now = array([[stress_list1[0]],
                                  [stress_list1[1]],
                                  [stress_list1[2]],
                                  [stress_list1[3]],
                                  [stress_list1[4]],
                                  [stress_list1[5]]])
        else:
            s1 = strain_matrix(latt_system, up)[0]
            s2 = strain_matrix(latt_system, up)[1]
            stress_list1 = array(stress_set_dict[up][0])
            stress_list2 = array(stress_set_dict[up][1])

            eplisons_now = array([[s1[0][0],  s1[1][1],  s1[2][2], 2*s1[0][1],  0.,   0.,   0.],
                                  [s1[1][1],  s1[0][0],  s1[2][2],-2*s1[0][1],  0.,   0.,   0.],
                                  [0.,    0., s1[0][0]+s1[1][1],      0.,   s1[2][2], 0.,   0.],
                                  [0.,    0.,     0.,     0.,     0.,     2*s1[1][2],       0.],
                                  [0.,    0.,     0.,     0.,     0.,     2*s1[0][2],       0.],
                                  [0.,    0.,     0., s1[0][0]-s1[1][1],   0.,   0.,2*s1[0][1]],
                                  [s2[0][0],  s2[1][1],  s2[2][2], 2*s2[0][1],  0.,   0.,   0.],
                                  [s2[1][1],  s2[0][0],  s2[2][2],-2*s2[0][1],  0.,   0.,   0.],
                                  [0.,    0., s2[0][0]+s2[1][1],      0.,   s2[2][2], 0.,   0.],
                                  [0.,    0.,     0.,     0.,     0.,     2*s2[1][2],       0.],
                                  [0.,    0.,     0.,     0.,     0.,     2*s2[0][2],       0.],
                                  [0.,    0.,     0., s2[0][0]-s2[1][1],   0.,   0.,2*s2[0][1]]])

            stresses_now = array([[stress_list1[0]],
                                  [stress_list1[1]],
                                  [stress_list1[2]],
                                  [stress_list1[3]],
                                  [stress_list1[4]],
                                  [stress_list1[5]],
                                  [stress_list2[0]],
                                  [stress_list2[1]],
                                  [stress_list2[2]],
                                  [stress_list2[3]],
                                  [stress_list2[4]],
                                  [stress_list2[5]]])
        if n == 0:
            eplisons = deepcopy(eplisons_now)
            stresses = deepcopy(stresses_now)
        else:
            eplisons = vstack((eplisons, eplisons_now))
            stresses = vstack((stresses, stresses_now))

    cij = linalg.lstsq(eplisons, stresses, rcond=None)[0] * conversion_factor

    c11 = cij[0][0]
    c12 = cij[1][0]
    c13 = cij[2][0]
    c16 = cij[3][0]
    c33 = cij[4][0]
    c44 = cij[5][0]
    c66 = cij[6][0]

    M = c11+c12+2*c33-4*c13
    C2 = (c11+c12)*c33-2*c13*c13

    B_v = (2*(c11+c12)+c33+4*c13)/9.
    G_v = (M+3*c11-3*c12+12*c44+6*c66)/30.
    B_r = C2/M
    G_r = 15./(18*B_v/C2+6/(c11-c12)+6/c44+3/c66)

    B_vrh = (B_v+B_r)/2.
    G_vrh = (G_v+G_r)/2.
    E = 9*B_vrh*G_vrh/(3*B_vrh+G_vrh)
    v = (3*B_vrh-2*G_vrh)/(2*(3*B_vrh+G_vrh))

    elastic_constants_dict['c11'] = c11
    elastic_constants_dict['c12'] = c12
    elastic_constants_dict['c13'] = c13
    elastic_constants_dict['c16'] = c16
    elastic_constants_dict['c33'] = c33
    elastic_constants_dict['c44'] = c44
    elastic_constants_dict['c66'] = c66

    elastic_constants_dict['B_v'] = B_v
    elastic_constants_dict['B_r'] = B_r
    elastic_constants_dict['G_v'] = G_v
    elastic_constants_dict['G_r'] = G_r

    elastic_constants_dict['B_vrh'] = B_vrh
    elastic_constants_dict['G_vrh'] = G_vrh
    elastic_constants_dict['E'] = E
    elastic_constants_dict['v'] = v
    
    return elastic_constants_dict


def Orthorombic(latt_system, elastic_constants_dict, stress_set_dict, conversion_factor):
    for n, up in enumerate(stress_set_dict.keys()):
        if up == 0:
            s1 = strain_matrix(latt_system, up)[0]
            stress_list1 = array(stress_set_dict[up][0])

            eplisons_now = array([[s1[0][0], s1[1][1], s1[2][2], 0., 0., 0., 0., 0., 0.],
                                  [0., s1[0][0], 0., s1[1][1], s1[2][2], 0., 0., 0., 0.],
                                  [0., 0., s1[0][0], 0., s1[1][1], s1[2][2], 0., 0., 0.],
                                  [0.,   0.,   0.,   0.,   0.,  0., 2*s1[1][2],  0., 0.],
                                  [0.,   0.,   0.,   0.,   0.,  0.,  0.,2*s1[0][2],  0.],
                                  [0.,   0.,   0.,   0.,   0.,  0.,  0.,  0.,2*s1[0][1]]])

            stresses_now = array([[stress_list1[0]],
                                  [stress_list1[1]],
                                  [stress_list1[2]],
                                  [stress_list1[3]],
                                  [stress_list1[4]],
                                  [stress_list1[5]]])
        else:        
            s1 = strain_matrix(latt_system, up)[0]
            s2 = strain_matrix(latt_system, up)[1]
            s3 = strain_matrix(latt_system, up)[2]
            stress_list1 = array(stress_set_dict[up][0])
            stress_list2 = array(stress_set_dict[up][1])
            stress_list3 = array(stress_set_dict[up][2])

            eplisons_now = array([[s1[0][0], s1[1][1], s1[2][2], 0., 0., 0., 0., 0., 0.],
                                  [0., s1[0][0], 0., s1[1][1], s1[2][2], 0., 0., 0., 0.],
                                  [0., 0., s1[0][0], 0., s1[1][1], s1[2][2], 0., 0., 0.],
                                  [0.,   0.,   0.,   0.,   0.,  0., 2*s1[1][2],  0., 0.],
                                  [0.,   0.,   0.,   0.,   0.,  0.,  0.,2*s1[0][2],  0.],
                                  [0.,   0.,   0.,   0.,   0.,  0.,  0.,  0.,2*s1[0][1]],
                                  [s2[0][0], s2[1][1], s2[2][2], 0., 0., 0., 0., 0., 0.],
                                  [0., s2[0][0], 0., s2[1][1], s2[2][2], 0., 0., 0., 0.],
                                  [0., 0., s2[0][0], 0., s2[1][1], s2[2][2], 0., 0., 0.],
                                  [0.,   0.,   0.,   0.,   0.,   0.,2*s2[1][2],  0., 0.],
                                  [0.,   0.,   0.,   0.,   0.,   0.,  0.,2*s2[0][2], 0.],
                                  [0.,   0.,   0.,   0.,   0.,   0.,  0., 0.,2*s2[0][1]],
                                  [s3[0][0], s3[1][1], s3[2][2], 0., 0., 0., 0., 0., 0.],
                                  [0., s3[0][0], 0., s3[1][1], s3[2][2], 0., 0., 0., 0.],
                                  [0., 0., s3[0][0], 0., s3[1][1], s3[2][2], 0., 0., 0.],
                                  [0.,   0.,   0.,   0.,   0.,   0., 2*s3[1][2], 0., 0.],
                                  [0.,   0.,   0.,   0.,   0.,   0., 0., 2*s3[0][2], 0.],
                                  [0.,   0.,   0.,   0.,   0.,   0., 0., 0., 2*s3[0][1]]])

            stresses_now = array([[stress_list1[0]],
                                  [stress_list1[1]],
                                  [stress_list1[2]],
                                  [stress_list1[3]],
                                  [stress_list1[4]],
                                  [stress_list1[5]],
                                  [stress_list2[0]],
                                  [stress_list2[1]],
                                  [stress_list2[2]],
                                  [stress_list2[3]],
                                  [stress_list2[4]],
                                  [stress_list2[5]],
                                  [stress_list3[0]],
                                  [stress_list3[1]],
                                  [stress_list3[2]],
                                  [stress_list3[3]],
                                  [stress_list3[4]],
                                  [stress_list3[5]]])
        if n == 0:
            eplisons = deepcopy(eplisons_now)
            stresses = deepcopy(stresses_now)
        else:
            eplisons = vstack((eplisons, eplisons_now))
            stresses = vstack((stresses, stresses_now))

    cij = linalg.lstsq(eplisons, stresses, rcond=None)[0] * conversion_factor
    c11 = cij[0][0]
    c12 = cij[1][0]
    c13 = cij[2][0]
    c22 = cij[3][0]
    c23 = cij[4][0]
    c33 = cij[5][0]
    c44 = cij[6][0]
    c55 = cij[7][0]
    c66 = cij[8][0]

    D = c13*(c12*c23-c13*c22)+c23*(c12*c13-c23*c11)+c33*(c11*c22-c12*c12)
    B_v = (c11+c22+c33+2*(c12+c13+c23))/9.
    G_v = (c11+c22+c33+3*(c44+c55+c66)-(c12+c13+c23))/15.
    B_r = D/(c11*(c22+c33-2*c23)+c22*(c33-2*c13)-2*c33*c12+c12*(2*c23-c12)+c13*(2*c12-c13)+c23*(2*c13-c23))
    G_r = 15/(4*(c11*(c22+c33+c23)+c22*(c33+c13)+c33*c12-c12*(c23+c12)-c13*(c12+c13)-c23*(c13+c23))/D+3*(1/c44+1/c55+1/c66))

    B_vrh = (B_v+B_r)/2.
    G_vrh = (G_v+G_r)/2.
    E = 9*B_vrh*G_vrh/(3*B_vrh+G_vrh)
    v = (3*B_vrh-2*G_vrh)/(2*(3*B_vrh+G_vrh))

    elastic_constants_dict['c11'] = c11
    elastic_constants_dict['c12'] = c12
    elastic_constants_dict['c13'] = c13
    elastic_constants_dict['c22'] = c22
    elastic_constants_dict['c23'] = c23
    elastic_constants_dict['c33'] = c33
    elastic_constants_dict['c44'] = c44
    elastic_constants_dict['c55'] = c55
    elastic_constants_dict['c66'] = c66

    elastic_constants_dict['B_v'] = B_v
    elastic_constants_dict['B_r'] = B_r
    elastic_constants_dict['G_v'] = G_v
    elastic_constants_dict['G_r'] = G_r
    elastic_constants_dict['B_vrh'] = B_vrh
    elastic_constants_dict['G_vrh'] = G_vrh
    elastic_constants_dict['E'] = E
    elastic_constants_dict['v'] = v

    return elastic_constants_dict


def Monoclinic(latt_system, elastic_constants_dict, stress_set_dict, conversion_factor):
    for n, up in enumerate(stress_set_dict.keys()):
        if up == 0:
            s1 = strain_matrix(latt_system, up)[0]
            stress_list1 = array(stress_set_dict[up][0])

            eplisons_now = array([[s1[0][0], s1[1][1], s1[2][2], 2*s1[0][2], 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                                  [0., s1[0][0], 0., 0., s1[1][1], s1[2][2], 2*s1[0][2], 0., 0., 0., 0., 0., 0.],
                                  [0., 0., s1[0][0], 0.,  0., s1[1][1], 0.,s1[2][2], 2*s1[0][2], 0., 0., 0., 0.],
                                  [0., 0.,  0.,  0., 0.,  0.,  0.,  0.,  0.,  2*s1[1][2],   2*s1[0][1],  0., 0.],
                                  [0., 0.,  0.,s1[0][0], 0.,  0.,s1[1][1],  0., s1[2][2],  0., 0.,2*s1[0][2],0.],
                                  [0., 0.,  0., 0.,  0., 0.,  0.,   0.,  0.,  0.,   2*s1[1][2], 0.,  2*s1[0][1]]])

            stresses_now = array([[stress_list1[0]],
                                  [stress_list1[1]],
                                  [stress_list1[2]],
                                  [stress_list1[3]],
                                  [stress_list1[4]],
                                  [stress_list1[5]]])
        else:        
            s1 = strain_matrix(latt_system, up)[0]
            s2 = strain_matrix(latt_system, up)[1]
            s3 = strain_matrix(latt_system, up)[2]
            s4 = strain_matrix(latt_system, up)[3]
            #s5 = strain_matrix(latt_system, up)[4]
            stress_list1 = array(stress_set_dict[up][0])
            stress_list2 = array(stress_set_dict[up][1])
            stress_list3 = array(stress_set_dict[up][2])
            stress_list4 = array(stress_set_dict[up][3])
            #stress_list5 = array(stress_set_dict[up][4])

            eplisons_now = array([[s1[0][0], s1[1][1], s1[2][2], 2*s1[0][2], 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                                  [0., s1[0][0], 0., 0., s1[1][1], s1[2][2], 2*s1[0][2], 0., 0., 0., 0., 0., 0.],
                                  [0., 0., s1[0][0], 0.,  0., s1[1][1], 0.,s1[2][2], 2*s1[0][2], 0., 0., 0., 0.],
                                  [0., 0.,  0.,  0., 0.,  0.,  0.,  0.,  0.,  2*s1[1][2],   2*s1[0][1],  0., 0.],
                                  [0., 0.,  0.,s1[0][0], 0.,  0.,s1[1][1],  0., s1[2][2],  0., 0.,2*s1[0][2],0.],
                                  [0., 0.,  0., 0.,  0., 0.,  0.,   0.,  0.,  0.,   2*s1[1][2], 0.,  2*s1[0][1]],
                                  [s2[0][0], s2[1][1], s2[2][2], 2*s2[0][2], 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                                  [0., s2[0][0], 0., 0., s2[1][1], s2[2][2], 2*s2[0][2], 0., 0., 0., 0., 0., 0.],
                                  [0., 0., s2[0][0], 0.,  0., s2[1][1], 0.,s2[2][2], 2*s2[0][2], 0., 0., 0., 0.],
                                  [0., 0.,  0.,  0., 0.,  0.,  0.,  0.,  0.,  2*s2[1][2],   2*s2[0][1],  0., 0.],
                                  [0., 0.,  0.,s2[0][0], 0.,  0.,s2[1][1],  0., s2[2][2],  0., 0.,2*s2[0][2],0.],
                                  [0., 0.,  0., 0.,  0., 0.,  0.,   0.,  0.,  0.,   2*s2[1][2], 0.,  2*s2[0][1]],
                                  [s3[0][0], s3[1][1], s3[2][2], 2*s3[0][2], 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                                  [0., s3[0][0], 0., 0., s3[1][1], s3[2][2], 2*s3[0][2], 0., 0., 0., 0., 0., 0.],
                                  [0., 0., s3[0][0], 0.,  0., s3[1][1], 0.,s3[2][2], 2*s3[0][2], 0., 0., 0., 0.],
                                  [0., 0.,  0.,  0., 0.,  0.,  0.,  0.,  0.,  2*s3[1][2],   2*s3[0][1],  0., 0.],
                                  [0., 0.,  0.,s3[0][0], 0.,  0.,s3[1][1],  0., s3[2][2],  0., 0.,2*s3[0][2],0.],
                                  [0., 0.,  0., 0.,  0., 0.,  0.,   0.,  0.,  0.,   2*s3[1][2], 0.,  2*s3[0][1]],
                                  [s4[0][0], s4[1][1], s4[2][2], 2*s4[0][2], 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                                  [0., s4[0][0], 0., 0., s4[1][1], s4[2][2], 2*s4[0][2], 0., 0., 0., 0., 0., 0.],
                                  [0., 0., s4[0][0], 0.,  0., s4[1][1], 0.,s4[2][2], 2*s4[0][2], 0., 0., 0., 0.],
                                  [0., 0.,  0.,  0., 0.,  0.,  0.,  0.,  0.,  2*s4[1][2],   2*s4[0][1],  0., 0.],
                                  [0., 0.,  0.,s4[0][0], 0.,  0.,s4[1][1],  0., s4[2][2],  0., 0.,2*s4[0][2],0.],
                                  [0., 0.,  0., 0.,  0., 0.,  0.,   0.,  0.,  0.,   2*s4[1][2], 0.,  2*s4[0][1]]])

            stresses_now = array([[stress_list1[0]],
                                  [stress_list1[1]],
                                  [stress_list1[2]],
                                  [stress_list1[3]],
                                  [stress_list1[4]],
                                  [stress_list1[5]],
                                  [stress_list2[0]],
                                  [stress_list2[1]],
                                  [stress_list2[2]],
                                  [stress_list2[3]],
                                  [stress_list2[4]],
                                  [stress_list2[5]],
                                  [stress_list3[0]],
                                  [stress_list3[1]],
                                  [stress_list3[2]],
                                  [stress_list3[3]],
                                  [stress_list3[4]],
                                  [stress_list3[5]],
                                  [stress_list4[0]],
                                  [stress_list4[1]],
                                  [stress_list4[2]],
                                  [stress_list4[3]],
                                  [stress_list4[4]],
                                  [stress_list4[5]]])
        if n == 0:
            eplisons = deepcopy(eplisons_now)
            stresses = deepcopy(stresses_now)
        else:
            eplisons = vstack((eplisons, eplisons_now))
            stresses = vstack((stresses, stresses_now))
            
    cij = linalg.lstsq(eplisons, stresses, rcond=None)[0] * conversion_factor
    c11 = cij[0][0]
    c12 = cij[1][0]
    c13 = cij[2][0]
    c15 = cij[3][0]
    c22 = cij[4][0]
    c23 = cij[5][0]
    c25 = cij[6][0]
    c33 = cij[7][0]
    c35 = cij[8][0]
    c44 = cij[9][0]
    c46 = cij[10][0]
    c55 = cij[11][0]
    c66 = cij[12][0]

    a = c33*c55-c35*c35
    b = c23*c55-c25*c35
    c = c13*c35-c15*c33
    d = c13*c55-c15*c35
    e = c13*c25-c15*c23
    f = c11*(c22*c55-c25*c25)-c12*(c12*c55-c15*c25)+c15*(c12*c25-c15*c22)+c25*(c23*c35-c25*c33)
    g = c11*c22*c33-c11*c23*c23-c22*c13*c13-c33*c12*c12+2*c12*c13*c23
    O = 2*(c15*c25*(c33*c12-c13*c23)+c15*c35*(c22*c13-c12*c23)+c25*c35*(c11*c23-c12*c13))-(c15*c15*(c22*c33-c23*c23)+c25*c25*(c11*c33-c13*c13)+c35*c35*(c11*c22-c12*c12))+g*c55

    B_v = (c11+c22+c33+2*(c12+c13+c23))/9.
    G_v = (c11+c22+c33+3*(c44+c55+c66)-(c12+c13+c23))/15.
    B_r = O/(a*(c11+c22-2*c12)+b*(2*c12-2*c11-c23)+c*(c15-2*c25)+d*(2*c12+2*c23-c13-2*c22)+2*e*(c25-c15)+f)
    G_r = 15/(4*(a*(c11+c22+c12)+b*(c11-c12-c23)+c*(c15+c25)+d*(c22-c12-c23-c13)+e*(c15-c25)+f)/O+3*(g/O+(c44+c66)/(c44*c66-c46*c46)))

    B_vrh = (B_v+B_r)/2.
    G_vrh = (G_v+G_r)/2.
    E = 9*B_vrh*G_vrh/(3*B_vrh+G_vrh)
    v = (3*B_vrh-2*G_vrh)/(2*(3*B_vrh+G_vrh))

    elastic_constants_dict['c11'] = c11
    elastic_constants_dict['c12'] = c12
    elastic_constants_dict['c13'] = c13
    elastic_constants_dict['c15'] = c15
    elastic_constants_dict['c22'] = c22
    elastic_constants_dict['c23'] = c23
    elastic_constants_dict['c25'] = c25
    elastic_constants_dict['c33'] = c33
    elastic_constants_dict['c35'] = c35
    elastic_constants_dict['c44'] = c44
    elastic_constants_dict['c46'] = c46
    elastic_constants_dict['c55'] = c55
    elastic_constants_dict['c66'] = c66

    elastic_constants_dict['B_v'] = B_v
    elastic_constants_dict['B_r'] = B_r
    elastic_constants_dict['G_v'] = G_v
    elastic_constants_dict['G_r'] = G_r
    elastic_constants_dict['B_vrh'] = B_vrh
    elastic_constants_dict['G_vrh'] = G_vrh
    elastic_constants_dict['E'] = E
    elastic_constants_dict['v'] = v


    return elastic_constants_dict


def Triclinic(latt_system, elastic_constants_dict, stress_set_dict, conversion_factor):
    for n, up in enumerate(stress_set_dict.keys()):
        if up == 0:
            s1 = strain_matrix(latt_system, up)[0]
            stress_list1 = array(stress_set_dict[up][0])

            eplisons_now = array([[s1[0][0], s1[1][1], s1[2][2],2*s1[1][2],2*s1[0][2],2*s1[0][1], 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                                  [0., s1[0][0], 0., 0., 0., 0., s1[1][1], s1[2][2],2*s1[1][2],2*s1[0][2],2*s1[0][1], 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                                  [0., 0., s1[0][0], 0., 0., 0., 0., s1[1][1], 0., 0., 0., s1[2][2],2*s1[1][2],2*s1[0][2],2*s1[0][1], 0., 0., 0., 0., 0., 0.],
                                  [0., 0., 0., s1[0][0], 0., 0., 0., 0., s1[1][1], 0., 0., 0., s1[2][2], 0., 0.,2*s1[1][2],2*s1[0][2],2*s1[0][1], 0., 0., 0.],
                                  [0., 0., 0., 0., s1[0][0], 0., 0., 0., 0., s1[1][1], 0., 0., 0., s1[2][2], 0., 0.,2*s1[1][2], 0.,2*s1[0][2],2*s1[0][1], 0.],
                                  [0., 0., 0., 0., 0., s1[0][0], 0., 0., 0., 0., s1[1][1], 0., 0., 0., s1[2][2], 0.,0.,2*s1[1][2], 0.,2*s1[0][2], 2*s1[0][1]]])

            stresses_now = array([[stress_list1[0]],
                                     [stress_list1[1]],
                                     [stress_list1[2]],
                                     [stress_list1[3]],
                                     [stress_list1[4]],
                                     [stress_list1[5]]])
        else:        
            s1 = strain_matrix(latt_system, up)[0]
            s2 = strain_matrix(latt_system, up)[1]
            s3 = strain_matrix(latt_system, up)[2]
            s4 = strain_matrix(latt_system, up)[3]
            s5 = strain_matrix(latt_system, up)[4]
            s6 = strain_matrix(latt_system, up)[5]
            stress_list1 = array(stress_set_dict[up][0])
            stress_list2 = array(stress_set_dict[up][1])
            stress_list3 = array(stress_set_dict[up][2])
            stress_list4 = array(stress_set_dict[up][3])
            stress_list5 = array(stress_set_dict[up][4])
            stress_list6 = array(stress_set_dict[up][5])

            eplisons_now = array([[s1[0][0], s1[1][1], s1[2][2],2*s1[1][2],2*s1[0][2],2*s1[0][1], 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                                 [0., s1[0][0], 0., 0., 0., 0., s1[1][1], s1[2][2],2*s1[1][2],2*s1[0][2],2*s1[0][1], 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                                 [0., 0., s1[0][0], 0., 0., 0., 0., s1[1][1], 0., 0., 0., s1[2][2],2*s1[1][2],2*s1[0][2],2*s1[0][1], 0., 0., 0., 0., 0., 0.],
                                 [0., 0., 0., s1[0][0], 0., 0., 0., 0., s1[1][1], 0., 0., 0., s1[2][2], 0., 0.,2*s1[1][2],2*s1[0][2],2*s1[0][1], 0., 0., 0.],
                                 [0., 0., 0., 0., s1[0][0], 0., 0., 0., 0., s1[1][1], 0., 0., 0., s1[2][2], 0., 0.,2*s1[1][2], 0.,2*s1[0][2],2*s1[0][1], 0.],
                                 [0., 0., 0., 0., 0., s1[0][0], 0., 0., 0., 0., s1[1][1], 0., 0., 0., s1[2][2], 0.,0.,2*s1[1][2], 0.,2*s1[0][2], 2*s1[0][1]],
                                 [s2[0][0], s2[1][1], s2[2][2],2*s2[1][2],2*s2[0][2],2*s2[0][1], 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                                 [0., s2[0][0], 0., 0., 0., 0., s2[1][1], s2[2][2],2*s2[1][2],2*s2[0][2],2*s2[0][1], 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                                 [0., 0., s2[0][0], 0., 0., 0., 0., s2[1][1], 0., 0., 0., s2[2][2],2*s2[1][2],2*s2[0][2],2*s2[0][1], 0., 0., 0., 0., 0., 0.],
                                 [0., 0., 0., s2[0][0], 0., 0., 0., 0., s2[1][1], 0., 0., 0., s2[2][2], 0., 0.,2*s2[1][2],2*s2[0][2],2*s2[0][1], 0., 0., 0.],
                                 [0., 0., 0., 0., s2[0][0], 0., 0., 0., 0., s2[1][1], 0., 0., 0., s2[2][2], 0., 0.,2*s2[1][2], 0.,2*s2[0][2],2*s2[0][1], 0.],
                                 [0., 0., 0., 0., 0., s2[0][0], 0., 0., 0., 0., s2[1][1], 0., 0., 0., s2[2][2], 0.,0.,2*s2[1][2], 0.,2*s2[0][2], 2*s2[0][1]],
                                 [s3[0][0], s3[1][1], s3[2][2],2*s3[1][2],2*s3[0][2],2*s3[0][1], 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                                 [0., s3[0][0], 0., 0., 0., 0., s3[1][1], s3[2][2],2*s3[1][2],2*s3[0][2],2*s3[0][1], 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                                 [0., 0., s3[0][0], 0., 0., 0., 0., s3[1][1], 0., 0., 0., s3[2][2],2*s3[1][2],2*s3[0][2],2*s3[0][1], 0., 0., 0., 0., 0., 0.],
                                 [0., 0., 0., s3[0][0], 0., 0., 0., 0., s3[1][1], 0., 0., 0., s3[2][2], 0., 0.,2*s3[1][2],2*s3[0][2],2*s3[0][1], 0., 0., 0.],
                                 [0., 0., 0., 0., s3[0][0], 0., 0., 0., 0., s3[1][1], 0., 0., 0., s3[2][2], 0., 0.,2*s3[1][2], 0.,2*s3[0][2],2*s3[0][1], 0.],
                                 [0., 0., 0., 0., 0., s3[0][0], 0., 0., 0., 0., s3[1][1], 0., 0., 0., s3[2][2], 0.,0.,2*s3[1][2], 0.,2*s3[0][2], 2*s3[0][1]],
                                 [s4[0][0], s4[1][1], s4[2][2],2*s4[1][2],2*s4[0][2],2*s4[0][1], 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                                 [0., s4[0][0], 0., 0., 0., 0., s4[1][1], s4[2][2],2*s4[1][2],2*s4[0][2],2*s4[0][1], 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                                 [0., 0., s4[0][0], 0., 0., 0., 0., s4[1][1], 0., 0., 0., s4[2][2],2*s4[1][2],2*s4[0][2],2*s4[0][1], 0., 0., 0., 0., 0., 0.],
                                 [0., 0., 0., s4[0][0], 0., 0., 0., 0., s4[1][1], 0., 0., 0., s4[2][2], 0., 0.,2*s4[1][2],2*s4[0][2],2*s4[0][1], 0., 0., 0.],
                                 [0., 0., 0., 0., s4[0][0], 0., 0., 0., 0., s4[1][1], 0., 0., 0., s4[2][2], 0., 0.,2*s4[1][2], 0.,2*s4[0][2],2*s4[0][1], 0.],
                                 [0., 0., 0., 0., 0., s4[0][0], 0., 0., 0., 0., s4[1][1], 0., 0., 0., s4[2][2], 0.,0.,2*s4[1][2], 0.,2*s4[0][2], 2*s4[0][1]],
                                 [s5[0][0], s5[1][1], s5[2][2],2*s5[1][2],2*s5[0][2],2*s5[0][1], 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                                 [0., s5[0][0], 0., 0., 0., 0., s5[1][1], s5[2][2],2*s5[1][2],2*s5[0][2],2*s5[0][1], 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                                 [0., 0., s5[0][0], 0., 0., 0., 0., s5[1][1], 0., 0., 0., s5[2][2],2*s5[1][2],2*s5[0][2],2*s5[0][1], 0., 0., 0., 0., 0., 0.],
                                 [0., 0., 0., s5[0][0], 0., 0., 0., 0., s5[1][1], 0., 0., 0., s5[2][2], 0., 0.,2*s5[1][2],2*s5[0][2],2*s5[0][1], 0., 0., 0.],
                                 [0., 0., 0., 0., s5[0][0], 0., 0., 0., 0., s5[1][1], 0., 0., 0., s5[2][2], 0., 0.,2*s5[1][2], 0.,2*s5[0][2],2*s5[0][1], 0.],
                                 [0., 0., 0., 0., 0., s5[0][0], 0., 0., 0., 0., s5[1][1], 0., 0., 0., s5[2][2], 0.,0.,2*s5[1][2], 0.,2*s5[0][2], 2*s5[0][1]],
                                 [s6[0][0], s6[1][1], s6[2][2],2*s6[1][2],2*s6[0][2],2*s6[0][1], 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                                 [0., s6[0][0], 0., 0., 0., 0., s6[1][1], s6[2][2],2*s6[1][2],2*s6[0][2],2*s6[0][1], 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                                 [0., 0., s6[0][0], 0., 0., 0., 0., s6[1][1], 0., 0., 0., s6[2][2],2*s6[1][2],2*s6[0][2],2*s6[0][1], 0., 0., 0., 0., 0., 0.],
                                 [0., 0., 0., s6[0][0], 0., 0., 0., 0., s6[1][1], 0., 0., 0., s6[2][2], 0., 0.,2*s6[1][2],2*s6[0][2],2*s6[0][1], 0., 0., 0.],
                                 [0., 0., 0., 0., s6[0][0], 0., 0., 0., 0., s6[1][1], 0., 0., 0., s6[2][2], 0., 0.,2*s6[1][2], 0.,2*s6[0][2],2*s6[0][1], 0.],
                                 [0., 0., 0., 0., 0., s6[0][0], 0., 0., 0., 0., s6[1][1], 0., 0., 0., s6[2][2], 0.,0.,2*s6[1][2], 0.,2*s6[0][2], 2*s6[0][1]]])

            stresses_now = array([[stress_list1[0]],
                                  [stress_list1[1]],
                                  [stress_list1[2]],
                                  [stress_list1[3]],
                                  [stress_list1[4]],
                                  [stress_list1[5]],
                                  [stress_list2[0]],
                                  [stress_list2[1]],
                                  [stress_list2[2]],
                                  [stress_list2[3]],
                                  [stress_list2[4]],
                                  [stress_list2[5]],
                                  [stress_list3[0]],
                                  [stress_list3[1]],
                                  [stress_list3[2]],
                                  [stress_list3[3]],
                                  [stress_list3[4]],
                                  [stress_list3[5]],
                                  [stress_list4[0]],
                                  [stress_list4[1]],
                                  [stress_list4[2]],
                                  [stress_list4[3]],
                                  [stress_list4[4]],
                                  [stress_list4[5]],
                                  [stress_list5[0]],
                                  [stress_list5[1]],
                                  [stress_list5[2]],
                                  [stress_list5[3]],
                                  [stress_list5[4]],
                                  [stress_list5[5]],
                                  [stress_list6[0]],
                                  [stress_list6[1]],
                                  [stress_list6[2]],
                                  [stress_list6[3]],
                                  [stress_list6[4]],
                                  [stress_list6[5]]])
        if n == 0:
            eplisons = deepcopy(eplisons_now)
            stresses = deepcopy(stresses_now)
        else:
            eplisons = vstack((eplisons, eplisons_now))
            stresses = vstack((stresses, stresses_now))
            
    cij = linalg.lstsq(eplisons, stresses, rcond=None)[0] * conversion_factor

    c11 = cij[0][0]
    c12 = cij[1][0]
    c13 = cij[2][0]
    c14 = cij[3][0]
    c15 = cij[4][0]
    c16 = cij[5][0]
    c22 = cij[6][0]
    c23 = cij[7][0]
    c24 = cij[8][0]
    c25 = cij[9][0]
    c26 = cij[10][0]
    c33 = cij[11][0]
    c34 = cij[12][0]
    c35 = cij[13][0]
    c36 = cij[14][0]
    c44 = cij[15][0]
    c45 = cij[16][0]
    c46 = cij[17][0]
    c55 = cij[18][0]
    c56 = cij[19][0]
    c66 = cij[20][0]

    elastic_constants_dict['c11'] = c11
    elastic_constants_dict['c12'] = c12
    elastic_constants_dict['c13'] = c13
    elastic_constants_dict['c14'] = c14
    elastic_constants_dict['c15'] = c15
    elastic_constants_dict['c16'] = c16
    elastic_constants_dict['c22'] = c22
    elastic_constants_dict['c23'] = c23
    elastic_constants_dict['c24'] = c24
    elastic_constants_dict['c25'] = c25
    elastic_constants_dict['c26'] = c26
    elastic_constants_dict['c33'] = c33
    elastic_constants_dict['c34'] = c34
    elastic_constants_dict['c35'] = c35
    elastic_constants_dict['c36'] = c36
    elastic_constants_dict['c44'] = c44
    elastic_constants_dict['c45'] = c45
    elastic_constants_dict['c46'] = c46
    elastic_constants_dict['c55'] = c55
    elastic_constants_dict['c56'] = c56
    elastic_constants_dict['c66'] = c66

    return elastic_constants_dict


def isotropy(pos_optimized, latt_system, elastic_constants_dict, stress_set_dict, conversion_factor):
    for n, up in enumerate(stress_set_dict.keys()):
        s = strain_matrix(latt_system, up)[0]
        #print(s)

        stress_list = array(stress_set_dict[up][0])

        eplisons_now = array([[s[0][0],   s[1][1]],
                              [s[1][1],   s[0][0]],
                              [s[0][1],  -s[0][1]]])

        stresses_now = array([[stress_list[0]],
                              [stress_list[1]],
                              [stress_list[5]]])
        if n == 0:
            eplisons = deepcopy(eplisons_now)
            stresses = deepcopy(stresses_now)
        else:
            eplisons = vstack((eplisons, eplisons_now))
            stresses = vstack((stresses, stresses_now))

    cij = linalg.lstsq(eplisons, stresses, rcond=None)[0] * conversion_factor
    lenth_angl = pos_optimized.get_cell_lengths_and_angles()

    elastic_constants_dict['c11'] = cij[0][0] * lenth_angl[2] / 10.
    elastic_constants_dict['c12'] = cij[1][0] * lenth_angl[2] / 10.

    return elastic_constants_dict


def tetragonal(pos_optimized, latt_system, elastic_constants_dict, stress_set_dict, conversion_factor):
    for n, up in enumerate(stress_set_dict.keys()):
        s = strain_matrix(latt_system, up)[0]

        stress_list = array(stress_set_dict[up][0])

        eplisons_now = array([[s[0][0], s[1][1], 0.],
                              [s[1][1], s[0][0], 0.],
                              [0.,   0.,  2*s[0][1]]])

        stresses_now = array([[stress_list[0]],
                              [stress_list[1]],
                              [stress_list[5]]])
        if n == 0:
            eplisons = deepcopy(eplisons_now)
            stresses = deepcopy(stresses_now)
        else:
            eplisons = vstack((eplisons, eplisons_now))
            stresses = vstack((stresses, stresses_now))

    cij = linalg.lstsq(eplisons, stresses, rcond=None)[0] * conversion_factor
    lenth_angl = pos_optimized.get_cell_lengths_and_angles()

    elastic_constants_dict['c11'] = cij[0][0] * lenth_angl[2] / 10.
    elastic_constants_dict['c12'] = cij[1][0] * lenth_angl[2] / 10.
    elastic_constants_dict['c66'] = cij[2][0] * lenth_angl[2] / 10.

    return elastic_constants_dict


def orthotropy(pos_optimized, latt_system, elastic_constants_dict, stress_set_dict, conversion_factor):
    for n, up in enumerate(stress_set_dict.keys()):
        if up == 0:
            s1 = strain_matrix(latt_system, up)[0]
            stress_list1 = array(stress_set_dict[up][0])
            eplisons_now = array([[s1[0][0], s1[1][1], 0., 0.],
                                  [0., s1[0][0], s1[1][1], 0.],
                                  [0.,   0.,  0.,  2*s1[0][1]]])
            
            stresses_now = array([[stress_list1[0]],
                                  [stress_list1[1]],
                                  [stress_list1[5]]])
        else:
            s1 = strain_matrix(latt_system, up)[0]
            s2 = strain_matrix(latt_system, up)[1]

            stress_list1 = array(stress_set_dict[up][0])
            stress_list2 = array(stress_set_dict[up][1])

            eplisons_now = array([[s1[0][0], s1[1][1], 0., 0.],
                                  [0., s1[0][0], s1[1][1], 0.],
                                  [0.,   0.,  0.,  2*s1[0][1]],
                                  [s2[0][0], s2[1][1], 0., 0.],
                                  [0., s2[0][0], s2[1][1], 0.],
                                  [0.,   0.,  0.,  2*s2[0][1]]])

            stresses_now = array([[stress_list1[0]],
                                  [stress_list1[1]],
                                  [stress_list1[5]],
                                  [stress_list2[0]],
                                  [stress_list2[1]],
                                  [stress_list2[5]]])
        if n == 0:
            eplisons = deepcopy(eplisons_now)
            stresses = deepcopy(stresses_now)
        else:
            eplisons = vstack((eplisons, eplisons_now))
            stresses = vstack((stresses, stresses_now))

    cij = linalg.lstsq(eplisons, stresses, rcond=None)[0] * conversion_factor
    lenth_angl = pos_optimized.get_cell_lengths_and_angles()

    elastic_constants_dict['c11'] = cij[0][0] * lenth_angl[2] / 10.
    elastic_constants_dict['c12'] = cij[1][0] * lenth_angl[2] / 10.
    elastic_constants_dict['c22'] = cij[2][0] * lenth_angl[2] / 10.
    elastic_constants_dict['c66'] = cij[3][0] * lenth_angl[2] / 10.

    return elastic_constants_dict


def anisotropy(pos_optimized, latt_system, elastic_constants_dict, stress_set_dict, conversion_factor):
    for n, up in enumerate(stress_set_dict.keys()):
        if up == 0:
            s1 = strain_matrix(latt_system, up)[0]
            stress_list1 = array(stress_set_dict[up][0])

            eplisons_now = array([[s1[0][0], s1[1][1], 2*s1[0][1], 0., 0., 0.],
                                  [0., s1[0][0], 0., s1[1][1], 2*s1[0][1], 0.],
                                  [0., 0., s1[0][0], 0., s1[1][1], 2*s1[0][1]]])

            stresses_now = array([[stress_list1[0]],
                                  [stress_list1[1]],
                                  [stress_list1[5]]])
        else:
            s1 = strain_matrix(latt_system, up)[0]
            s2 = strain_matrix(latt_system, up)[1]
            s3 = strain_matrix(latt_system, up)[2]

            stress_list1 = array(stress_set_dict[up][0])
            stress_list2 = array(stress_set_dict[up][1])
            stress_list3 = array(stress_set_dict[up][2])

            eplisons_now = array([[s1[0][0], s1[1][1], 2*s1[0][1], 0., 0., 0.],
                                  [0., s1[0][0], 0., s1[1][1], 2*s1[0][1], 0.],
                                  [0., 0., s1[0][0], 0., s1[1][1], 2*s1[0][1]],
                                  [s2[0][0], s2[1][1], 2*s2[0][1], 0., 0., 0.],
                                  [0., s2[0][0], 0., s2[1][1], 2*s2[0][1], 0.],
                                  [0., 0., s2[0][0], 0., s2[1][1], 2*s2[0][1]],
                                  [s3[0][0], s3[1][1], 2*s3[0][1], 0., 0., 0.],
                                  [0., s3[0][0], 0., s3[1][1], 2*s3[0][1], 0.],
                                  [0., 0., s3[0][0], 0., s3[1][1], 2*s3[0][1]]])

            stresses_now = array([[stress_list1[0]],
                                  [stress_list1[1]],
                                  [stress_list1[5]],
                                  [stress_list2[0]],
                                  [stress_list2[1]],
                                  [stress_list2[5]],
                                  [stress_list3[0]],
                                  [stress_list3[1]],
                                  [stress_list3[5]]])
        if n == 0:
            eplisons = deepcopy(eplisons_now)
            stresses = deepcopy(stresses_now)
        else:
            eplisons = vstack((eplisons, eplisons_now))
            stresses = vstack((stresses, stresses_now))

    cij = linalg.lstsq(eplisons, stresses, rcond=None)[0] * conversion_factor
    lenth_angl = pos_optimized.get_cell_lengths_and_angles()

    elastic_constants_dict['c11'] = cij[0][0] * lenth_angl[2] / 10.
    elastic_constants_dict['c12'] = cij[1][0] * lenth_angl[2] / 10.
    elastic_constants_dict['c16'] = cij[2][0] * lenth_angl[2] / 10.
    elastic_constants_dict['c22'] = cij[3][0] * lenth_angl[2] / 10.
    elastic_constants_dict['c26'] = cij[4][0] * lenth_angl[2] / 10.
    elastic_constants_dict['c66'] = cij[5][0] * lenth_angl[2] / 10.

    return elastic_constants_dict


def any1D(pos_optimized, latt_system, elastic_constants_dict, stress_set_dict, convertion_factor):
    stresses = []
    strains = []
        

    for key in stress_set_dict.keys():
        s = strain_matrix(latt_system, key)[0]
        current_stress_list = array(stress_set_dict[key][0])
                              
        current_strains = array([[s[0][0], s[0][1], s[0][2]],  # longitudinal
                                     [s[1][0], s[1][1], s[1][2]],  # circumferential
                                     [s[2][0], s[2][1], s[2][2]]]) # torsional

        current_stresses = np.array([[current_stress_list[0]],  # Longitudinal stress (σzz)
                                     [current_stress_list[1]],  # Circumferential stress (σrr)
                                     [current_stress_list[2]]])  # Shear stress (τrz)
                              
                              
                              
        strains.append(current_strains)
        stresses.append(current_stresses)
        
            
    strains = vstack(strains)
    stresses = vstack(stresses)
            
            
    cij = linalg.lstsq(strains, stresses, rcond=None)[0] * convertion_factor
    c33 = cij[2][0]
    c23 = cij[1][0]


    lenth_angl = pos_optimized.get_cell_lengths_and_angles()

    elastic_constants_dict['c33'] = c33 #* lenth_angl[0] * lenth_angl[1] / 100.
    elastic_constants_dict['c23'] = c23 #* lenth_angl[0] * lenth_angl[1] / 100.


#    elastic_constants_dict['E'] = E
#    elastic_constants_dict['v'] = v

    return elastic_constants_dict


def any1DModify(pos_optimized, latt_system, elastic_constants_dict, stress_set_dict, conversion_factor):
    stresses = []
    strains = []

    
    for key in stress_set_dict.keys():
        s = strain_matrix(latt_system, key)[0]
        current_stress_list = array(stress_set_dict[key][0])

#        current_strains = array([[s[1][1], s[2][2]],
#                                 [s[2][2], s[1][1]]])
#        current_stresses = array([[current_stress_list[1]],
#                                  [current_stress_list[2]]])


        current_strains = array([[s[0][0], s[0][1], s[0][2]],  # longitudinal
                                     [s[1][0], s[1][1], s[1][2]],  # circumferential
                                     [s[2][0], s[2][1], s[2][2]]]) # torsional

        current_stresses = np.array([[current_stress_list[0]],  # Longitudinal stress (σzz)
                                     [current_stress_list[1]],  # Circumferential stress (σrr)
                                     [current_stress_list[2]]])  # Shear stress (τrz)

        strains.append(current_strains)
        stresses.append(current_stresses)



    strains = vstack(strains)
    stresses = vstack(stresses)


    cij = linalg.lstsq(strains, stresses, rcond=None)[0] * conversion_factor
    print(cij)
    #exit(1)
    c22 = cij[2][0]
    c23 = cij[1][0]
    # E = c22 - c23*c23/c22
    # v = c23/c22

    lenth_angl = pos_optimized.get_cell_lengths_and_angles()

    elastic_constants_dict['c22'] = c22 #* lenth_angl[0] * lenth_angl[1] / 100.
    elastic_constants_dict['c23'] = c23 #* lenth_angl[0] * lenth_angl[1] / 100.


#    elastic_constants_dict['E'] = E
#    elastic_constants_dict['v'] = v

    return elastic_constants_dict



def true1D(pos_optimized, latt_system, elastic_constants_dict, stress_set_dict, conversion_factor):
    stresses = []
    strains = []

    for key in stress_set_dict.keys():
        s = strain_matrix(latt_system, key)[0]  # Getting the strain matrix
        current_stress_list = array(stress_set_dict[key][0])

        # Considering only the deformation in the 1D system
        current_strains = array([s[0][0]])  
        current_stresses = array([current_stress_list[0]])  # Considering only the relevant stress component

        strains.append(current_strains)
        stresses.append(current_stresses)

    strains = vstack(strains)
    stresses = vstack(stresses)

    cij = linalg.lstsq(strains, stresses, rcond=None)[0] * conversion_factor


    c11 = cij[0][0]  # Only one component in 1D with longitudinal strain

    length = pos_optimized.get_cell_lengths_and_angles()[2]

    elastic_constants_dict['c11'] = c11 #* length / 10.

    return elastic_constants_dict





# Example usage
#original_value = 100  # Replace this with the original value you want to decrease
#decrease_percent = 5  # Replace this with the desired percentage decrease
#decreased_value = percent_decrease(original_value, decrease_percent)
#print(f"The {decrease_percent} percent decreased value is: {decreased_value}")


import numpy as np
from numpy.linalg import lstsq


import numpy as np

import numpy as np

def NanotubeJax(pos_optimized, latt_system, elastic_constants_dict, stress_set_dict, conversion_factor):
    stresses = []
    strains = []

    # Loop over stress_set_dict to calculate strains and stresses for each key (strain type)
    for key in stress_set_dict.keys():
        s = strain_matrix(latt_system, key)[0]  # Getting the strain matrix
        current_stress_list = np.array(stress_set_dict[key])

        # Assuming the stress list is of size 3
        if current_stress_list.size != 3:
            raise ValueError("Error: Stress list should be of size 3.")
        
        current_stresses = current_stress_list.reshape(-1, 1)  # reshaping stress list to column vector

        # Select appropriate strain matrix based on key
        if key == "Longitudinal":
            current_strains = np.array([[s[2][2], 0, 0],
                                        [0, 0, 0],
                                        [0, 0, 0]])
        elif key == "Circumferential":
            current_strains = np.array([[s[0][0], 0, 0],
                                        [0, s[1][1], 0],
                                        [0, 0, 0]])
        elif key == "Torsional":
            current_strains = np.array([[0, 0, s[0][2]],
                                        [0, 0, 0],
                                        [s[2][0], 0, 0]])
        else:
            raise ValueError("Error: Unrecognized strain type.")

        strains.append(current_strains)
        stresses.append(current_stresses)

    # Stack the strains and stresses to form the matrices
    strains = np.hstack(strains)
    stresses = np.hstack(stresses)

    if strains.shape[0] != stresses.shape[0]:
        raise ValueError("Error: Incompatible dimensions for strains and stresses.")

    # Compute the elastic constants using the least squares method
    cij = np.linalg.lstsq(strains, stresses, rcond=None)[0]

    c33 = cij[0][0]  # Longitudinal elastic constant
    c22 = cij[1][1]  # Circumferential elastic constant
    c23 = cij[2][2]  # Torsional elastic constant

    # Obtain length and angle from pos_optimized
    length, angle = pos_optimized.get_cell_lengths_and_angles()[:2]

    # Store the elastic constants in the elastic_constants_dict
    elastic_constants_dict['c33'] = c33 * length * angle / 100.  # Longitudinal elastic constant
    elastic_constants_dict['c22'] = c22 * length * angle / 100.  # Circumferential elastic constant
    elastic_constants_dict['c23'] = c23 * length * angle / 100.  # Torsional elastic constant

    return elastic_constants_dict







import numpy as np

def NanotubeWatch(pos_optimized, latt_system, elastic_constants_dict, stress_set_dict, conversion_factor):
    stresses = []
    strains = []

    # Loop over stress_set_dict to calculate strains and stresses for each key (strain type)
    for key in stress_set_dict.keys():
        s = strain_matrix(latt_system, key)[0]  # Getting the strain matrix
        current_stress_list = np.array(stress_set_dict[key][0])

        # Considering deformation in the system
        # Fix for the shear strain component
        current_strains = np.array([[s[2][2], s[0][0], s[0][2] if s[0][2]!=0 else s[2][0]],  # Longitudinal, circumferential, and shear strains (εz, εr, γ)
                                    [s[0][0], s[1][1], s[1][2]],
                                    [s[0][2], s[1][2], s[2][2]]])

        current_stresses = np.array([[current_stress_list[0]],  # Longitudinal stress (σzz)
                                     [current_stress_list[1]],  # Circumferential stress (σrr)
                                     [current_stress_list[2]]])  # Shear stress (τrz)

        strains.append(current_strains)
        stresses.append(current_stresses)

    # Stack the strains and stresses to form the matrices
    strains = np.append(strains)
    stresses = np.append(stresses)

    # Compute the elastic constants using the least squares method
    cij = lstsq(strains, stresses, rcond=None)[0] * conversion_factor

    c33 = cij[0][0]  # Longitudinal elastic constant
    c22 = cij[1][1]  # Circumferential elastic constant
    c23 = cij[2][2]  # Shear (Torsional) elastic constant

    # Obtain length and angle from pos_optimized
    length, angle = pos_optimized.get_cell_lengths_and_angles()[:2]

    # Store the elastic constants in the elastic_constants_dict
    elastic_constants_dict['c33'] = c33 * length * angle / 100.  # Longitudinal elastic constant
    elastic_constants_dict['c22'] = c22 * length * angle / 100.  # Circumferential elastic constant
    elastic_constants_dict['c23'] = c23 * length * angle / 100.  # Torsional elastic constant

    return elastic_constants_dict



def Nanotube(pos_optimized, latt_system, elastic_constants_dict, stress_set_dict, conversion_factor):
    stresses = []
    strains = []
    for key in stress_set_dict.keys():
        strain_matrices = strain_matrix(latt_system, key)  # Getting the strain matrices
        current_stress_list = array(stress_set_dict[key][0])

        for s in strain_matrices:
            current_strains = array([[s[0][0], s[0][1], s[0][2]],  # longitudinal
                                     [s[1][0], s[1][1], s[1][2]],  # circumferential
                                     [s[2][0], s[2][1], s[2][2]]]) # torsional
            strains.append(current_strains)

        for i in range(3):
            current_stresses = array([[current_stress_list[0]],  # Longitudinal stress (σzz)
                                      [current_stress_list[1]],  # Circumferential stress (σrr)
                                      [current_stress_list[2]]]) # Shear stress (τrz)
            stresses.append(current_stresses)

    strains = vstack(strains)
    stresses = vstack(stresses)

    cij = linalg.lstsq(strains, stresses, rcond=None)[0] * conversion_factor
    #print(cij)
    #exit(1)

    c33 = cij[0][0]  # Longitudinal elastic constant
    c22 = cij[1][0]  # Circumferential elastic constant
    c23 = cij[2][0]  # Shear (Torsional) elastic constant

    length = pos_optimized.get_cell_lengths_and_angles()
    cross_sectional_area = pi*length[0]*length[1]/4.

    elastic_constants_dict['c33'] = c33 #* cross_sectional_area/length[2]  / 10.  # Longitudinal elastic constant
    elastic_constants_dict['c22'] = c22 #* cross_sectional_area/length[2]  / 10.  # Circumferential elastic constant
    elastic_constants_dict['c23'] = c23 #* cross_sectional_area/length[2]  / 10.  # Torsional elastic constant

    return elastic_constants_dict



def NanotubeDifferentdim(pos_optimized, latt_system, elastic_constants_dict, stress_set_dict, conversion_factor):
    stresses = []
    strains = []
    for key in stress_set_dict.keys():
        strain_matrices = strain_matrix(latt_system, key)  # Getting the strain matrices
        current_stress_list = array(stress_set_dict[key][0])

        for s in strain_matrices:
            current_strains = array([[s[0][0], s[0][1], s[0][2]],  # longitudinal
                                     [s[1][0], s[1][1], s[1][2]],  # circumferential
                                     [s[2][0], s[2][1], s[2][2]]]) # torsional

            strains.append(current_strains)

        current_stresses = array([[current_stress_list[0]],  # Longitudinal stress (σzz)
                                  [current_stress_list[1]],  # Circumferential stress (σrr)
                                  [current_stress_list[2]]]) # Shear stress (τrz)

        stresses.append(current_stresses)


    strains = vstack(strains)
    stresses = vstack(stresses)

    cij = linalg.lstsq(strains, stresses, rcond=None)[0] * conversion_factor
    
    print(strains)
    exit(1)

    c33 = cij[2][0]  # Longitudinal elastic constant
    c22 = cij[0][0]  # Circumferential elastic constant
    c23 = cij[1][0]  # Shear (Torsional) elastic constant

    length = pos_optimized.get_cell_lengths_and_angles()
    cross_sectional_area = pi*length[0]*length[1]/4.

    elastic_constants_dict['c33'] = c33 * cross_sectional_area/length[2]  / 10.  # Longitudinal elastic constant
    elastic_constants_dict['c22'] = c22 * cross_sectional_area/length[2]  / 10.  # Circumferential elastic constant
    elastic_constants_dict['c23'] = c23 * cross_sectional_area/length[2]  / 10.  # Torsional elastic constant

    return elastic_constants_dict



def Nanotube08042023(pos_optimized, latt_system, elastic_constants_dict, stress_set_dict, conversion_factor):
    # Initialize lists to store the strain and stress components for all scenarios
    strain_components = []
    stress_components = []

    # Loop over all strain scenarios
    for key in stress_set_dict.keys():
        # Compute the strain matrices for this scenario
        strain_matrices = strain_matrix(latt_system, key)
        
        # Loop over each strain matrix and the corresponding stress component
        for i, strain_matrix in enumerate(strain_matrices):
            # Append the strain component (i.e., the average strain for this scenario and direction) to the list
            strain_component = np.trace(strain_matrix) / 3
            strain_components.append(strain_component)
            
            # Append the stress component (i.e., the average stress for this scenario and direction) to the list
            stress_component = stress_set_dict[key][0][i]
            stress_components.append(stress_component)

    # Convert the lists to arrays and reshape them into column vectors
    #strains = np.array(strain_components).reshape(-1, 1)
    #stresses = np.array(stress_components).reshape(-1, 1)
    strains = vstack(strain_components)
    stresses = vstack(stress_components)

    # Solve for the elastic constants
    cij = np.linalg.lstsq(strains, stresses, rcond=None)[0] * conversion_factor

    print(cij)
    #exit(1)
    c33 = cij[2][0]  # Longitudinal elastic constant
    c22 = cij[0][0]  # Circumferential elastic constant
    c23 = cij[1][0]  # Shear (Torsional) elastic constant
    
    length = pos_optimized.get_cell_lengths_and_angles()
    cross_sectional_area = pi*length[0]*length[1]/4.

    elastic_constants_dict['c33'] = c33 * cross_sectional_area/length[2]  / 10.  # Longitudinal elastic constant
    elastic_constants_dict['c22'] = c22 * cross_sectional_area/length[2]  / 10.  # Circumferential elastic constant
    elastic_constants_dict['c23'] = c23 * cross_sectional_area/length[2]  / 10.  # Torsional elastic constant

    return elastic_constants_dict





def NanotubeThisisOkay(pos_optimized, latt_system, elastic_constants_dict, stress_set_dict, conversion_factor):
    stresses = []
    strains = []
    for key in stress_set_dict.keys():
        s = strain_matrix(latt_system, key)[0]  # Getting the strain matrix
        current_stress_list = array(stress_set_dict[key][0])

        current_strains = array([[s[0][0], s[0][1], s[0][2]],  # longitudinal
                                 [s[1][0], s[1][1], s[1][2]],  # circumferential
                                 [s[2][0], s[2][1], s[2][2]]]) # torsional


        current_stresses = array([[current_stress_list[0]],  # Longitudinal stress (σzz)
                                  [current_stress_list[1]],  # Circumferential stress (σrr)
                                  [current_stress_list[2]]]) # Shear stress (τrz)

        strains.append(current_strains)
        stresses.append(current_stresses)

    strains = vstack(strains)
    stresses = vstack(stresses)

    cij = linalg.lstsq(strains, stresses, rcond=None)[0] * conversion_factor
    print(cij)
    #exit(1)
    c33 = cij[2][0]  # Longitudinal elastic constant
    c22 = cij[0][0]  # Circumferential elastic constant
    c23 = cij[1][0]  # Shear (Torsional) elastic constant
    
    length = pos_optimized.get_cell_lengths_and_angles()
    cross_sectional_area = pi*length[0]*length[1]/4.

    elastic_constants_dict['c33'] = c33 * cross_sectional_area/length[2]  / 10.  # Longitudinal elastic constant
    elastic_constants_dict['c22'] = c22 * cross_sectional_area/length[2]  / 10.  # Circumferential elastic constant
    elastic_constants_dict['c23'] = c23 * cross_sectional_area/length[2]  / 10.  # Torsional elastic constant

    return elastic_constants_dict




def NanotubeModifying(pos_optimized, latt_system, elastic_constants_dict, stress_set_dict, conversion_factor):
    stresses = []
    strains = []

    # Loop over stress_set_dict to calculate strains and stresses for each key (strain type)
    for key in stress_set_dict.keys():
        s = strain_matrix(latt_system, key)[0]  # Getting the strain matrix

        current_stress_list = np.array(stress_set_dict[key][0])


        # Considering deformation in the system
 #       current_strains = np.array([[s[2][2], s[0][0], s[0][2]],  # Longitudinal, circumferential, and shear strains (εz, εr, γ)
 #                                   [s[0][0], s[1][1], s[1][2]],
 #                                   [s[0][2], s[1][2], s[2][2]]])

#        current_strains = np.array([[s[0][0], s[0][1]/2, s[0][2]/2],
#                                   [s[1][0]/2, s[1][1], s[1][2]/2],
#                                   [s[2][0]/2, s[2][1]/2, s[2][2]]])



#        current_strains = np.array([[s[2][2], s[0][0], s[0][2] if s[0][2]!=0 else s[2][0]],  # Longitudinal, circumferential, and shear strains (εz, εr, γ)
#                                    [s[0][0], s[1][1], s[1][2]],
#                                    [s[0][2], s[1][2], s[2][2]]])


#current_strains = np.array([[s_2[0][0], 0., s_3[0][2]],  # εr, 0, γrz
#                            [0., s_2[1][1], 0.],  # 0, εθ, 0
#                            [s_3[2][0], 0., s_1[2][2]]]) # γzr, 0, εz

#        current_strains = np.array([[s[0][0], 0., s[0][2]],  # εr, 0, γrz
#                                   [0., s[1][1], 0.],  # 0, εθ, 0
#                                   [s[2][0], 0., s[2][2]]])  # γzr, 0, εz




        current_strains = array([[s[0][0], s[0][1], s[0][2]],  # longitudinal
                                 [s[1][0], s[1][1], s[1][2]],  # circumferential
                                 [s[2][0], s[2][1], s[2][2]]]) # torsional



# Definition below:
#epsilon =	[er 0 yrz]
#		[0 e\theta 0]
#		[yzr 0 ez]	

	
#After stacking it, it will look like as shown below
#current_strains = np.array([
#    [c_up, 0., -t_up],  # εr, 0, γrz
#    [0., c_up, 0.],  # 0, εθ, 0
#    [t_up, 0., up]  # γzr, 0, εz
#])


#current_strains = np.array([s[2][2], s[0][0], s[0][2] if s[0][2]!=0 else s[2][0]])




#        current_stresses = np.array([[current_stress_list[2]],  # Longitudinal stress (σzz)
#                                     [current_stress_list[1]],  # Circumferential stress (σrr)
#                                     [current_stress_list[3]]])  # Shear stress (τrz)


        current_stresses = array([[current_stress_list[0]],  # Longitudinal stress (σzz)
                                     [current_stress_list[1]],  # Circumferential stress (σrr)
                                     [current_stress_list[2]]]) # Shear stress (τrz)


        strains.append(current_strains)
        stresses.append(current_stresses)

      
    # Stack the strains and stresses to form the matrices
    strains = vstack(strains)
    stresses = vstack(stresses)
   # print( array(strains))
   # exit(1)
    #print("stresses  ", stresses)
    #print("This is currentnt ", strains)
    #exit(1)
    # Compute the elastic constants using the least squares method
    cij = linalg.lstsq(strains, stresses, rcond=None)[0] * conversion_factor
    #print("DDDDDDDDDDDDD ", cij)
    #exit(1)

# Check finite ones
#c33 = cij[0][0]  # Stiffness related to longitudinal strain (?zz)
#c44 = cij[0][2]  # Stiffness related to shear strain (?rz)
#c11 = cij[1][1]  # Stiffness related to circumferential strain (???)
#c66 = cij[2][2]  # Stiffness related to torsional strain (???)

    c33 = cij[2][0]  # Longitudinal elastic constant
    c22 = cij[0][0]  # Circumferential elastic constant
    c23 = cij[1][0]  # Shear (Torsional) elastic constant
    

    # get cell lengths and angles
    length = pos_optimized.get_cell_lengths_and_angles()
    cross_sectional_area = pi*length[0]*length[1]/4.

    # Store the elastic constants in the elastic_constants_dict
    elastic_constants_dict['c33'] = c33 * cross_sectional_area/length[2]  / 100.  # Longitudinal elastic constant
    elastic_constants_dict['c22'] = c22 * cross_sectional_area/length[2]  / 100.  # Circumferential elastic constant
    elastic_constants_dict['c23'] = c23 * cross_sectional_area/length[2]  / 100.  # Torsional elastic constant

    return elastic_constants_dict


def Nanotubeworked(pos_optimized, latt_system, elastic_constants_dict, stress_set_dict, conversion_factor):
    stresses = []
    strains = []

    # Loop over stress_set_dict to calculate strains and stresses for each key (strain type)
    for key in stress_set_dict.keys():
        s = strain_matrix(latt_system, key)[0]  # Getting the strain matrix
        current_stress_list = np.array(stress_set_dict[key][0])

        # Considering deformation in the system
        current_strains = np.array([[s[1][1], s[2][2]],  # Strains from the strain matrix
                                    [s[2][2], s[1][1]]])

        current_stresses = np.array([[current_stress_list[1]],  # Stresses from the stress list
                                     [current_stress_list[2]]])

        strains.append(current_strains)
        stresses.append(current_stresses)

    # Stack the strains and stresses to form the matrices
    strains = np.vstack(strains)
    stresses = np.vstack(stresses)

    # Compute the elastic constants using the least squares method
    cij = lstsq(strains, stresses, rcond=None)[0] * conversion_factor

    c22 = cij[0][0]  # Circumferential elastic constant
    c33 = cij[1][1]  # Circumferential elastic constant
    c23 = cij[1][0]  # Torsional elastic constant

    # Compute c33 (Longitudinal elastic constant) using E = c22 - c23*c23/c22
    #c33 = c22 - c23 * c23 / c22

    # Obtain length and angle from pos_optimized
    length, angle = pos_optimized.get_cell_lengths_and_angles()[:2]

    # Store the elastic constants in the elastic_constants_dict
    elastic_constants_dict['c22'] = c22 * length * angle / 100.  # Circumferential elastic constant
    elastic_constants_dict['c23'] = c23 * length * angle / 100.  # Torsional elastic constant
    elastic_constants_dict['c33'] = c33 * length * angle / 100.  # Longitudinal elastic constant

    return elastic_constants_dict





def NanotubeOld(pos_optimized, latt_system, elastic_constants_dict, stress_set_dict, conversion_factor):
    stresses = []
    strains = []

    for key in stress_set_dict.keys():
        s = strain_matrix(latt_system, key)  # Getting the strain matrix
        strains.append(s)
        current_stress_list = array(stress_set_dict[key][0])

        # Considering deformation in the system
#        current_strains = array([s[2][2], (s[0][0] + s[1][1])/2, s[0][2]])  # Longitudinal, circumferential, and torsional strains
#        current_stresses = array([current_stress_list[2], current_stress_list[1], current_stress_list[4]])  # Relevant stress components


        current_stresses = array([current_stress_list[2], current_stress_list[1], current_stress_list[4]])  # Relevant stress components
        stresses.append(current_stresses)


#        strains.append(current_strains)
#        stresses.append(current_stresses)

#    strains = vstack(strains)
#    stresses = vstack(stresses)


    strains = vstack(strains)
    stresses = vstack(stresses)



    cij = lstsq(strains, stresses, rcond=None)[0] * conversion_factor

    c33 = cij[0] #[0]  # Longitudinal elastic constant
    c22 = cij[1] #[1]  # Circumferential elastic constant
    c23 = cij[2] #[2]  # Torsional elastic constant

    length = pos_optimized.get_cell_lengths_and_angles()[2]

    elastic_constants_dict['c33'] = c33 * length / 100.
    elastic_constants_dict['c22'] = c22 * length / 100.
    elastic_constants_dict['c23'] = c23 * length / 100.

    return elastic_constants_dict



def calc_elastic_constants(pos_optimized, latt_system, elastic_constants_dict, stress_set_dict):
    # kB to GPa
    conversion_factor = -0.1
    strains_matrix = indict['strains_matrix'][0]
    if strains_matrix != 'asess':
        if latt_system == 'Cubic':
            elastic_constants_dict = Cubic(latt_system, elastic_constants_dict, stress_set_dict, conversion_factor)

        elif latt_system == 'Hexagonal':
            elastic_constants_dict = Hexagonal(latt_system, elastic_constants_dict, stress_set_dict, conversion_factor)
                
        elif latt_system == 'Trigonal1':
            elastic_constants_dict = Trigonal1(latt_system, elastic_constants_dict, stress_set_dict, conversion_factor)

        elif latt_system == 'Trigonal2':
            elastic_constants_dict = Trigonal2(latt_system, elastic_constants_dict, stress_set_dict, conversion_factor)

        elif latt_system == 'Tetragonal1':
            elastic_constants_dict = Tetragonal1(latt_system, elastic_constants_dict, stress_set_dict, conversion_factor)

        elif latt_system == 'Tetragonal2':
            elastic_constants_dict = Tetragonal2(latt_system, elastic_constants_dict, stress_set_dict, conversion_factor)

        elif latt_system == 'Orthorombic':
            elastic_constants_dict = Orthorombic(latt_system, elastic_constants_dict, stress_set_dict, conversion_factor)

        elif latt_system == 'Monoclinic':
            elastic_constants_dict = Monoclinic(latt_system, elastic_constants_dict, stress_set_dict, conversion_factor)
                
        elif latt_system == 'Triclinic':
            elastic_constants_dict = Triclinic(latt_system, elastic_constants_dict, stress_set_dict, conversion_factor)
        
        elif latt_system == 'isotropy':
            elastic_constants_dict = isotropy(pos_optimized, latt_system, elastic_constants_dict, stress_set_dict, conversion_factor)

        elif latt_system == 'tetragonal':
            elastic_constants_dict = tetragonal(pos_optimized, latt_system, elastic_constants_dict, stress_set_dict, conversion_factor)

        elif latt_system == 'orthotropy':
            elastic_constants_dict = orthotropy(pos_optimized, latt_system, elastic_constants_dict, stress_set_dict, conversion_factor)

        elif latt_system == 'anisotropy':
            elastic_constants_dict = anisotropy(pos_optimized, latt_system, elastic_constants_dict, stress_set_dict, conversion_factor)                        

        elif latt_system == 'any1D':
            elastic_constants_dict = any1D(pos_optimized, latt_system, elastic_constants_dict, stress_set_dict, conversion_factor)  
        elif latt_system == 'true1D':
            elastic_constants_dict = true1D(pos_optimized, latt_system, elastic_constants_dict, stress_set_dict, conversion_factor)    
        elif latt_system == 'Nanotube':
            elastic_constants_dict = Nanotube(pos_optimized, latt_system, elastic_constants_dict, stress_set_dict, conversion_factor)                      
        else:
            print('Crystal system is not parsed correctly!!!')
            exit(1)

    elif strains_matrix == 'asess':
        if indict['dimensional'][0] == '3D':
            elastic_constants_dict = Triclinic(latt_system, elastic_constants_dict, stress_set_dict, conversion_factor)
        elif indict['dimensional'][0] == '2D':
            elastic_constants_dict = anisotropy(pos_optimized, latt_system, elastic_constants_dict, stress_set_dict, conversion_factor)                        
        elif indict['dimensional'][0] == '1D':
            elastic_constants_dict = any1D(pos_optimized, latt_system, elastic_constants_dict, stress_set_dict, conversion_factor)                        

    return elastic_constants_dict
