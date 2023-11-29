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
from numpy import array
from numpy.linalg import det

def criteria(elastic_constants_dict, latt_system):
    stable = False
    condition1 = True
    condition2 = True
    condition3 = True
    condition4 = True
    condition5 = True
    condition6 = True
    condition7 = True
    #for cij in elastic_constants_dict.keys():
    #    globals()[cij] = elastic_constants_dict[cij]
    if 'c11' in elastic_constants_dict:
        c11 = elastic_constants_dict['c11']
    if 'c12' in elastic_constants_dict:
        c12 = elastic_constants_dict['c12']
    if 'c13' in elastic_constants_dict:
        c13 = elastic_constants_dict['c13']
    if 'c14' in elastic_constants_dict:
        c14 = elastic_constants_dict['c14']
    if 'c15' in elastic_constants_dict:
        c15 = elastic_constants_dict['c15']
    if 'c16' in elastic_constants_dict:
        c16 = elastic_constants_dict['c16']
    if 'c22' in elastic_constants_dict:
        c22 = elastic_constants_dict['c22']
    if 'c23' in elastic_constants_dict:
        c23 = elastic_constants_dict['c23']
    if 'c24' in elastic_constants_dict:
        c24 = elastic_constants_dict['c24']
    if 'c25' in elastic_constants_dict:
        c25 = elastic_constants_dict['c25']
    if 'c26' in elastic_constants_dict:
        c26 = elastic_constants_dict['c26']
    if 'c33' in elastic_constants_dict:
        c33 = elastic_constants_dict['c33']
    if 'c34' in elastic_constants_dict:
        c34 = elastic_constants_dict['c34']
    if 'c35' in elastic_constants_dict:
        c35 = elastic_constants_dict['c35']
    if 'c36' in elastic_constants_dict:
        c36 = elastic_constants_dict['c36']
    if 'c44' in elastic_constants_dict:
        c44 = elastic_constants_dict['c44']
    if 'c45' in elastic_constants_dict:
        c45 = elastic_constants_dict['c45']
    if 'c46' in elastic_constants_dict:
        c46 = elastic_constants_dict['c46']
    if 'c55' in elastic_constants_dict:
        c55 = elastic_constants_dict['c55']
    if 'c56' in elastic_constants_dict:
        c56 = elastic_constants_dict['c56']
    if 'c66' in elastic_constants_dict:
        c66 = elastic_constants_dict['c66']

    # PHYSICAL REVIEW B 90, 224104 (2014)
    if latt_system == 'Cubic':
        condition1 = c11 - c12 > 0
        condition2 = c11 + 2*c12 > 0
        condition3 = c44 > 0

    elif latt_system == 'Hexagonal' or latt_system == 'Tetragonal1':
        condition1 = c11 > abs(c12)
        condition2 = 2*c13*c13 < c33*(c11+c12)
        condition3 = c44 > 0
        condition4 = c11 > c12

    elif latt_system == 'Trigonal1':
        condition1 = c11 > abs(c12)
        condition2 = 2*c13*c13 < c33*(c11+c12)
        condition3 = c44 > 0
        condition4 = 2*c14*c14 < c44*(c11-c12)

    elif latt_system == 'Trigonal2':
        condition1 = c11 > abs(c12)
        condition2 = 2*c13*c13 < c33*(c11+c12)
        condition3 = c44 > 0
        condition4 = 2*c14*c14 + 2*c15*c15 < c44*(c11-c12)

    elif latt_system == 'Tetragonal2':
        condition1 = c11 > abs(c12)
        condition2 = 2*c13*c13 < c33*(c11+c12)
        condition3 = c44 > 0
        condition4 = 2*c16*c16 < c66*(c11-c12)

    elif latt_system == 'Orthorombic':
        condition1 = c11 > 0
        condition2 = c11 * c22 > c12 * c12
        condition3 = c11*c22*c33 + 2*c12*c13*c23 - c11*c23*c23 - c22*c13*c13 - c33*c12*c12 > 0
        condition4 = c44 > 0 and c55 > 0 and c66 > 0

    elif latt_system == 'Monoclinic':
        # PHYSICAL REVIEW B 76, 054115 (2007)
        g = c11*c22*c33-c11*c23*c23-c22*c13*c13-c33*c12*c12+2*c12*c13*c23
        condition1 = c11 > 0 and c22 > 0 and c33 > 0 and c44 > 0 and c55 > 0 and c66 > 0
        condition2 = c11 + c22 + c33 + 2*(c12+c13+c23) > 0
        condition3 = c33*c55 > c35*c35
        condition4 = c44*c66 > c46*c46
        condition5 = c22 + c33 - 2*c23 > 0
        condition6 = c22*(c33*c55-c35*c35)+2*c23*c25*c35-c23*c23*c55-c25*c25*c33 > 0
        condition7 = 2*(c15*c25*(c33*c12-c13*c23)+c15*c35*(c22*c13-c12*c23)+c25*c35*(c11*c23-c12*c13)) - \
                     (c15*c15*(c22*c33-c23*c23)+c25*c25*(c11*c33-c13*c13)+c35*c35*(c11*c22-c12*c12)) + c55*g > 0
            
    elif latt_system == 'Triclinic':
        condition1 = det(array([c11, c12, c13, c14, c15, c16],
                                  [c12, c22, c23, c24, c25, c26],
                                  [c13, c23, c33, c34, c35, c36],
                                  [c14, c24, c34, c44, c45, c46],
                                  [c15, c25, c35, c45, c55, c56],
                                  [c16, c26, c36, c46, c56, c66])) > 0

    # 2D Mater. 6 (2019) 048001
    elif latt_system == 'isotropy':
        condition1 = c11 > 0
        condition2 = c11 > abs(c12)

    elif latt_system == 'tetragonal':
        condition1 = c11 > 0
        condition2 = c66 > 0
        condition3 = c11 > abs(c12)

    elif latt_system == 'orthotropy':
        condition1 = c11 > 0
        condition2 = c66 > 0
        condition3 = c11*c22 > c12*c12

    elif latt_system == 'anisotropy':
        condition1 = c11 > 0
        condition2 = c11*c22 > c12*c12
        condition3 = det(array([[c11, c12, c16],
                        [c12, c22, c26],
                        [c16, c26, c66]])) > 0

#1D crystal
    elif latt_system == 'Nanotube':
        condition1 = elastic_constants_dict['c33'] > 0 
        condition2 = elastic_constants_dict['c33'] > elastic_constants_dict['c23']

    else:
        print('Crystal system is not parsed correctly!!!')
        exit(1)

    stable = condition1 and condition2 and condition3 and condition4 and condition5 and condition6 and condition7

    return stable
