"""
  Elastool -- Elastic toolkit for zero and finite-temperature elastic constants and mechanical properties calculations

  Copyright (C) 2019-2024 by Zhong-Li Liu and Chinedu Ekuma

  This program is free software; you can redistribute it and/or modify it under the
  terms of the GNU General Public License as published by the Free Software Foundation
  version 3 of the License.

  This program is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
  PARTICULAR PURPOSE.  See the GNU General Public License for more details.

  E-mail: zl.liu@163.com or cekuma1@gmail.com
"""
from read_input import indict
from numpy import array, cos, radians, pi,sin

p_decr = lambda original_value, decrease_percent: original_value * (1 - decrease_percent / 100)

def strain_matrix(latt_system, up0):
    strain_matrix_list = []
    if indict['strains_matrix'][0] == 'ohess':
        up = up0
        if indict['dimensional'][0] == '3D':
            if latt_system == 'Cubic':
                strain_matrix_1 = array([[up, 0., 0.],
                                         [0., 0., up / 2.],
                                         [0., up / 2., 0.]])
                strain_matrix_list.append(strain_matrix_1)

            elif latt_system == 'Hexagonal':
                strain_matrix_1 = array([[up, 0., 0.],
                                         [0., 0., 0.],
                                         [0., 0., 0.]])
                strain_matrix_2 = array([[0., 0., 0.],
                                         [0., 0., up / 2.],
                                         [0., up / 2., up]])
                strain_matrix_list.append(strain_matrix_1)
                strain_matrix_list.append(strain_matrix_2)
                
            elif latt_system == 'Nanoribbon':
                strain_matrix_1 = array([[0., 0., 0.],
                         [0., up, 0.],
                         [0., 0., 0.]])
                strain_matrix_2 = array([[0., 0., up / 2.],
                         [0., 0., 0.],
                         [up / 2., 0., up]])

                strain_matrix_list.append(strain_matrix_1)
                strain_matrix_list.append(strain_matrix_2)
                   
                         
            elif latt_system == 'Trigonal1' or latt_system == 'Trigonal2':
                strain_matrix_1 = array([[up, 0., 0.],
                                         [0., 0., 0.],
                                         [0., 0., 0.]])
                strain_matrix_2 = array([[0., 0., 0.],
                                         [0., 0., up / 2.],
                                         [0., up / 2., up]])
                strain_matrix_list.append(strain_matrix_1)
                strain_matrix_list.append(strain_matrix_2)

            elif latt_system == 'Tetragonal1' or latt_system == 'Tetragonal2':
                strain_matrix_1 = array([[up, 0., 0.],
                                         [0., 0., 0.],
                                         [0., 0., 0.]])
                strain_matrix_2 = array([[0., up / 2., 0.],
                                         [up / 2., 0., up / 2.],
                                         [0., up / 2., up]])
                strain_matrix_list.append(strain_matrix_1)
                strain_matrix_list.append(strain_matrix_2)

            elif latt_system == 'Orthorombic':
                strain_matrix_1 = array([[up, 0., 0.],
                                         [0., 0., 0.],
                                         [0., 0., 0.]])
                strain_matrix_2 = array([[0., 0., 0.],
                                         [0., up, 0.],
                                         [0., 0., 0.]])
                strain_matrix_3 = array([[0., up / 2., up / 2.],
                                         [up / 2., 0., up / 2.],
                                         [up / 2., up / 2., up]])
                strain_matrix_list.append(strain_matrix_1)
                strain_matrix_list.append(strain_matrix_2)
                strain_matrix_list.append(strain_matrix_3)

            elif latt_system == 'Monoclinic':
                strain_matrix_1 = array([[up, 0., 0.],
                                         [0., 0., 0.],
                                         [0., 0., 0.]])
                strain_matrix_2 = array([[0., 0., 0.],
                                         [0., up, 0.],
                                         [0., 0., 0.]])
                strain_matrix_3 = array([[0., 0., 0.],
                                         [0., 0., up / 2.],
                                         [0., up / 2., up]])
                strain_matrix_4 = array([[0., up / 2., up / 2.],
                                         [up / 2., 0., 0.],
                                         [up / 2., 0., 0.]])
                strain_matrix_list.append(strain_matrix_1)
                strain_matrix_list.append(strain_matrix_2)
                strain_matrix_list.append(strain_matrix_3)
                strain_matrix_list.append(strain_matrix_4)

            elif latt_system == 'Triclinic':
                strain_matrix_1 = array([[up, 0., 0.],
                                         [0., 0., 0.],
                                         [0., 0., 0.]])
                strain_matrix_2 = array([[0., 0., 0.],
                                         [0., up, 0.],
                                         [0., 0., 0.]])
                strain_matrix_3 = array([[0., 0., 0.],
                                         [0., 0., 0.],
                                         [0., 0., up]])
                strain_matrix_4 = array([[0., 0., 0.],
                                         [0., 0., up / 2.],
                                         [0., up / 2., 0.]])
                strain_matrix_5 = array([[0., 0., up / 2.],
                                         [0., 0., 0.],
                                         [up / 2., 0., 0.]])
                strain_matrix_6 = array([[0., up / 2., 0.],
                                         [up / 2., 0., 0.],
                                         [0., 0., 0.]])
                strain_matrix_list.append(strain_matrix_1)
                strain_matrix_list.append(strain_matrix_2)
                strain_matrix_list.append(strain_matrix_3)
                strain_matrix_list.append(strain_matrix_4)
                strain_matrix_list.append(strain_matrix_5)
                strain_matrix_list.append(strain_matrix_6)

        elif indict['dimensional'][0] == '2D':
            # 2D in-plane strains: in xy plane
            if latt_system == 'isotropy':
                strain_matrix_1 = array([[up, 0., 0.],
                                         [0., 0., 0.],
                                         [0., 0., 0.]])
                strain_matrix_list.append(strain_matrix_1)

            elif latt_system == 'tetragonal':
                strain_matrix_1 = array([[up, up / 2, 0.],
                                         [up / 2, 0., 0.],
                                         [0., 0., 0.]])
                strain_matrix_list.append(strain_matrix_1)

            elif latt_system == 'orthotropy':
                strain_matrix_1 = array([[up, up / 2, 0.],
                                         [up / 2, 0., 0.],
                                         [0., 0., 0.]])
                strain_matrix_2 = array([[0., 0., 0.],
                                         [0., up, 0.],
                                         [0., 0., 0.]])
                strain_matrix_list.append(strain_matrix_1)
                strain_matrix_list.append(strain_matrix_2)

            elif latt_system == 'anisotropy':
                strain_matrix_1 = array([[up, 0., 0.],
                                         [0., 0., 0.],
                                         [0., 0., 0.]])
                strain_matrix_2 = array([[0., 0., 0.],
                                         [0., up, 0.],
                                         [0., 0., 0.]])
                strain_matrix_3 = array([[0., up / 2, 0.],
                                         [up / 2, 0., 0.],
                                         [0., 0., 0.]])
                strain_matrix_list.append(strain_matrix_1)
                strain_matrix_list.append(strain_matrix_2)
                strain_matrix_list.append(strain_matrix_3)

        elif indict['dimensional'][0] == '1D':
            if latt_system == 'Nanotube':
                strain_matrix_1 = array([[0, 0., 0.],
                                         [0., 0., up/2], 
                                         [0., up/2, up]])
                strain_matrix_list.append(strain_matrix_1)

    elif indict['strains_matrix'][0] == 'asess':
        up = up0
        if indict['dimensional'][0] == '3D':
            strain_matrix_1 = up0 * array([[1., 0., 0.],
                                           [0., 0., 0.],
                                           [0., 0., 0.]])
            strain_matrix_2 = up0 * array([[0., 0., 0.],
                                           [0., 1., 0.],
                                           [0., 0., 0.]])
            strain_matrix_3 = up0 * array([[0., 0., 0.],
                                           [0., 0., 0.],
                                           [0., 0., 1.]])
            strain_matrix_4 = up0 * array([[0., 0, 0.],
                                           [0., 0., 1. / 2.],
                                           [0., 1. / 2., 0.]])
            strain_matrix_5 = up0 * array([[0., 0., 1. / 2.],
                                           [0., 0., 0.],
                                           [1. / 2., 0., 0.]])
            strain_matrix_6 = up0 * array([[0., 1. / 2., 0.],
                                           [1. / 2., 0., 0.],
                                           [0., 0., 0.]])
            strain_matrix_list = [
                strain_matrix_1,
                strain_matrix_2,
                strain_matrix_3,
                strain_matrix_4,
                strain_matrix_5,
                strain_matrix_6]

        elif indict['dimensional'][0] == '2D':
            strain_matrix_1 = up0 * array([[1., 0., 0.],
                                           [0., 0., 0.],
                                           [0., 0., 0.]])
            strain_matrix_2 = up0 * array([[0., 0., 0.],
                                           [0., 1., 0.],
                                           [0., 0., 0.]])
            strain_matrix_3 = up0 * array([[0., 1. / 2., 0.],
                                           [1. / 2., 0., 0.],
                                           [0., 0., 0.]])
            strain_matrix_list = [
                strain_matrix_1,
                strain_matrix_2,
                strain_matrix_3]

        elif indict['dimensional'][0] == '1D':
            if latt_system == 'Nanotube':
                strain_matrix_1 = array([[0, 0., 0.],
                                         [0., up/2, 0.], 
                                         [0., 0., up]])
                strain_matrix_list.append(strain_matrix_1)



    elif indict['strains_matrix'][0] == 'ulics':
        up = up0
        if indict['dimensional'][0] == '3D':
            #up = 10. ** (-3.)
            up = up0 / 6.
            strain_matrix_1 = up * array([[1., 6. / 2., 5. / 2.],
                                          [6. / 2., 2., 4. / 2.],
                                          [5. / 2., 4. / 2., 3.]])
            strain_matrix_2 = up * array([[-2., -5. / 2., 6. / 2.],
                                          [-5. / 2., 1., -3. / 2.],
                                          [6. / 2., -3. / 2., 4.]])
            strain_matrix_3 = up * array([[3., -4. / 2., 2. / 2.],
                                          [-4. / 2., -5., 6. / 2.],
                                          [2. / 2., 6. / 2., -1.]])
            strain_matrix_4 = up * array([[-4., -2. / 2., -3. / 2.],
                                          [-2. / 2., -6., 1. / 2.],
                                          [-3. / 2., 1. / 2., 5.]])
            strain_matrix_5 = up * array([[5., -3. / 2., -1. / 2.],
                                          [-3. / 2., 4., -2. / 2.],
                                          [-1. / 2., -2. / 2., 6.]])
            strain_matrix_6 = up * array([[-6., 1. / 2., -4. / 2.],
                                          [1. / 2., 3., 5. / 2.],
                                          [-4. / 2., 5. / 2., -2.]])
            strain_matrix_list = [
                strain_matrix_1,
                strain_matrix_2,
                strain_matrix_3,
                strain_matrix_4,
                strain_matrix_5,
                strain_matrix_6]

        elif indict['dimensional'][0] == '2D':
            up = up0 / 3.
            strain_matrix_1 = up * array([[1., 3. / 2., 0.],
                                          [3. / 2., 2., 0.],
                                          [0., 0., 0.]])
            strain_matrix_2 = up * array([[2., 1. / 2., 0.],
                                          [1. / 2, 3., 0.],
                                          [0., 0., 0.]])
            strain_matrix_3 = up * array([[3., 1., 0.],
                                          [1., 1., 0.],
                                          [0., 0., 0.]])
            strain_matrix_list = [
                strain_matrix_1,
                strain_matrix_2,
                strain_matrix_3]

        elif indict['dimensional'][0] == '1D':
            if latt_system == 'Nanotube':
                strain_matrix_1 = array([[0, 0., 0.],
                                         [0., up/2, 0.], 
                                         [0., 0., up]])
                strain_matrix_list.append(strain_matrix_1)

    return strain_matrix_list
