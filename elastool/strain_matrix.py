"""
  Elastool -- Elastic toolkit for zero and finite-temperature elastic constants and mechanical properties calculations

  Copyright (C) 2019-2024 by Zhong-Li Liu and Chinedu Ekuma

  This program is free software; you can redistribute it and/or modify it under the
  terms of the GNU General Public License as published by the Free Software Foundation
  version 3 of the License.

  This program is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
  PARTICULAR PURPOSE.  See the GNU General Public License for more details.


Some notes for the 1D case
For any1D case actually defines a 3x3 strain matrix, just like the 2D case. However, it only has non-zero entries in the diagonal, which suggests that strain is being applied independently in the three dimensions. This matrix definition might make sense in the context of a 1D system embedded in a 3D space, where the non-zero entries denote strains in directions orthogonal to the 1D system. Here, up indicates the degree of deformation along the y-axis and 2*up along the z-axis. The 1D system is not strained along its length (x-axis), which is why the corresponding entry is zero.


Contrarily, the true1D strain matrix I initially suggested only considers deformation along the 1D system itself (without considering a possible embedding in higher-dimensional space). The strain is a scalar, i.e., a 1x1 matrix, since there is only one dimension in which the system can be deformed.
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
            if latt_system == 'any1D':
                strain_matrix_1 = array([[0, 0., 0.],
                                         [0., -0.05*up, 0.], #Add 1% of'up" as shearing strain
                                         [0., 0., up]])
                strain_matrix_list.append(strain_matrix_1)

            elif latt_system == 'true1D':
                strain_matrix_1 = up * array([[1.]])
                strain_matrix_2 = up * array([[2.]])
                strain_matrix_3 = up * array([[3.]])
                strain_matrix_list = [
                    strain_matrix_1, strain_matrix_2, strain_matrix_3]
            elif latt_system == 'Nanotube':
                strain_matrix_1 = array([[0., 0., 0.],
                                         [0., -0.05*up, 0.],
                                         [0., 0., up]])  # Longitudinal strain

#angle = pos_optimized.get_cell_lengths_and_angles()[2]

                c_up = up #p_decr(up,5)

		#c_tt = -c_rr * (1 + epsilon)
		#c_rt = c_rr * gamma
                strain_matrix_2 = array([[c_up, 0. , 0.], # εr, 0, 0 
                                         [0., c_up*(1.+0.01), 0.], # 0, εθ, 0
                                         [0., 0., -2 * c_up * cos(radians(60)) ** 2]])  # 0, 0, γrz  # Circumferential strain # Circumferential strain # -2 * c_up * cos(radians(0)) ** 2



                t_up = up*pi / 6.  # For example, 30 degrees twist #p_decr(up,10)

            #    strain_matrix_3 = array([[0., -t_up , 0.],
            #                             [t_up, 0., 0.],
            #                             [0., 0., 0.]])  # Torsional strain

#assuming the twist is happening about the z-axis (the long axis of the nanotube)
                strain_matrix_3 = array([[0., 0., -t_up],
                                        [0., 0., 0.],
                                        [t_up,0., 0.]])  # Torsional strain

                strain_matrix_list.append(strain_matrix_1)
                strain_matrix_list.append(strain_matrix_2)
                strain_matrix_list.append(strain_matrix_3)

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
            if latt_system == 'any1D':
                strain_matrix_1 = array([[0, 0., 0.],
                                         [0., -0.05*up, 0.], #Add 1% of'up" as shearing strain
                                         [0., 0., up]])
                strain_matrix_list.append(strain_matrix_1)

            elif latt_system == 'true1D':
                strain_matrix_1 = up * array([[1.]])
                strain_matrix_2 = up * array([[2.]])
                strain_matrix_3 = up * array([[3.]])
                strain_matrix_list = [
                    strain_matrix_1, strain_matrix_2, strain_matrix_3]
            elif latt_system == 'Nanotube':
                strain_matrix_1 = array([[0., 0., 0.],
                                         [0., -0.05*up, 0.],
                                         [0., 0., up]])  # Longitudinal strain

#angle = pos_optimized.get_cell_lengths_and_angles()[2]

                c_up = up #p_decr(up,5)

		#c_tt = -c_rr * (1 + epsilon)
		#c_rt = c_rr * gamma
                strain_matrix_2 = array([[c_up, 0. , 0.], # εr, 0, 0 
                                         [0., c_up*(1.+0.01), 0.], # 0, εθ, 0
                                         [0., 0., -2 * c_up * cos(radians(60)) ** 2]])  # 0, 0, γrz  # Circumferential strain # Circumferential strain # -2 * c_up * cos(radians(0)) ** 2



                t_up = up*pi / 6.  # For example, 30 degrees twist #p_decr(up,10)

            #    strain_matrix_3 = array([[0., -t_up , 0.],
            #                             [t_up, 0., 0.],
            #                             [0., 0., 0.]])  # Torsional strain

#assuming the twist is happening about the z-axis (the long axis of the nanotube)
                strain_matrix_3 = array([[0., 0., -t_up],
                                        [0., 0., 0.],
                                        [t_up,0., 0.]])  # Torsional strain

                strain_matrix_list.append(strain_matrix_1)
                strain_matrix_list.append(strain_matrix_2)
                strain_matrix_list.append(strain_matrix_3)
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
            if latt_system == 'any1D':
                strain_matrix_1 = array([[0, 0., 0.],
                                         [0., -0.05*up, 0.], #Add 1% of'up" as shearing strain
                                         [0., 0., up]])
                strain_matrix_list.append(strain_matrix_1)

            elif latt_system == 'true1D':
                strain_matrix_1 = up * array([[1.]])
                strain_matrix_2 = up * array([[2.]])
                strain_matrix_3 = up * array([[3.]])
                strain_matrix_list = [
                    strain_matrix_1, strain_matrix_2, strain_matrix_3]
            elif latt_system == 'Nanotube':
                strain_matrix_1 = array([[0., 0., 0.],
                                         [0., -0.05*up, 0.],
                                         [0., 0., up]])  # Longitudinal strain

#angle = pos_optimized.get_cell_lengths_and_angles()[2]

                c_up = up #p_decr(up,5)

		#c_tt = -c_rr * (1 + epsilon)
		#c_rt = c_rr * gamma
                strain_matrix_2 = array([[c_up, 0. , 0.], # εr, 0, 0 
                                         [0., c_up*(1.+0.01), 0.], # 0, εθ, 0
                                         [0., 0., -2 * c_up * cos(radians(60)) ** 2]])  # 0, 0, γrz  # Circumferential strain # Circumferential strain # -2 * c_up * cos(radians(0)) ** 2



                t_up = up*pi / 6.  # For example, 30 degrees twist #p_decr(up,10)

            #    strain_matrix_3 = array([[0., -t_up , 0.],
            #                             [t_up, 0., 0.],
            #                             [0., 0., 0.]])  # Torsional strain

#assuming the twist is happening about the z-axis (the long axis of the nanotube)
                strain_matrix_3 = array([[0., 0., -t_up],
                                        [0., 0., 0.],
                                        [t_up,0., 0.]])  # Torsional strain

                strain_matrix_list.append(strain_matrix_1)
                strain_matrix_list.append(strain_matrix_2)
                strain_matrix_list.append(strain_matrix_3)

    return strain_matrix_list
