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
from numpy import identity, dot
from read_input import indict
from strain_matrix import strain_matrix


def deform_cell_ulics(latt_system, cell, up):
    deformed_cell_list = []
    deform_matrix_list = []
    identity_matrix = identity(3)
    if indict['dimensional'][0] == '3D':
        deform_matrix_1 = identity_matrix + strain_matrix(latt_system, up)[0]
        deform_matrix_2 = identity_matrix + strain_matrix(latt_system, up)[1]
        deform_matrix_3 = identity_matrix + strain_matrix(latt_system, up)[2]
        deform_matrix_4 = identity_matrix + strain_matrix(latt_system, up)[3]
        deform_matrix_5 = identity_matrix + strain_matrix(latt_system, up)[4]
        deform_matrix_6 = identity_matrix + strain_matrix(latt_system, up)[5]
                                                                                                                    
        if latt_system == 'Cubic':
            deformed_cell = dot(cell, deform_matrix_1)
            deformed_cell_list.append(deformed_cell)

        elif latt_system == 'Hexagonal':
            deform_matrix_list.append(deform_matrix_1)
            deform_matrix_list.append(deform_matrix_2)
            for deform_matrix in deform_matrix_list:
                deformed_cell = dot(cell, deform_matrix)
                deformed_cell_list.append(deformed_cell)
                
        elif latt_system == 'Trigonal1':
            deform_matrix_list.append(deform_matrix_1)
            deform_matrix_list.append(deform_matrix_2)
            for deform_matrix in deform_matrix_list:
                deformed_cell = dot(cell, deform_matrix)
                deformed_cell_list.append(deformed_cell)

        elif latt_system == 'Trigonal2':
            deform_matrix_list.append(deform_matrix_1)
            deform_matrix_list.append(deform_matrix_2)
            for deform_matrix in deform_matrix_list:
                deformed_cell = dot(cell, deform_matrix)
                deformed_cell_list.append(deformed_cell)

        elif latt_system == 'Tetragonal1':
            deform_matrix_list.append(deform_matrix_1)
            deform_matrix_list.append(deform_matrix_2)
            for deform_matrix in deform_matrix_list:
                deformed_cell = dot(cell, deform_matrix)
                deformed_cell_list.append(deformed_cell)

        elif latt_system == 'Tetragonal2':
            deform_matrix_list.append(deform_matrix_1)
            deform_matrix_list.append(deform_matrix_2)
            for deform_matrix in deform_matrix_list:
                deformed_cell = dot(cell, deform_matrix)
                deformed_cell_list.append(deformed_cell)

        elif latt_system == 'Orthorombic':
            deform_matrix_list.append(deform_matrix_1)
            deform_matrix_list.append(deform_matrix_2)
            deform_matrix_list.append(deform_matrix_3)
            for deform_matrix in deform_matrix_list:
                deformed_cell = dot(cell, deform_matrix)
                deformed_cell_list.append(deformed_cell)   

        elif latt_system == 'Monoclinic':
            deform_matrix_list.append(deform_matrix_1)
            deform_matrix_list.append(deform_matrix_2)
            deform_matrix_list.append(deform_matrix_3)
            deform_matrix_list.append(deform_matrix_4)
            deform_matrix_list.append(deform_matrix_5)
            for deform_matrix in deform_matrix_list:
                deformed_cell = dot(cell, deform_matrix)
                deformed_cell_list.append(deformed_cell) 

        elif latt_system == 'Triclinic':
            deform_matrix_list.append(deform_matrix_1)
            deform_matrix_list.append(deform_matrix_2)
            deform_matrix_list.append(deform_matrix_3)
            deform_matrix_list.append(deform_matrix_4)
            deform_matrix_list.append(deform_matrix_5)
            deform_matrix_list.append(deform_matrix_6)
            for deform_matrix in deform_matrix_list:
                deformed_cell = dot(cell, deform_matrix)
                deformed_cell_list.append(deformed_cell)
                
        else:
            print('Crystal system is not parsed correctly!!!')
            exit(1)

    elif indict['dimensional'][0] == '2D':
        deform_matrix_1 = identity_matrix + strain_matrix(latt_system, up)[0]
        deform_matrix_2 = identity_matrix + strain_matrix(latt_system, up)[1]
        deform_matrix_3 = identity_matrix + strain_matrix(latt_system, up)[2]

        if latt_system == 'isotropy' or latt_system == 'tetragonal':
            deformed_cell = dot(cell, deform_matrix_1)
            deformed_cell_list.append(deformed_cell)

        elif latt_system == 'orthotropy':
            deform_matrix_list.append(deform_matrix_1)
            deform_matrix_list.append(deform_matrix_2)
            for deform_matrix in deform_matrix_list:
                deformed_cell = dot(cell, deform_matrix)
                deformed_cell_list.append(deformed_cell)

        elif latt_system == 'anisotropy':
            deform_matrix_list.append(deform_matrix_1)
            deform_matrix_list.append(deform_matrix_2)
            deform_matrix_list.append(deform_matrix_3)
            for deform_matrix in deform_matrix_list:
                deformed_cell = dot(cell, deform_matrix)
                deformed_cell_list.append(deformed_cell)

        else:
            print('Crystal system is not parsed correctly!!!')
            exit(1)

    elif indict['dimensional'][0] == '1D':
        if latt_system == 'Nanotube':
            deform_matrix = identity_matrix + strain_matrix(latt_system, up)[0]
            deformed_cell = dot(cell, deform_matrix)
            deformed_cell_list.append(deformed_cell)
            
    return deformed_cell_list
