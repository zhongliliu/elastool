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

from ase.io import vasp
from spglib import standardize_cell
from ase import io, Atoms
from read_input import indict


def make_conventional_cell(file_name_read):
    if file_name_read.split('.')[-1] == 'vasp':
        input_file_format = 'vasp'
    elif file_name_read.split('.')[-1] == 'cif':
        input_file_format = 'cif'
    else:
        print('Please provide correct file name: .vasp or .cif')
        exit(1)

    file_read = open(file_name_read, 'r')

    if input_file_format == 'vasp':
        pos = vasp.read_vasp(file_read)
    elif input_file_format == 'cif':
        pos = io.read(file_read, format='cif')

    if indict['dimensional'][0] == '3D':
        if indict['if_conventional_cell'][0] == 'yes':
            to_primitive = False
        elif indict['if_conventional_cell'][0] == 'no':
            to_primitive = True

        pos_std = standardize_cell(pos, to_primitive=to_primitive)
        pos_conv = Atoms(pos_std[2], scaled_positions=pos_std[1], cell=pos_std[0])

    else:
        pos_conv = pos

    return pos_conv
    
