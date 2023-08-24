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

from os import mkdir, chdir
from os.path import isdir
from copy import deepcopy
from ase.io import vasp
from read_input import indict
from vasp_run import vasp_run


def relax_atoms_pos(pos_optimized, cell_strain, cwd):
    if not isdir('RELAX'):
        mkdir('RELAX')
    chdir('RELAX')

    #print()
    #print(cell_strain)
    pos_strain = deepcopy(pos_optimized)
    pos_strain.set_cell(cell_strain, scale_atoms=True)
    vasp.write_vasp('POSCAR', pos_strain, vasp5=True, direct=True)

    # relax atoms positions under fixed strains
    kpoints_file_name = 'KPOINTS-static'
    pos_conv_strain = vasp_run('relax', kpoints_file_name, cwd)
    chdir('..')

    return pos_conv_strain
    
