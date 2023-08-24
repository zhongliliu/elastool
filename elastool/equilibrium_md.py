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
from ase.io import vasp
from vasp_run import vasp_run
from read_input import indict


def equil_md(pos_optimized, cwd):
	if not isdir('NO_STRAIN_MD'):
		mkdir('NO_STRAIN_MD')
	chdir('NO_STRAIN_MD')
	repeat = [int(indict['repeat_num'][0]),int(indict['repeat_num'][1]),int(indict['repeat_num'][2])]
	#pos_supercell = make_supercell(pos_optimized, repeat, wrap=True, tol=1e-5)
	pos_supercell = pos_optimized.repeat(repeat)
	#vasp.write_vasp('POSCAR', pos_supercell, vasp5=True, direct=True, sort=True)
	vasp.write_vasp('POSCAR', pos_supercell, direct=True, sort=True, vasp5=True)
	kpoints_file_name = 'KPOINTS-dynamic'

	pos_optimized = vasp_run('NPT-MD', kpoints_file_name, cwd)
	chdir('..')

	return pos_optimized
    
