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
from ase.io import vasp
from read_input import indict
from vasp_run import vasp_run
from extract_mean_values import mean_stress


def calc_stress(pos_optimized, cell_strain, method_stress_statistics, stress_set_dict, num_last_samples, up, cwd):
    pos_strain = deepcopy(pos_optimized)
    pos_strain.set_cell(cell_strain, scale_atoms=True)

    if method_stress_statistics == 'dynamic':
        #pos_supercell = pos_strain.repeat((int(indict['repeat_num'][0]),int(indict['repeat_num'][1]),int(indict['repeat_num'][2])))        
        step = 'NVT-MD'
        tag = 'Total+kin.'
        kpoints_file_name = 'KPOINTS-dynamic'
        #vasp.write_vasp('POSCAR', pos_supercell, vasp5=True, direct=True)
        vasp.write_vasp('POSCAR', pos_strain, vasp5=True, sort=True, direct=True)
    else:
        vasp.write_vasp('POSCAR', pos_strain, vasp5=True, sort=True, direct=True)
        step ='fixed-volume-opt'
        tag = 'in kB'
        kpoints_file_name = 'KPOINTS-static'
        
    # calculate stresses via static or molecular dynamics at fixed pressure and/or temperature
    vasp_run(step, kpoints_file_name, cwd)

    run_mode = int(indict['run_mode'][0])
    if run_mode == 1 or run_mode == 3:
        stress = mean_stress('OUTCAR', num_last_samples, tag)
        stress_set_dict[up].append(stress)

    return stress_set_dict
