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
import os
from ase.io import vasp, cif,read, write
from os import system
from os.path import isfile
import subprocess
from time import sleep
from write_incar import write_incar
from read_input import indict
from extract_mean_values import get_pressure, mean_pressure




def vasp_run(step, kpoints_file_name, cwd):
    write_incar(step, cwd)
    os.system("cp %s/%s KPOINTS" % (cwd, kpoints_file_name))

    #structure_file = indict['structure_file'][0]
    structure_file = os.path.join(cwd, indict['structure_file'][0])
    
    if structure_file.endswith('.vasp'):
        pos = vasp.read_vasp(structure_file)
    elif structure_file.endswith('.cif'):
        #pos = cif.read_cif(structure_file,index=0)

        #atoms = read(structure_file)
        # Convert to POSCAR format
        #poscar_file = os.path.join(cwd, 'POSCAR')
        #write(poscar_file, atoms, format='vasp')
        pos = read(structure_file) #vasp.read_vasp(poscar_file)
    else:
        raise ValueError("Unsupported structure file format!")

    chem_symb = pos.get_chemical_symbols()
    
    # Determine unique elements
    ele_list = []
    for i in chem_symb:
        if i not in ele_list:
            ele_list.append(i)

    if os.path.isfile('POTCAR'):
        os.system("rm POTCAR")

    potential_directory = indict.get('potential_dir', ['./'])[0]
    

    for ele in ele_list:
        # Path possibilities
        potential_potcar_paths = [
            os.path.join(potential_directory, ele, "POTCAR"),
            os.path.join(potential_directory, ele + "_pv", "POTCAR"),
            os.path.join(potential_directory, ele + "_sv", "POTCAR"),
            os.path.join(potential_directory, ele + "_GW", "POTCAR")
        ]

        pot_file_path = next((path for path in potential_potcar_paths if os.path.exists(path)), None)

        if pot_file_path:
            os.system(f"cat {pot_file_path} >> POTCAR")
        else:
            raise Exception(f"POTCAR for {ele} not found in any of the expected directories!")

    for line in open('INCAR','r'):
        if 'PSTRESS' in line:
            # p_target = 0.1 * float(line.split()[2])
            break

    if int(indict['run_mode'][0]) == 1:
        # computcode = indict['parallel_submit_command'][3].split("/")[-1]
        para_sub_com = ''
        for i in range(len(indict['parallel_submit_command'])):
            para_sub_com += indict['parallel_submit_command'][i]
            para_sub_com += ' '

        #para_sub_com = 'yhrun -p paratera -N 1 -n 24 vasp_std'
        #print(para_sub_com)
        #start = datetime.datetime.now()

        go = subprocess.Popen(para_sub_com, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        while go.poll() is None:
            sleep(2)
            #if method_stress_statistics == 'dynamic' and 'NO_STRAIN_MD' in os.getcwd():
            #    p = get_pressure('OUTCAR')
            #    mean_press = 0.1 * mean_pressure('OUTCAR', 1)

            #    if mean_press is not np.nan and len(p) > 2:
            #        if abs(mean_press - p_target) < 0.2:
                        #os.system("killall -9 $(ps H -e -o cmd --sort=pcpu | tail -1)")
            #            time.sleep(10)
            #            os.system("killall -9 vasp544")

            #now = datetime.datetime.now()
            #if (now - start).seconds > float(indict['MaxHour'][0])*3600:
            #    os.system("killall -9 $(ps H -e -o cmd --sort=pcpu | tail -1)")
            #    break
        pos_optimized = vasp.read_vasp('CONTCAR')

    elif int(indict['run_mode'][0]) == 2:
        pos_optimized = pos

    elif int(indict['run_mode'][0]) == 3:
        pos_optimized = vasp.read_vasp('CONTCAR')

    return pos_optimized
