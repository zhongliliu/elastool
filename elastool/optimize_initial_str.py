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
from ase.io import vasp,write
from vasp_run import vasp_run
from read_input import indict
import os
import shutil


#def optimize_initial_str(pos_conv, cwd, tag):
#    if not isdir('OPT'):
#        mkdir('OPT')
#    chdir('OPT')

#    #vasp.write_vasp('POSCAR', pos_conv, vasp5=True, direct=True)
#    write('POSCAR', pos_conv, format='vasp', direct=True)
#    kpoints_file_name = 'KPOINTS-static'
#    pos_optimized = vasp_run(tag, kpoints_file_name, cwd)

#    chdir('..')

#    return pos_optimized
    
    
    

def optimize_initial_str(pos_conv, cwd, tag, fresh=False, max_retries=3):
#    if not os.path.isdir('OPT'):
#        os.mkdir('OPT')
#    os.chdir('OPT')

    opt_dir = os.path.join(cwd, 'OPT')
    
    if fresh:
        if os.path.isdir(opt_dir):
            shutil.rmtree(opt_dir)
        os.mkdir(opt_dir)
    else:
        if not os.path.isdir(opt_dir):
            os.mkdir(opt_dir)

    os.chdir(opt_dir)
    
    attempt = 0
    while attempt < max_retries:
        try:
            write('POSCAR', pos_conv, format='vasp', direct=True)
            kpoints_file_name = 'KPOINTS-static'
            pos_optimized = vasp_run(tag, kpoints_file_name, cwd)
            break  

        except Exception as e:
            attempt += 1
            print(f"An error occurred on attempt {attempt}: {e}")
            if attempt < max_retries and os.path.exists('CONTCAR') and os.path.getsize('CONTCAR') > 0:
                print("Attempting to restart from CONTCAR...")
                shutil.copy('CONTCAR', 'POSCAR')
            #else:
            #    print("Maximum retries reached or CONTCAR file not found/empty.")
            #    raise
        finally:
            os.chdir('..')

    return pos_optimized #if attempt < max_retries else None

