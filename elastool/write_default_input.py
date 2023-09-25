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

def write_default_input(method_stress_statistics, cwd):
    
    if method_stress_statistics == "static":
        # Static KPOINTS content
        kpoints_content = """KPOINT-values
0
Gamma
   5   5   1
0.0  0.0  0.0
"""
        # Write KPOINTS content to a file
        with open(os.path.join(cwd, "KPOINTS-static"), "w") as kpoints_file:
            kpoints_file.write(kpoints_content)

        # INCAR content for static
        incar_content = """########-------- Step: fixed-pressure-opt --------########
PREC    = Accurate
ENCUT   = 550
EDIFF   = 1e-6
EDIFFG  = -0.02

IBRION  = 2
ISIF    = 4
ISYM    = 2
NSW     = 250
ISMEAR  = 0
SIGMA   = 0.1
POTIM   = 0.1
PSTRESS = 0.001

NPAR    = 8
NSIM    = 4
ALGO    = Normal
IALGO   = 48
ISTART  = 0
#LPLANE  = .TRUE.
LCHARG  = .FALSE.
LWAVE   = .FALSE.
KLANE = .FALSE.
KPAR = 8
########-------- Step: fixed-volume-opt --------########
PREC    = Accurate
ENCUT   = 550
EDIFF   = 1e-6

IBRION  = 2
ISIF    = 2
ISYM    = 2
NSW     = 250
ISMEAR  = 0
SIGMA   = 0.1
POTIM   = 0.1

NPAR    = 8
NSIM    = 4
ALGO    = Normal
IALGO   = 48
ISTART  = 0

#LPLANE  = .TRUE.
LCHARG  = .FALSE.
LWAVE   = .FALSE.
KLANE = .FALSE.
KPAR = 8
"""
    
    elif method_stress_statistics == "dynamic":
        # Dynamic KPOINTS content
        kpoints_content_d = """KPOINT-values
0
Gamma
   5   5   1
0.0  0.0  0.0
"""
        # Write KPOINTS content to a file
        with open(os.path.join(cwd, "KPOINTS-dynamic"), "w") as kpoints_file_d:
            kpoints_file_d.write(kpoints_content_d)
# 

        kpoints_content_s = """KPOINT-values
0
Gamma
   5   5   1
0.0  0.0  0.0
"""
        # Write KPOINTS content to a file
        with open(os.path.join(cwd, "KPOINTS-static"), "w") as kpoints_file_s:
            kpoints_file_s.write(kpoints_content_s)

        # INCAR content for dynamic
        incar_content = """
#######-------- Step: fixed-pressure-opt --------########
PREC    = Accurate
ENCUT   = 550
EDIFF   = 1e-6
EDIFFG  = -0.001

IBRION  = 2
ISIF    = 2
ISYM    = 2
NSW     = 50
ISMEAR  = 2
SIGMA   = 0.2
POTIM   = 0.1
PSTRESS = 500
ISPIN     = 1
NPAR    = 11
NSIM    = 4
ALGO    = Normal
IALGO   = 48
ISTART  = 0

LPLANE  = .TRUE.
LCHARG  = .FALSE.
LWAVE   = .FALSE.
IWAVPR  = 11
########-------- Step: fixed-volume-opt --------########
PREC    = Accurate
ENCUT   = 550
EDIFF   = 1e-6
EDIFFG  = -0.001
ISPIN     = 1
IBRION  = 2
ISIF    = 2
ISYM    = 2
NSW     = 50
ISMEAR  = 2
SIGMA   = 0.2
POTIM   = 0.1

NPAR    = 11
NSIM    = 4
ALGO    = Normal
IALGO   = 48
ISTART  = 0

LPLANE  = .TRUE.
LCHARG  = .FALSE.
LWAVE   = .FALSE.
IWAVPR  = 11
PSTRESS = 500

########-------- Step: NPT-MD --------########
ENCUT   = 550
EDIFF   = 1E-4
ALGO    = Normal
IALGO   = 48
MAXMIX  = 40
IBRION  = 0
NSW     = 100
NBLOCK  = 1
KBLOCK  = 10
POTIM   = 2
ISYM    = 0

# NPT ensemble
ISIF    = 4
MDALGO  = 3
PSTRESS = 500
TEBEG   = 300
PMASS   = 5000
LANGEVIN_GAMMA = 10 10
LANGEVIN_GAMMA_L = 1 1

LREAL   = False
NELMIN  = 4
PREC    = Normal
ISTART  = 0
ISMEAR  = 2
SIGMA   = 0.2

NPAR    = 11
NCORE   = 1
NSIM    = 4
NWRITE  = 0

LCHARG  = .FALSE.
LPLANE  = .TRUE.
LWAVE   = .FALSE.
IWAVPR  = 11
ISPIN     = 1
########-------- Step: NVT-MD --------########
ENCUT   = 550
EDIFF   = 1E-4
ALGO    = Normal
IALGO   = 48
MAXMIX  = 40
IBRION  = 0
NSW     = 100
NBLOCK  = 1
KBLOCK  = 10
POTIM   = 2
ISYM    = 0

# NVT ensemble
ISIF    = 2
SMASS   = 2
MDALGO  = 2
TEBEG   = 300
PSTRESS = 500
LREAL   = False
NELMIN  = 4
PREC    = Normal
ISTART  = 0
ISMEAR  = 2
SIGMA   = 0.2

NPAR    = 11
NCORE   = 1
NSIM    = 4
NWRITE  = 0

LCHARG  = .FALSE.
LPLANE  = .TRUE.
LWAVE   = .FALSE.
IWAVPR  = 11
ISPIN     = 1
"""

    else:
        print("The given method_stress_statistics is not supported yet.")
        return
    
    # Write INCAR content to a file
    with open(os.path.join(cwd, "INCARs"), "w") as incar_file:
        incar_file.write(incar_content)

    #print_default_input_message()


def write_default_elastool_in(cwd):
    # Check if "elastool.in" already exists in cwd
    if not os.path.exists(os.path.join(cwd, "elastool.in")):
        
        elastool_in = """###############################################################################
### The input file to control the calculation details of elastic constants  ###
###############################################################################

# run mode: 0 for generating key default input files, 1 for automatic run, 2 for pre-processing, 3 for post-processing
# if 2, plz ensure the structure opt. is performed at fixed pressure or volume
# i.e. CONTCAR and OUTCAR files both exist in ./OPT directory.
run_mode = 0

# Define the dimensional of the system: 1D/2D/3D. 
dimensional = 2D

# the crystal structure file in vasp POSCAR (.vasp) or cif (.cif) format
structure_file = filename.vasp

# if use conventional cell, no for primitive cell, yes for conventional cell
if_conventional_cell = no

# static or dynamic, static for 0 K, dynamic for finite-temperature
method_stress_statistics = static

# strains matrix for solve all elastic constants, asess or ohess or ulics
strains_matrix = ohess

# strains list for deforming lattice cell, 0 will be neglected because of 
# the initial optimization, if method_statistics = dynamic, the first one is used
strains_list = -0.06 -0.03 0.03 0.06  

# repeat numbers of three lattice vectors in conventional lattice for making
# supercell of molecular dynamics simulations (method_statistics = dynamic)
repeat_num = kx ky kz

# last number of steps for sampling stresses used in the dynamic method
num_last_samples = 500

#Potential directory - specify the location of your POTCAR files
potential_dir = /user/potential/

# The parallel submiting commmd
parallel_submit_command = vasp_cmd > log.vasp
"""
        # Write to "elastool.in" file
        with open(os.path.join(cwd, "elastool.in"), "w") as elin:
            elin.write(elastool_in)




def display_help():
    help_text = """---------------------------
ElasTool Help Documentation
---------------------------

Copyright (C) 2019-2024 by Zhong-Li Liu and Chinedu Ekuma

**Disclaimer:**
This program is free software, and you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, version 3 of the License.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

---------------------------------------
**Configuration and Input Parameters:**
---------------------------------------

1. run_mode: Controls the operation mode of ElasTool.
   - `0`: Generates key default input files.
   - `1`: Automatic run.
   - `2`: Pre-processing (Ensure structure optimization is performed at fixed pressure or volume, i.e., CONTCAR and OUTCAR files both exist in ./OPT directory).
   - `3`: Post-processing.

2. dimensional: Define the dimensions of the system. Options are `1D`, `2D`, or `3D`.

3. structure_file: Specifies the crystal structure file in VASP POSCAR (.vasp) or CIF (.cif) format.

4. if_conventional_cell: Set to `yes` for conventional cell, `no` for primitive cell.

5. method_stress_statistics: Define the method of stress statistics. Options are `static` for 0 K, `dynamic` for finite-temperature.

6. strains_matrix: Specify the strains matrix to solve all elastic constants. Options are `asess`, `ohess`, or `ulics`.

7. strains_list: Define the strains list for deforming lattice cell. `0` will be neglected because of the initial optimization. If `method_statistics` is `dynamic`, the first one is used.

8. repeat_num: Define the repeat numbers of three lattice vectors in conventional lattice for making supercell of molecular dynamics simulations (if `method_statistics` is `dynamic`).

9. num_last_samples: Specify the last number of steps for sampling stresses used in the dynamic method.

10. potential_dir: Define the directory where your POTCAR files are located.

11. parallel_submit_command: Specify the parallel submitting command.

-------------------------
**Example Configuration:**
-------------------------

run_mode = 0
dimensional = 2D
structure_file = filename.vasp
if_conventional_cell = no
method_stress_statistics = dynamic
strains_matrix = ohess
strains_list = -0.06 -0.03 0.03 0.06
repeat_num = kx ky kz
num_last_samples = 500
potential_dir = /user/potential/
parallel_submit_command = vasp_cmd > log.vasp
```

Please ensure to replace placeholder values with your specific values, especially for parameters like `structure_file`, `repeat_num`, `potential_dir`, and `parallel_submit_command`. Always consult the relevant documentation and guidelines for your system and simulation needs.
    """
    print(help_text)





def print_default_input_message_0():
    print("╔════════════════════════════════════════════════════════════════════════════════╗")
    print("║                                                                                ║")
    print("║                       ♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥                       ║")
    print("║                  ♥♥♥  Default elastool.in template generated.  ♥♥♥             ║")
    print("║                 ♥♥ Modify and rerun elastool -0 to generate other ♥♥           ║")
    print("║                 ♥♥    important input files. Happy running :)    ♥♥            ║")
    print("║                       ♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥                       ║")
    print("║                                   Exiting...                                   ║")
    print("║                                                                                ║")
    print("╚════════════════════════════════════════════════════════════════════════════════╝")





def print_default_input_message_1():
    print("╔════════════════════════════════════════════════════════════════════════════════╗")
    print("║                                                                                ║")
    print("║                               ♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥                                  ║")
    print("║                           ♥♥                    ♥♥                             ║")
    print("║                        ♥ All default inputs written to files. ♥                ║")
    print("║                        ♥ Modify according to dimensionality ♥                  ║")
    print("║                        ♥ and other criteria ♥                                  ║")
    print("║                        ♥       Happy running :)        ♥                       ║")
    print("║                           ♥♥                    ♥♥                             ║")
    print("║                               ♥♥♥♥♥♥♥♥♥♥♥♥♥♥♥                                  ║")
    print("║                                Exiting...                                      ║")
    print("║                                                                                ║")
    print("╚════════════════════════════════════════════════════════════════════════════════╝")





