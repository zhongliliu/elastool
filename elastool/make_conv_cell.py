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
import numpy as np
from ase.geometry import cell_to_cellpar, cellpar_to_cell
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

def remove_spurious_distortion(pos):
    # Normalize and orthogonalize the cell vectors
    cell_params = cell_to_cellpar(pos.get_cell())
    new_cell = cellpar_to_cell(cell_params)
    pos.set_cell(new_cell, scale_atoms=True)

    # Adjust atom positions
    pos.wrap()

    pos.center()

    return pos
    

import numpy as np
from ase.io import read, write


def swap_axes_to_longest_c(pos):
    """
    Reorders the lattice vectors of the given structure (pos) such that the longest vector is always c-axis.

    Parameters:
    pos (ASE Atoms object): The atomic structure.

    Returns:
    ASE Atoms object: Updated structure with reordered lattice vectors.
    """
    a, b, c = pos.get_cell()
    lengths = [np.linalg.norm(a), np.linalg.norm(b), np.linalg.norm(c)]
    max_index = lengths.index(max(lengths))

    new_cell = [a, b, c]  
    new_cell[2], new_cell[max_index] = new_cell[max_index], new_cell[2]  

    if max_index == 0:
        new_cell[0], new_cell[1] = new_cell[1], new_cell[0]
    pos.set_cell(new_cell)

    return pos


    
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
    nanoribbon = False
    cell = pos.get_cell()
    a_length = np.linalg.norm(cell[0])
    b_length = np.linalg.norm(cell[1])
    c_length = np.linalg.norm(cell[2])
    
    is_conv = indict.get('if_conventional_cell', ['no'])[0]
    if indict['dimensional'][0] == '3D':       
        if is_conv == 'yes':
            to_primitive = False
        elif is_conv == 'no':
            to_primitive = True
            

        if max(a_length, b_length, c_length) > 12.5: #We define a cut to check for a potential nanoribbon structure
            #print("Lattice parameter suggests a nanoribbon structure")
            nanoribbon = True
            
        if nanoribbon:
            lengths = [a_length, b_length, c_length]
            max_index = np.argmax(lengths)
            # Swap the cell vectors to ensure the longest is along the c-axis
            pos = swap_axes_to_longest_c(pos)
            pos_conv = remove_spurious_distortion(pos)

        else:
            pos_std = standardize_cell(pos, to_primitive=to_primitive,symprec=1e-5)
            pos_conv = Atoms(pos_std[2], scaled_positions=pos_std[1], cell=pos_std[0],pbc=True)
            
            if not to_primitive:
                pmg_structure = AseAtomsAdaptor.get_structure(pos_conv)
                sga = SpacegroupAnalyzer(pmg_structure)
                conventional_pmg_structure = sga.get_conventional_standard_structure()
                pos_conv = AseAtomsAdaptor.get_atoms(conventional_pmg_structure)
            pos_conv = remove_spurious_distortion(pos_conv)

    else:
        pos_conv = remove_spurious_distortion(pos)
        #pos_conv = pos


    return pos_conv, nanoribbon

    
