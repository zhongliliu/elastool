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
from spglib import spglib
from read_input import indict

def findspg(atoms):
    spg0 = spglib.get_spacegroup(atoms, symprec=0.1)
    if spg0:
        spg1 = spg0.split()
        spg = [str(spg1[0]), int(spg1[1][1:-1])]
    else:
        spg = []
    # print spg0, spg

    return spg


def find_crystal_system(pos_conv, dimensional,tubestrain_type):
    
    if dimensional == '1D':
        latt_system = 'Nanotube'  # Default value
        dist_acc = 0.1
        lenth_angl = pos_conv.get_cell_lengths_and_angles()
        a = lenth_angl[0]
        b = lenth_angl[1]
        c = lenth_angl[2]


        if tubestrain_type == "Nanotube":
            latt_system = 'Nanotube'  # This is 1D embedded in 3D space

        else:
            print('ERROR: Define appropriate nanotube structure!!!\n')
 

    elif dimensional == '2D':
        dist_acc = 0.1
        angl_acc = 0.5

        lenth_angl = pos_conv.get_cell_lengths_and_angles()
        a = lenth_angl[0]
        b = lenth_angl[1]
        c = lenth_angl[2]
        gamma = lenth_angl[5]

        # The 2D lattice system is defined according to 2D Mater. 6 (2019) 048001
        if c > a and c > b:
            if abs(a - b) <= dist_acc:
                if abs(
                        gamma -
                        120) <= angl_acc or abs(
                        gamma -
                        60) <= angl_acc:  # The last part is for some 2D systems
                    latt_system = 'isotropy'  # This is 2D Hexagonal system
                elif abs(gamma - 90) <= angl_acc:
                    latt_system = 'tetragonal'
            else:
                if abs(gamma - 90) <= angl_acc:
                    latt_system = 'orthotropy'
                else:
                    latt_system = 'anisotropy'
        else:
            print('ERROR: the vacuum is not along the c axis!!!\n')
            print('Plz adjust the vacuum to be along the c axis!!!\n')
            exit(1)

    elif dimensional == '3D':
        spg = findspg(pos_conv)
        spg_num = spg[1]
        crystal_system = [
            [2, "Triclinic"],
            [15, "Monoclinic"],
            [74, "Orthorombic"],
            [88, "Tetragonal2"],
            [142, "Tetragonal1"],
            [148, "Trigonal2"],
            [167, "Trigonal1"],
            [194, "Hexagonal"],
            [230, "Cubic"]
        ]

        for l in crystal_system:
            if spg_num <= l[0]:
                latt_system = l[1]
                break

    return latt_system
