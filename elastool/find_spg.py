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
import os
import numpy as np
from pymatgen.core import Structure
from ase.build import cut
import matplotlib.pyplot as plt
import copy


class Bravais2D:
    """ 
    Class to analyze 2D bravais lattice type based on given structure.

    Args:
        struct: Pymatgen structure object.
        eps_r: Tolerance for comparing lengths.
        eps_a: Tolerance for comparing angles.
        numpoints: The number of desired points to plot. Must be a square number >= 9.
    """

    def __init__(self, struct, dist_acc=0.1, angl_acc=0.5, numpoints=16):
        self.struct = struct
        lattice = struct.get_cell_lengths_and_angles()
        self.a = lattice[0]
        self.b = lattice[1]
        self.c = lattice[2]
        self.dist_acc = dist_acc
        self.angl_acc = angl_acc
        self.centered = self._get_centered()
        self._numpoints = numpoints
        self.gamma = lattice[5]





    def _get_centered(self):
        copied_struct = copy.deepcopy(self.struct)
        primitive_struct = cut(copied_struct)
        if len(self.struct) != len(primitive_struct):
            return True
        else:
            return False

    @property
    def a_vec(self):
        return np.array([self.a, 0])

    @property
    def b_vec(self):
        return np.array([self.b * np.cos(self.gamma), self.b * np.sin(self.gamma)])

    @property
    def lattice_type(self):
        if self.c > self.a and self.c > self.b:
            if abs(self.a - self.b) <= self.dist_acc:
                if abs(self.gamma - 120) <= self.angl_acc or abs(self.gamma - 60) <= self.angl_acc:  
                    return 'Isotropy'  # This is 2D Hexagonal system
                elif abs(self.gamma - 90) <= self.angl_acc:
                    return 'Tetragonal' #Square
            else:
                if abs(self.gamma - 90) <= self.angl_acc:
                    return 'Orthotropy' #Rectangular
                else:
                    return 'Anisotropy'
        else:
            print('ERROR: the vacuum is not along the c axis!!!\n')
            print('Plz adjust the vacuum to be along the c axis!!!\n')
            exit(1)


    @property
    def unit_cell_area(self):
        return self.a * self.b * np.sin(self.gamma)

    @property
    def numpoints(self):
        val = round(self._numpoints**0.5, 0)**2
        if val == self._numpoints and val >= 9:
            return self._numpoints
        else:
            raise Exception("numpoints must be a square number >= 9.")


    def _find_points(self):
        """ Finds all the x and y coordinates of the lattice.
        :return: (list(list)) x, y
        """

        def f(start, stop, x_list, y_list):
            for j in np.arange(start, stop):
                for i in np.arange(start, stop):
                    vec = i*self.a_vec + j*self.b_vec
                    x_list.append(vec[0])
                    y_list.append(vec[1])
            return x_list, y_list

        p = int(self.numpoints**0.5)
        x, y = f(0, p, [], [])
        if self.centered:
            x, y = f(0.5, p - 1, x, y)
        return x, y


    def _unit_cell(self, x, y):
        """ Finds the x and y coordinates for the unit cell.
        :param x: (list) The x coordinates of the lattice points.
        :param y: (list) The y coordinates of the lattice points.
        :return: (list(list()) x, y
        """

        root = int(self.numpoints**0.5)
        return (x[0], x[1], x[1+root], x[root], x[0]), (y[0], y[1], y[1+root], y[root], y[0])


    def plot(self,filename='bravais_lattice.png'):
        """Creates a 2D scatter plot of the Bravais lattice."""
        x, y = self._find_points()
        fig, ax = plt.subplots()
        angle = f"{self.gamma:.2f}"
        title = f"Bravais Lattice: {self.lattice_type}\n|a| = {self.a:.3f}, |b| = {self.b:.3f}, \u03b8  =  {self.gamma:.3f}\u00b0"
        ax.set_title(title)
        ax.scatter(x, y, label="Lattice Points")
        ax.plot(*self._unit_cell(x, y), color="darkorange", label="Unit Cell")
        plt.legend(loc="best")
        plt.savefig(filename, format='png', dpi=300)
        #plt.show()
        plt.close(fig)

def findspg(atoms):
    spg0 = spglib.get_spacegroup(atoms, symprec=0.1)
    if spg0:
        spg1 = spg0.split()
        spg = [str(spg1[0]), int(spg1[1][1:-1])]
    else:
        spg = []
    # print spg0, spg

    return spg


def find_crystal_system(pos_conv, dimensional,tubestrain_type,plotparameters):
    
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


        if plotparameters:
            #print("Before:", pos_conv)
            bravais = Bravais2D(pos_conv)
            #print("After:", pos_conv)
            bravais_lattice_type = bravais.lattice_type
            print(f"Bravais Lattice Type: {bravais_lattice_type}")
            bravais.plot()

    elif dimensional == '3D':
#        cell = pos_conv.get_cell()
#        a_length = np.linalg.norm(cell[0])
#        b_length = np.linalg.norm(cell[1])
#        c_length = np.linalg.norm(cell[2])
#        if b_length > 4 * a_length or b_length > 4 * c_length or c_length > 4 * a_length:
#             print("Lattice parameter suggests a nanoribbon structure")
#             spg_num = 231
#             crystal_system = [
#                 [231, "Nanoribbon"]
#             ]       
#        else:
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
