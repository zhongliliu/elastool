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
from math import pi, sqrt
from numpy import cross, linalg


def sound_velocity(elastic_constants_dict, cwd, dimensional, latt_system):
    pos = vasp.read_vasp('%s/OPT/CONTCAR' % cwd)
    # The Planck's constant in m^2 Kg s^-1
    h = 6.626E-34
    # The reduced Planck's constant hbar in m^2 Kg s^-1
    hbar = 1.0546E-34
    # The Boltzmann constant in m^2 Kg s^-2 K-1
    k = 1.381E-23
    # The Avogadro's constant
    Na = 6.02E+23
    # The total volume of the system cell 
    volume = pos.get_volume()
    # The total mass of all atoms in the material, in AMU
    M = sum(pos.get_masses())
    # The number of atoms in the system (material)
    n = pos.get_global_number_of_atoms()
    if dimensional == '3D':
        try:
            # The density in Kg/m^3
            rho = M * 1E-3 / Na / volume / 1E-30

            B = elastic_constants_dict['B_vrh']
            G = elastic_constants_dict['G_vrh']
            V_s = 1E-3 * G**0.5 * \
                (1E+9 / (M * 1E-3 / volume / 1E-30 / (6.02E+23)))**0.5
            V_b = 1E-3 * B**0.5 * \
                (1E+9 / (M * 1E-3 / volume / 1E-30 / (6.02E+23)))**0.5
            V_p = 1E-3 * (B + 4. * G / 3.)**0.5 * (1E+9 / \
                          (M * 1E-3 / volume / 1E-30 / (6.02E+23)))**0.5
            V_m = ((2. / V_s**3. + 1. / V_p**3) / 3.)**(-1. / 3.)
            T_D = (h / k) * (3. * n / 4. / pi * Na * rho / M / 1E-3)**(1. / 3.) * V_m * 1E+3

            elastic_constants_dict['V_s'] = V_s
            elastic_constants_dict['V_b'] = V_b
            elastic_constants_dict['V_p'] = V_p
            elastic_constants_dict['V_m'] = V_m
            elastic_constants_dict['T_D'] = T_D
        except BaseException:
            pass
    elif dimensional == '2D':
        c11 = elastic_constants_dict['c11']
        c12 = elastic_constants_dict['c12']
#        cii = 0
        if 'c22' in elastic_constants_dict:
            c22 = elastic_constants_dict['c22']
#            cii = (c11+c22)/2.0
        else:
            c22 = elastic_constants_dict['c11']  # Lattice a=b

        # Thomas et al., RSC Adv., 2018, 8, 27283
        Y_2Da = (c11 * c22 - c12**2) / c11
        Y_2Db = (c11 * c22 - c12**2) / c22
        v_a = c12 / c11
        v_b = c12 / c22

        # Note B is the in-plane stiffness (2D analogy of the bulk modulus)
        B_a = Y_2Da / 2.0 / (1 - v_a)
        G_a = Y_2Da / 2.0 / (1 + v_a)  # Note you can use (C_ii-C_ij)/2
        B_b = Y_2Db / 2.0 / (1 - v_b)
        # Note you can use (C_ii-C_ij)/2        cell = pos.get_cell()
        G_b = Y_2Db / 2.0 / (1 + v_b)

        # PRB 94, 245420 (2016)
        # Calculate sound velocities in km/s

        cell = pos.get_cell()
        # The 2D density in Kg/m^2
        area = linalg.norm(cross(cell[0], cell[1]))
        rho_2D = M * 1E-3 / Na / area / 1E-20


        # 1E-3*sqrt(Y_2D*(1-v)/rho_2D/(1+v)/(1-2*v))
        V_la = 1E-3 * sqrt(abs(B_a + G_a) / rho_2D)
        V_sa = 1E-3 * sqrt(abs(G_a) / rho_2D)

        # 1E-3*sqrt(Y_2D*(1-v_b)/rho_2D/(1+v_b)/(1-2*v_b)) 
        V_lb = 1E-3 * sqrt(abs(B_b + G_b) / rho_2D)
        V_sb = 1E-3 * sqrt(abs(G_b) / rho_2D)

        V_ma = (1.0 / 2.0 * (1.0 / V_sa**2.0 + 1 / V_la**2.0))**(-1.0 / 2.0) #(1.0 / 3.0 * (2.0 / V_sa**3.0 + 1 / V_la**3.0))**(-1.0 / 3.0)
        T_Da = (hbar / k) * (4. * pi * n / area / 1E-20)**(1. / 2.) * V_ma * 1E+3

        V_mb = (1.0 / 2.0 * (1.0 / V_sb**2.0 + 1 / V_lb**2.0))**(-1.0 / 2.0)
        T_Db = (hbar / k) * (4. * pi * n / area / 1E-20)**(1. / 2.) * V_mb * 1E+3

        elastic_constants_dict['V_la'] = V_la
        elastic_constants_dict['V_lb'] = V_lb
        elastic_constants_dict['V_sa'] = V_sa
        elastic_constants_dict['V_sb'] = V_sb
        elastic_constants_dict['V_ma'] = V_ma
        elastic_constants_dict['V_mb'] = V_mb
        elastic_constants_dict['T_Da'] = T_Da
        elastic_constants_dict['T_Db'] = T_Db
        elastic_constants_dict['Y_2Da'] = Y_2Da
        elastic_constants_dict['Y_2Db'] = Y_2Db
        elastic_constants_dict['G_a'] = G_a
        elastic_constants_dict['G_b'] = G_b
        elastic_constants_dict['v_a'] = v_a
        elastic_constants_dict['v_b'] = v_b
        elastic_constants_dict['B_a'] = B_a
        elastic_constants_dict['B_b'] = B_b

    elif dimensional == '1D':
        cell = pos.get_cell()
        # We have1D structure is oriented along cell[2]
        length = linalg.norm(cell[2])
        rho_1D = M * 1E-3 / Na / length  / 1E-10
        mass_density = M * 1E-3 / Na / volume / 1E-30 # kg/m^3
        
       # M = pos.get_masses().mean()  # Mean atomic mass in the unit cell
       # n = len(pos)  # Number of atoms in the unit cell
       # Na = 6.02214076E+23  # Avogadro's number
 
        if latt_system == 'Nanotube': 

            c33 = elastic_constants_dict['c33']
            c23 = elastic_constants_dict['c23']
            Y = c33-c23
            v = -c23 / c33
            K = (c33-2*c23)/3.
            G = (c33-c23)/2.
            V = 1E-3* sqrt(abs(c33*1E+9) / mass_density)
            T_D = (h / k) * (3. * n / 4. / pi * Na * mass_density / M / 1E-3)**(1. / 3.) * V * 1E+3 
            R_f = sqrt(abs(c33) * 1E+9 / mass_density) / (2 * pi * length* 1E-10 )/1E+9

            elastic_constants_dict['Y'] = Y
            elastic_constants_dict['B'] = K # Bulk modulus
            elastic_constants_dict['G'] = G #Shear modulus
            #elastic_constants_dict['v'] = v
            elastic_constants_dict['V'] = V
            elastic_constants_dict['T_D'] = T_D
            #elastic_constants_dict['R_f'] = R_f



    return elastic_constants_dict


if __name__ == '__main__':
    elastic_constants_dict = {'B_vrh': 34.99, 'G_vrh': 17.20}
    cwd = '.'
    elastic_constants_dict = sound_velocity(elastic_constants_dict, cwd)
    print(elastic_constants_dict)
