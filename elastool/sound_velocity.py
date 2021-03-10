"""
  Elastool -- Elastic toolkit for finite-temperature elastic constants calculations

  Copyright (C) 2019-2020 by Zhong-Li Liu

  This program is free software; you can redistribute it and/or modify it under the
  terms of the GNU General Public License as published by the Free Software Foundation
  version 3 of the License.

  This program is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
  PARTICULAR PURPOSE.  See the GNU General Public License for more details.

  E-mail: zl.liu@163.com
"""
from ase.io import vasp
from math import pi


def sound_velocity(elastic_constants_dict, cwd):
    try:
        pos = vasp.read_vasp('%s/OPT/CONTCAR'%cwd)
        # The Planck's constant in m^2 Kg s^-1
        h = 6.626E-34
        # The Boltzmann constant in m^2 Kg s^-2 K-1
        k = 1.381E-23
        # The Avogadro's Constant
        Na = 6.02E+23
        # The number of atoms in the supercell (molecule)
        n = pos.get_global_number_of_atoms()
        # The total volume of the supercell (molecule)
        volume = pos.get_volume()
        # The total mass of all atoms in the supercell (molecule), in AMU
        M = sum(pos.get_masses())
        # The density in Kg/m^3
        rho = M*1E-3/Na/volume/1E-30

        B = elastic_constants_dict['B_vrh']
        G = elastic_constants_dict['G_vrh']
        V_s = 1E-3*G**0.5*(1E+9/(M*1E-3/volume/1E-30/(6.02E+23)))**0.5
        V_b = 1E-3*B**0.5*(1E+9/(M*1E-3/volume/1E-30/(6.02E+23)))**0.5
        V_p = 1E-3*(B+4.*G/3.)**0.5*(1E+9/(M*1E-3/volume/1E-30/(6.02E+23)))**0.5
        V_m = ((2./V_s**3.+1./V_p**3)/3.)**(-1./3.)
        T_D = h/k*(3.*n/4./pi*Na*rho/M/1E-3)**(1./3.)*V_m*1E+3

        elastic_constants_dict['V_s'] = V_s
        elastic_constants_dict['V_b'] = V_b
        elastic_constants_dict['V_p'] = V_p
        elastic_constants_dict['V_m'] = V_m
        elastic_constants_dict['T_D'] = T_D
    except:
        pass

    return elastic_constants_dict


if __name__ == '__main__':
    elastic_constants_dict = {'B_vrh': 34.99, 'G_vrh': 17.20}
    cwd = '.'
    elastic_constants_dict= sound_velocity(elastic_constants_dict, cwd)
    print(elastic_constants_dict)
