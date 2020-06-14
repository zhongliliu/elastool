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


def sound_velocity(elastic_constants_dict, cwd):
    try:
        pos = vasp.read_vasp('%s/OPT/CONTCAR'%cwd)
        volume = pos.get_volume()/0.529177**3.0
        M = sum(pos.get_masses())
        B = elastic_constants_dict['B_vrh']
        G = elastic_constants_dict['G_vrh']

        V_s = 0.001*G**0.5*(10**9/(M*10**(-3)/volume/(0.529177*10**(-10))**3/(6.02*10**23)))**0.5
        V_b = 0.001*B**0.5*(10**9/(M*10**(-3)/volume/(0.529177*10**(-10))**3/(6.02*10**23)))**0.5
        V_p = 0.001*(B+4.*G/3.)**0.5*(10**9/(M*10**(-3)/volume/(0.529177*10**(-10))**3/(6.02*10**23)))**0.5
        V_m = ((2./V_s**3.+1/V_p**3)/3.)**(-1./3.)

        T_D = (6.626*10**(-34)/(1.381*10**(-23)))*(3./(volume*(.529177*10**(-10))**3.)/4./3.1415926)**(1./3.)*V_m*1000

        elastic_constants_dict['V_s'] = V_s
        elastic_constants_dict['V_b'] = V_b
        elastic_constants_dict['V_p'] = V_p
        elastic_constants_dict['V_m'] = V_m
        elastic_constants_dict['T_D'] = T_D
    except:
        pass

    return elastic_constants_dict
    