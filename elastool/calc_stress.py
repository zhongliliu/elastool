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

import numpy as np
from copy import deepcopy
from ase.io import vasp, write
from read_input import indict
from vasp_run import vasp_run
from extract_mean_values import mean_stress
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import CubicSpline
from scipy.integrate import simps
import re


strain_data = []
stress_data = []
energy_data = []

def calc_stress(pos_optimized, cell_strain, method_stress_statistics, stress_set_dict, num_last_samples, up, cwd):
    pos_strain = deepcopy(pos_optimized)
    pos_strain.set_cell(cell_strain, scale_atoms=True)

    if method_stress_statistics == 'dynamic':
        #pos_supercell = pos_strain.repeat((int(indict['repeat_num'][0]),int(indict['repeat_num'][1]),int(indict['repeat_num'][2])))        
        step = 'NVT-MD'
        tag = 'Total+kin.'
        kpoints_file_name = 'KPOINTS-dynamic'
        #vasp.write_vasp('POSCAR', pos_supercell, vasp5=True, direct=True)
        #vasp.write_vasp('POSCAR', pos_strain, vasp5=True, sort=True, direct=True)
        write('POSCAR', pos_strain, format='vasp', direct=True)
    else:
        #vasp.write_vasp('POSCAR', pos_strain, vasp5=True, sort=True, direct=True)
        write('POSCAR', pos_strain, format='vasp', direct=True)
        step ='fixed-volume-opt'
        tag = 'in kB'
        kpoints_file_name = 'KPOINTS-static'
        
    # calculate stresses via static or molecular dynamics at fixed pressure and/or temperature
    vasp_run(step, kpoints_file_name, cwd)

    run_mode = int(indict['run_mode'][0])
    if run_mode == 1 or run_mode == 3:
        stress = mean_stress('OUTCAR', num_last_samples, tag)
        total_energy = calculate_total_energy('OUTCAR')
        external_pressure = extract_external_pressure('OUTCAR')
        external_pressure = abs(external_pressure)

        stress_set_dict[up].append(stress)
        strain_data.append(up)
        energy_data.append(total_energy)
        stress_data.append(external_pressure)              


    return stress_set_dict, [strain_data,energy_data,stress_data]
    
    
    
    
    
def calculate_total_energy(outcar_file):
    total_energy = None  # Initialize the total energy variable

    try:
        with open(outcar_file, 'r') as file:
            lines = file.readlines()

            for line in reversed(lines):
                if 'free  energy' in line.lower():

                    total_energy = float(line.split()[-2])
                    break

    except FileNotFoundError:
        print(f"ERROR: {outcar_file} not found.")

    return total_energy


def extract_external_pressure(outcar_filename):
    external_pressure = None
    try:
        with open(outcar_filename, "r") as outcar_file:
            lines = outcar_file.readlines()

        for line in reversed(lines):
            match = re.search(r'external pressure =\s+([\d\.\-]+)\s+kB', line)
            if match:
                external_pressure = float(match.group(1))
                break

    except FileNotFoundError:
        print(f"File '{outcar_filename}' not found.")
    except Exception as e:
        print(f"An error occurred: {str(e)}")

    return external_pressure




def calculate_total_energy_density(strain_data, energy_data, stress_data, dimensional, pos_optimized):

    energy_data = np.array(energy_data)
    stress_data = np.array(stress_data) *-0.1 # Conversion factor for vasp
    
    unique_strain, unique_indices = np.unique(strain_data, return_index=True)
    unique_stress = stress_data[unique_indices]
    unique_energy = energy_data[unique_indices]


    if dimensional == "2D":
        lenth_angl = pos_optimized.get_cell_lengths_and_angles()
        stress_data = stress_data * lenth_angl[2] / 10. # N/m
    else:
        stress_data = stress_data *1E9 # N/m^2
    # Calculate the mean stress for each unique strain value
    mean_stress = []
    mean_energy = []
    for unique_value in unique_strain:
        indices = np.where(strain_data == unique_value)
        mean_stress.append(np.mean(stress_data[indices]))
        mean_energy.append(np.mean(energy_data[indices]))

    mean_stress = np.array(mean_stress)
    mean_energy = np.array(mean_energy)

    
    
    sorted_indices = np.argsort(np.array(unique_strain))
    sorted_strain = np.array(unique_strain)[sorted_indices]
    sorted_energy = mean_energy[sorted_indices]
    sorted_stress = mean_stress[sorted_indices]
    
    # The area under the stress-strain curve => the energy density

    try:
        if len(unique_strain) <= 3:

            interp_func_stress = np.poly1d(np.polyfit(sorted_strain, sorted_stress, 1))
            interp_func_energy = np.poly1d(np.polyfit(sorted_strain, sorted_energy, 1))

            fit_strain = np.linspace(min(sorted_strain)+1e-6, max(sorted_strain)-1e-6, 100)

            derivative_stress = np.polyder(interp_func_stress)
            #energy_density_at_specific_strain = interp_func_stress(sorted_strain) * derivative_stress(sorted_strain)
            total_energy_density = simps(interp_func_stress(fit_strain) * derivative_stress(fit_strain), fit_strain )
            
            total_strain_energy = simps(interp_func_energy(fit_strain), fit_strain)
            #total_energy_density = np.trapz(interp_func_stress, fit_strain) 

        else:
            spline_energy = UnivariateSpline(sorted_strain, sorted_energy)
            spline_stress = UnivariateSpline(sorted_strain, sorted_stress)


            #derivative_stress = np.gradient(spline_stress, sorted_strain)
            derivative_stress = spline_stress.derivative(n=1)

            energy_density = sorted_strain * derivative_stress(sorted_strain)
            

            # Integrate the energy density using Simpson's rule
            total_energy_density = simps(energy_density, sorted_strain)

            # Calculate the strain energy using interpolation
            #total_strain_energy = interpolation.integral(min(strain_data), max(strain_data)) * 1e-9
            total_strain_energy = spline_energy.integral(min(sorted_strain), max(sorted_strain)) 

    except Exception as e:
        print(f"Using simple area under the stress-strain curve: {str(e)}")
        total_energy_density = np.trapz(sorted_stress, sorted_strain) 
    return abs(total_strain_energy), abs(total_energy_density)

