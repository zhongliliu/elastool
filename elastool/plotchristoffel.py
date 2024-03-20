"""
  Elastool -- Elastic toolkit for zero and finite-temperature elastic constants and mechanical properties calculations

  Copyright (C) 2019-2024 by Chinedu Ekuma

  This program is free software; you can redistribute it and/or modify it under the
  terms of the GNU General Public License as published by the Free Software Foundation
  version 3 of the License.

  This program is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
  PARTICULAR PURPOSE.  See the GNU General Public License for more details.

  E-mail: cekuma1@gmail.com

"""

import christoffel
import numpy as np
import os

def run_christoffel_simulation(stiffness_tensor, density, dim,latticesystem, directions=None, num_theta=180, total_num_phi=720):
    """
    Run Christoffel simulation with given parameters.

    :param stiffness_tensor: The stiffness tensor for the material.
    :param density: The density of the material.
    :param dim: Dimension of the material ("2D" or "3D").
    :param directions: Optional. List of directions for simulation.
    :param num_theta: Optional. Number of theta divisions for the simulation grid.
    :param total_num_phi: Optional. Total number of phi divisions for the simulation grid.
    """
    
    # Default directions if not provided
    directions_file = 'phase_directions.txt'
    print("**********************************************************************************")
    print("Generating plots. All plots are in the property_plots folder!!!")
    if directions is None:
        if os.path.isfile(directions_file):
            with open(directions_file, 'r') as file:
                directions = [list(map(float, line.split())) for line in file.readlines()]
            print("Using user defined directions from directions.txt")
            
        elif dim == "3D":
            print("Setting default 3D directions")
            directions = [
                [1, 0, 0],
                [0, 1, 0],
                [0, 0, 1],
                [1, 1, 0],
                [0, 1, 1],
                [1, 0, 1],
                [1, 1, 1],
                [1, -1, -1],
                [1, -1, 1]
            ]
        elif dim == "2D":
            print("Setting default 2D directions")            
            directions = [
                [1, 0, 0],
                [0, 1, 0],
                [1, 1, 0],
                [-1, 1, 0],
                [1, -1, 0],
                [-1, -1, 0]
            ]
                                        
                
    directions = np.array(directions)

    # Read special directions if present
    if directions.size > 0:
        try:
            directions = directions.reshape((-1, 3))
        except ValueError:
            print("Error in reshaping directions")
            directions = []

    # Create Christoffel object
    chris = christoffel.Christoffel(stiffness_tensor, density, dim,latticesystem)
    print('Christoffel object created.\n')

    #Dump the data that has been read
    statusfile = open('sound.out', 'w')

    if dim == "2D":
        statusfile.write('Density: {0:.2e} kg/m^2\n\n'.format(chris.density))
        statusfile.write('Stiffness tensor in N/m:\n')
    elif dim == "3D":
        statusfile.write('Density: {0:.2f} kg/m^3\n\n'.format(chris.density))
        statusfile.write('Stiffness tensor in GPa:\n')
    else:
        raise ValueError("Invalid dimension")
            
            
    for i in range(6):
        statusfile.write('{0:8.2f} {1:8.2f} {2:8.2f} {3:8.2f} {4:8.2f} {5:8.2f}\n'.format(*(chris.stiffness2D[i])))


    if dim == "2D":
        statusfile.write('\nStiffness constant: {0:.2f} N/m\n'.format(chris.bulk))
        statusfile.write('Shear modulus: {0:.2f} N/m\n'.format(chris.shear))
    elif dim == "3D":
        statusfile.write('\nBulk modulus: {0:.2f} GPa\n'.format(chris.bulk))
        statusfile.write('Shear modulus: {0:.2f} GPa\n'.format(chris.shear))
    else:
        raise ValueError("Invalid dimension")
        
    statusfile.write('Isotropic primary velocity: {0:.2f} km/s\n'.format(chris.iso_P))
    statusfile.write('Isotropic secondary velocity: {0:.2f} km/s\n'.format(chris.iso_S))

    statusfile.close()

    #Prepare the output files
    SS_outfile = open('slow_secondary.dat', 'w')
    SS_outfile.write('# Slow secondary\n')
    SS_outfile.write('# Theta (rad) # Phi (rad) # Cube x y z ')
    SS_outfile.write('# Phase (km/s) # Phase (rel) # Phase pol x y z ')
    SS_outfile.write('# Group (km/s) # Group (rel) # Group x y z (km/s) ')
    SS_outfile.write('# Power flow angle (deg) # Enhancement factor\n')

    FS_outfile = open('fast_secondary.dat', 'w')
    FS_outfile.write('# Fast secondary\n')
    FS_outfile.write('# Theta (rad) # Phi (rad) # Cube x y z ')
    FS_outfile.write('# Phase (km/s) # Phase (rel) # Phase pol x y z ')
    FS_outfile.write('# Group (km/s) # Group (rel) # Group x y z (km/s) ')
    FS_outfile.write('# Power flow angle (deg) # Enhancement factor\n')

    P_outfile = open('primary.dat', 'w')
    P_outfile.write('# Primary\n')
    P_outfile.write('# Theta (rad) # Phi (rad) # Cube x y z ')
    P_outfile.write('# Phase (km/s) # Phase (rel) # Phase pol x y z ')
    P_outfile.write('# Group (km/s) # Group (rel) # Group x y z (km/s) ')
    P_outfile.write('# Power flow angle (deg) # Enhancement factor\n')

    all_out = [SS_outfile, FS_outfile, P_outfile]

    #Prepare min/max stuff
    maxvel = [1e-8, 1e-8, 1e-8]
    minvel = [1e8, 1e8, 1e8]
    maxdir = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
    mindir = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
    maxangle = [[0.0, 0.0], [0.0, 0.0], [0.0, 0.0]]
    minangle = [[0.0, 0.0], [0.0, 0.0], [0.0, 0.0]]

    #Get the isotropic velocities
    iso_vel = chris.get_isotropic()

    #print('Calculating phase velocities in default directions.')
   
    if os.path.isfile(directions_file):
        print('Calculating phase velocities using user defined directions.')
    else:
       print(f'Calculating phase velocities in default directions\nFor specific directions, modify phase_directions.dir ')
       
    #Loop over user-defined directions
    direction_file = open('directions.dat', 'w')
    direction_file.write('# The group velocity along various directions\n')
    for q in directions:
        chris.set_direction_cartesian(q)
        velocity = chris.get_phase_velocity()

        direction_file.write('Direction: ' + str(q[0]) + ' ' + str(q[1]) + ' ' + str(q[2]) + '\n')
        direction_file.write(' Primary: ' + str(velocity[2]) + ' km/s\n')
        direction_file.write(' Secondary 1: ' + str(velocity[1]) + ' km/s\n')
        direction_file.write(' Secondary 2: ' + str(velocity[0]) + ' km/s\n')
        direction_file.write(' Secondary avg: ' + str(0.5*(velocity[0]+velocity[1])) + ' km/s\n\n')
    direction_file.close()
    print('The directions.dat file has been successfully written.')

    print(f"Looping over surface along the following directions\n {directions}")
    #Loop over surface
    for theta in range(num_theta):
        ##print(str(theta+1) + '/' + str(num_theta))
        num_phi = total_num_phi
        #This commented bit scales num_phi by sin(theta) to get a constant sampling density
        #Sadly, GNUPlot doesn't play nice with grids that aren't rectangular
        #Still, if you want to cut your runtime down significantly, use the next two lines
        #sin_theta = np.sin(0.5*np.pi*theta / (num_theta - 1))
        #num_phi = int(np.ceil(sin_theta*total_num_phi)) + 1
        for phi in range(num_phi):
            theta_rad = 0.5*np.pi*theta/(num_theta-1)
            phi_rad = 2.0*np.pi*phi/(num_phi-1)
            chris.set_direction_spherical(theta_rad, phi_rad)
            q = chris.direction
            cubepos = q/np.linalg.norm(q, ord=np.inf)

            #Calculate everything and store in variables
            enhancefac = chris.get_enhancement()
            velocity = chris.get_phase_velocity()
            polarization = chris.get_eigenvec()
            group_vel = chris.get_group_velocity()
            group_abs = chris.get_group_abs()
            pf_angle = chris.get_powerflow()
 
            #Check for fastest and slowest velocities
            for i in range(3):
                if velocity[i] > maxvel[i]:
                    maxvel[i] = velocity[i]
                    maxdir[i] = q
                    maxangle[i] = [chris.theta, chris.phi]
                if velocity[i] < minvel[i]:
                    minvel[i] = velocity[i]
                    mindir[i] = q
                    minangle[i] = [chris.theta, chris.phi]
                rel_iso_vel = 100*(group_abs[i]/iso_vel[i]-1.0)
                if dim == "2D":
                    group_vel[i] *= 1E-6
                    group_abs[i] *= 1E-6
                    rel_iso_vel *=  1E-6
                    rel_iso_vel  *= rel_iso_vel*1E-9
                #And start writing data
                all_out[i].write('{0:.6f} {1:.6f}  {2: 8.6f} {3: 8.6f} {4: 8.6f}  '.format(chris.theta, chris.phi, *cubepos))
                all_out[i].write('{0: 10.6f}  {1: 8.6f}  {2: 8.6f} {3: 8.6f} {4: 8.6f}  '.format(velocity[i], 100*(velocity[i]/iso_vel[i]-1.0), *(polarization[i])))
                all_out[i].write('{0: 10.6f}  {1: 8.6f}  {2: 10.6f} {3: 10.6f} {4: 10.6f}  '.format(group_abs[i], rel_iso_vel, *(group_vel[i])))
                all_out[i].write('{0: 8.4f}  {1:.7g}\n'.format(180*pf_angle[i]/np.pi, enhancefac[i]))

        SS_outfile.write('\n')
        FS_outfile.write('\n')
        P_outfile.write('\n')
    SS_outfile.close()
    FS_outfile.close()
    P_outfile.close()

    print('Calculating anisotropy.')

    #Rounding gives nicer output
    maxdir = np.around(maxdir, 10)
    mindir = np.around(mindir, 10)

    #Write anisotropy data
    anisotropy_file = open('anisotropy.dat', 'w')
    anisotropy_file.write('Max Primary: ' + str(maxvel[2]) + ' km/s')
    anisotropy_file.write(' at direction (' + str(maxangle[2][0]) + ', ' + str(maxangle[2][1]) + ')')
    anisotropy_file.write(' - (' + str(maxdir[2][0]) + ', ' + str(maxdir[2][1]) + ', ' + str(maxdir[2][2]) + ')\n')
    anisotropy_file.write('Min Primary: ' + str(minvel[2]) + ' km/s')
    anisotropy_file.write(' at direction (' + str(minangle[2][0]) + ', ' + str(minangle[2][1]) + ')')
    anisotropy_file.write(' - (' + str(mindir[2][0]) + ', ' + str(mindir[2][1]) + ', ' + str(mindir[2][2]) + ')\n\n')

    anisotropy_file.write('Max Secondary 1: ' + str(maxvel[1]) + ' km/s')
    anisotropy_file.write(' at direction (' + str(maxangle[1][0]) + ', ' + str(maxangle[1][1]) + ')')
    anisotropy_file.write(' - (' + str(maxdir[1][0]) + ', ' + str(maxdir[1][1]) + ', ' + str(maxdir[1][2]) + ')\n')
    anisotropy_file.write('Min Secondary 1: ' + str(minvel[1]) + ' km/s')
    anisotropy_file.write(' at direction (' + str(minangle[1][0]) + ', ' + str(minangle[1][1]) + ')')
    anisotropy_file.write(' - (' + str(mindir[1][0]) + ', ' + str(mindir[1][1]) + ', ' + str(mindir[1][2]) + ')\n\n')

    anisotropy_file.write('Max Secondary 2: ' + str(maxvel[0]) + ' km/s')
    anisotropy_file.write(' at direction (' + str(maxangle[0][0]) + ', ' + str(maxangle[0][1]) + ')')
    anisotropy_file.write(' - (' + str(maxdir[0][0]) + ', ' + str(maxdir[0][1]) + ', ' + str(maxdir[0][2]) + ')\n')
    anisotropy_file.write('Min Secondary 2: ' + str(minvel[0]) + ' km/s')
    anisotropy_file.write(' at direction (' + str(minangle[0][0]) + ', ' + str(minangle[0][1]) + ')')
    anisotropy_file.write(' - (' + str(mindir[0][0]) + ', ' + str(mindir[0][1]) + ', ' + str(mindir[0][2]) + ')\n\n')

    anisotropy_file.write('Primary isotropic velocity: ' + str(chris.iso_P) + ' km/s\n')
    anisotropy_file.write('Secondary isotropic velocity: ' + str(chris.iso_S) + ' km/s\n')
    anisotropy_file.write('Primary Anisotropy: ' + str(100*(1.0 - minvel[2]/maxvel[2])) + '\n')
    anisotropy_file.write('Secondary Anisotropy: ' + str(100*(1.0 - minvel[0]/maxvel[1])) + '\n')

    anisotropy_file.close()

    print('Calculating phase velocities DONE!\nNow generating the phase plots')


# Example usage
# stiffness_tensor = [[...], [...], ...]  # Define your stiffness tensor heree
# run_christoffel_simulation(stiffness_tensor, density, "3D")

