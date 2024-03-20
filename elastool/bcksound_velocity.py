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
from math import pi, sqrt,exp
from numpy import cross, linalg,log, expm1
import numpy as np
from material_analysis import MaterialAnalysis
from calc_elastic_elate import visualize_elastic_anisotropy,print_elastic_tensor
from elate import ELATE, makePolarPlotPosNeg
from plotchristoffel import run_christoffel_simulation
from scipy.integrate import quad

def sound_velocity(elastic_constants_dict, elastic_tensor, cwd, dimensional,latt_system,plotparameters,elateparameters):
    pos = vasp.read_vasp('%s/OPT/CONTCAR' % cwd)
    # The Planck's constant in m^2 Kg s^-1
    h = 6.626E-34
    # The reduced Planck's constant hbar in m^2 Kg s^-1
    hbar = 1.0546E-34
    # The Boltzmann constant in m^2 Kg s^-2 K-1
    k = 1.381E-23
    # The Avogadro's constant
    Na = 6.023E+23
    AMU_to_grams = 1.66053906660e-27  # Conversion factor from AMU to kg
    # The total volume of the system cell 
    volume = pos.get_volume()
    # The total mass of all atoms in the material, in AMU
    M = sum(pos.get_masses())
    # The number of atoms in the system (material)
    n = pos.get_global_number_of_atoms()
    molar_mass = M * AMU_to_grams*Na # Kg/mol
    mean_mass = M * AMU_to_grams/ n  # in kg
    total_mass = M * 1E-3 / Na  # Kg
    hardnessvalues = None
    start = 10; interval = 10; end = 2000
    
    if dimensional == '3D':
        try:
            # The density in Kg/m^3
            rho = M * 1E-3 / Na / volume / 1E-30
            
            B = elastic_constants_dict['B_vrh']
            G = elastic_constants_dict['G_vrh']
            V_s = 1E-3 * G**0.5 * \
                (1E+9 / rho)**0.5
            V_b = 1E-3 * B**0.5 * \
                (1E+9 / rho )**0.5
            V_p = 1E-3 * (B + 4. * G / 3.)**0.5 * (1E+9 / rho )**0.5
            V_m = ((2. / V_s**3. + 1. / V_p**3) / 3.)**(-1. / 3.)
            T_D = (h / k) * (3. * n / 4. / pi * Na * rho / M / 1E-3)**(1. / 3.) * V_m * 1E+3

            
            if T_D > end:
               end = T_D +200
            
            
            analysis_instance_tmp = MaterialAnalysis(elastic_tensor,rho)
            dir_1 = np.array([1, 0, 0]) 
            dir_2 = np.array([1, 1, 0]) 
            dir_3 = np.array([0, 1, 0]) 
            dir_4 = np.array([1, 1, 1])
            
            lc_1 = analysis_instance_tmp.linear_compressibility(dir_1) 
            lc_2 = analysis_instance_tmp.linear_compressibility(dir_2)
            lc_3 = analysis_instance_tmp.linear_compressibility(dir_3)
            lc_4 = analysis_instance_tmp.linear_compressibility(dir_4)
            
            c11 = elastic_constants_dict['c11']
            c12 = elastic_constants_dict['c12']
            E = elastic_constants_dict['E']
            v = elastic_constants_dict['v']
            kl = (c11 + 8.*c12)/(7.*c11-2*c12)
            lame_1 = v * E/((1.+ v)*(1.-2.*v))
            lame_2 = E/2./(1.+v)
            

            cell = pos.get_cell()
            thickness  = np.linalg.norm(cell[2]) * 1E-10
            t_c = V_m*1.E3/thickness

            thermalcalculator = ThermalConductivityCalculator(dimensional)

            gamma = 3.*(1.+v)/(2.-3.*v)/2. #Grüneisen parameter or x_g = (V_p )/V_s; gamma_2 = 3./2.* (3.*x_g**2 - 4)/(x_g**2 +2) Technical Physics, 2009, Vol. 54, No. 9, pp. 1398–1401

            
            delta_atoms = volume**(1./3.)/n 
            #mean_mass
            K_Slack = thermalcalculator.slack_simple_model(M/n, T_D, delta_atoms, gamma, n, 300) 
            K_Clarke = thermalcalculator.clarke_model(n,elastic_constants_dict['E']*1E+9,rho, total_mass)  
            K_Cahill =  thermalcalculator.cahill_thermal_conductivity(n,volume*1E-30,V_p*1E+3, V_s*1E+3, V_b*1E+3) 

            K_Slack_int = thermalcalculator.slack_low_temp(V_m*1E3, T_D,T=300,t_c=t_c)
 
            temperatures, conductivities, conductivities_slack_simple, conductivities_cahill, capacities, entropy,debye_entropy,free_energy = thermalcalculator.compute_thermal_conductivity_over_temperature_range(V_m*1E3, T_D, M/n, \
                                                                                                      delta_atoms, gamma, n, volume*1E-30, start=start, interval=interval, end=end, t_c=t_c)
            
            filename = "thermodynamic_data.dat"
            data = np.column_stack((temperatures, conductivities, conductivities_slack_simple, conductivities_cahill, capacities,entropy,debye_entropy, free_energy))
            np.savetxt(filename, data, fmt='%.4f', delimiter=' ', \
            header='#T(K), Kl_SlackLow(W/mK), Kl_SlackSim(W/mK), Kl_CahillLow(W/mK), Cv(J/molK), S(J/molK), S_Debye(J/molK), H(eV)', comments='')

            #print("thermalcalculator ", thermalcalculator.slack_low_temp_gamma(V_m*1E3, T_D, volume*1E-30,gamma))
                        
            elastic_constants_dict['V_t'] = V_s
            elastic_constants_dict['V_l'] = V_b
            elastic_constants_dict['V_p'] = V_p
            elastic_constants_dict['V_D'] = V_m
            elastic_constants_dict['T_D'] = T_D
            elastic_constants_dict['Q']   = kl 
            elastic_constants_dict['M_1']   = lame_1
            elastic_constants_dict['M_2']   = lame_2
            elastic_constants_dict['S'] = gamma
            elastic_constants_dict['P']   = B/G  
            elastic_constants_dict['K_Clarke']   = K_Clarke
            elastic_constants_dict['K_Cahill']   = K_Cahill
            elastic_constants_dict['k_Cahill_lowT'] = thermalcalculator.cahill_thermalconduct(V_b*1E3, 300, T_D,n,volume*1E-30)
            elastic_constants_dict['k_Slack']   = K_Slack
            elastic_constants_dict['k_Slack_lowT']   = K_Slack_int
            elastic_constants_dict['h_Cv']   = thermalcalculator.constant_volume_hc(T_D, n,300)
            elastic_constants_dict['s_Entropy']   = thermalcalculator.entropy(T_D, n, 300)
            elastic_constants_dict['s_DebyeEntropy']   = thermalcalculator.debye_entropy(T_D, n, 300)
            elastic_constants_dict['C[100]']   = lc_1*1E3
            elastic_constants_dict['C[110]']   = lc_2*1E3
            elastic_constants_dict['C[010]']   = lc_3*1E3
            elastic_constants_dict['C[111]']   = lc_4*1E3
            elastic_constants_dict['C_avg']   =  (lc_1+lc_2+lc_3+lc_4)*1E3/4.
            
            elastic_constants_dict['D']   = ductility_test(B/G)

            np.savetxt("elastic_tensor.dat", elastic_tensor, fmt='%.4f', header='Elastic tensor in Voigt notation')
            rho_str = f"{rho:.8f}"  # Format to 8 decimal places
            data_to_save = np.array([rho_str, dimensional])
            
            np.savetxt("massdensity_dim.dat", [data_to_save], fmt='%s', header='Mass density in Kg/m^3, Dimension')

            vol = volume*1E-30
            hardnessvalues = calculate_hardness(B, G, E, v, vol,dim=dimensional)

            print_elastic_tensor(elastic_tensor=elastic_tensor, dim=dimensional, density=rho)

             
            if plotparameters[0]:
                
                analysis_instance = MaterialAnalysis(elastic_tensor,rho,plot=True,plotly=plotparameters[1])
                analysis_instance.plot_linear_compressibility_3D()
                #analysis_instance.plot_poisson_3Dprojection_3D()
                analysis_instance.plot_orientation_dependent_3D()
                analysis_instance.plot_contour_polar_3D()
                analysis_instance.heatmap_3D()
                analysis_instance.plot_moduli_heatmaps()
                analysis_instance.plot_directional_poisson_3D()
                analysis_instance.plot_poisson_3Dprojection_3D()
                analysis_instance.plot_shear_3Dprojection_3D()
                analysis_instance.plot_E_3Dprojection_3D()
                run_christoffel_simulation(elastic_tensor, rho, dimensional,latt_system)
                #makePolarPlotPosNeg(func=E)
                
                


            if "print" not in [item.lower() for item in elateparameters]:
                for elateparam in elateparameters:
                    if elateparam.lower() != "none" and elateparam.lower() != "print":
                        print("The followings are being plotted with Elate ", elateparameters)
                        visualize_elastic_anisotropy(elastic_tensor=elastic_tensor, density=rho, plot="3D", elastic_calc=elateparam)
                        visualize_elastic_anisotropy(elastic_tensor=elastic_tensor, density=rho, plot="3D_slice", elastic_calc=elateparam)
                        

        except BaseException:
            pass
    elif dimensional == '2D':
        try:
            c11 = elastic_constants_dict['c11']
            c12 = elastic_constants_dict['c12']
            if 'c22' in elastic_constants_dict:
                c22 = elastic_constants_dict['c22']
            else:
                c22 = elastic_constants_dict['c11']  # Lattice a=b
                #elastic_constants_dict['c66'] = (c11-c12)/2. 
                #elastic_constants_dict['c22'] = c11
            # Thomas et al., RSC Adv., 2018, 8, 27283
            Y_2Da = (c11 * c22 - c12**2) / c11
            Y_2Db = (c11 * c22 - c12**2) / c22
            v_a = c12 / c11
            v_b = c12 / c22
            # Note B is the in-plane stiffness (2D analogy of the bulk modulus)
            B_a = Y_2Da / 2.0 / (1 - v_a)
            G_a = Y_2Da / 2.0 / (1 + v_a)  # Note you can use (C_ii-C_ij)/2
            B_b = Y_2Db / 2.0 / (1 - v_b)
            # Note you can use (C_ii-C_ij)/2     
            G_b = Y_2Db / 2.0 / (1 + v_b)
            # PRB 94, 245420 (2016)
            # Calculate sound velocities in km/s

            cell = pos.get_cell()
            # The 2D density in Kg/m^2
            area = linalg.norm(cross(cell[0], cell[1]))
            vol = area*1E-20 # Note this is the area not volume
            rho_2D = M * 1E-3 / Na / area / 1E-20
            height = linalg.norm(cell[2])*1E-10 # m
            # 1E-3*sqrt(Y_2D*(1-v)/rho_2D/(1+v)/(1-2*v))
            V_la = 1E-3 * sqrt(abs(B_a + G_a) / rho_2D)
            V_sa = 1E-3 * sqrt(abs(G_a) / rho_2D)
            V_lb = 1E-3 * sqrt(abs(B_b + G_b) / rho_2D)
            V_sb = 1E-3 * sqrt(abs(G_b) / rho_2D)
            V_ma = (1.0 / 2.0 * (1.0 / V_sa**2.0 + 1 / V_la**2.0))**(-1.0 / 2.0) #(1.0 / 3.0 * (2.0 / V_sa**3.0 + 1 / V_la**3.0))**(-1.0 / 3.0)
            T_Da = (hbar / k) * (4. * pi * n / area / 1E-20)**(1. / 2.) * V_ma * 1E+3

            V_mb = (1.0 / 2.0 * (1.0 / V_sb**2.0 + 1 / V_lb**2.0))**(-1.0 / 2.0)
            T_Db = (hbar / k) * (4. * pi * n / area / 1E-20)**(1. / 2.) * V_mb * 1E+3
            layermodulus = 0.25*(c11+c22+2.*c12)
            hardnessvalues = calculate_hardness((B_a+B_b)/2., (G_a+G_b)/2., (Y_2Da+Y_2Db)/2., (v_a+v_b)/2., vol,dim=dimensional)
            
            if (T_Da +T_Db) > end:
               end = (T_Da +T_Db) +200            
            
            thermalcalculator = ThermalConductivityCalculator(dimension=dimensional)



            v_avg = (v_a + v_b)/2; E_avg = (Y_2Da + Y_2Db)/2.
            gamma = 3.*(1.+v_avg)/(2.-3.*v_avg)/2. #Grüneisen parameter
            delta_atoms = area**(1./2.)/n 

            K_Slack = thermalcalculator.slack_simple_model(M/n, (T_Da+T_Db)/2., delta_atoms, gamma, n, 300) 
            K_Cahill =  thermalcalculator.cahill_thermal_conductivity(n,area*1E-20,((V_la+V_lb)/2)*1E+3, ((V_sa+V_sb)/2)*1E+3) #k/2.48 * (rho_2D*Na/M/1E-3) *(V_ma)*1E3 
            K_Clarke = thermalcalculator.clarke_model(n,(Y_2Da+Y_2Db)/2.,rho_2D, total_mass)
            
            t_c = 1.0E13 #t_c is the scattering rate
            K_Slack_int = thermalcalculator.slack_low_temp(((V_ma+V_mb)/2.)*1E3, (T_Da+T_Db)/2.,T=300,t_c=t_c)

            temperatures, conductivities, conductivities_slack_simple, conductivities_cahill,capacities, entropy,debye_entropy,free_energy = thermalcalculator.compute_thermal_conductivity_over_temperature_range(((V_ma+V_mb)/2.)*1E3,\
                                                                                                     (T_Da+T_Db)/2., M/n, \
                                                                                                      delta_atoms, gamma, n, vol,start=start, interval=interval, end=end, t_c=t_c)
            
            #print("thermalcalculator ", thermalcalculator.slack_low_temp_gamma((V_ma+V_mb)/2., (T_Da+T_Db)/2., vol,gamma))
            
            filename = "thermodynamic_data.dat"
            data = np.column_stack((temperatures, conductivities, conductivities_slack_simple, conductivities_cahill, capacities,entropy,debye_entropy, free_energy))
            np.savetxt(filename, data, fmt='%.4f', delimiter=' ', \
            header='#T(K), Kl_SlackLow(W/mK), Kl_SlackSim(W/mK), Kl_CahillLow(W/mK), Cv(J/molK), S(J/molK), S_Debye(J/molK), H(eV)', comments='')


            Rf = sqrt( (Y_2Da + Y_2Db)/2./rho_2D/height**2)/2./pi/1E+9  #Resonance_F = sqrt(E/rho/height**2)/2/pi, h here is the thickness
            
            kl = (c11 + 8.*c12)/(7.*c11-2*c12)
            lame_1 = v_avg * E_avg/((1.+ v_avg)*(1.-2.*v_avg))
            lame_2 = E_avg/2./(1.+v_avg)

            analysis_instance_tmp = MaterialAnalysis(elastic_tensor,rho_2D)
            dir_1 = np.array([1, 0, 0]) 
            dir_2 = np.array([1, 1, 0]) 
            dir_3 = np.array([0, 1, 0]) 
            lc_1 = analysis_instance_tmp.linear_compressibility(dir_1) #(1.-elastic_constants_dict['v'])/elastic_constants_dict['E']  #
            lc_2 = analysis_instance_tmp.linear_compressibility(dir_2)
            lc_3 = analysis_instance_tmp.linear_compressibility(dir_3)

            
                        
            elastic_constants_dict['V_l[10]'] = V_la
            elastic_constants_dict['V_l[01]'] = V_lb
            elastic_constants_dict['V_s[10]'] = V_sa
            elastic_constants_dict['V_s[01]'] = V_sb
            elastic_constants_dict['V_D[10]'] = V_ma
            elastic_constants_dict['V_D[01]'] = V_mb
            elastic_constants_dict['T_D[10]'] = T_Da
            elastic_constants_dict['T_D[01]'] = T_Db
            elastic_constants_dict['Y[10]'] = Y_2Da
            elastic_constants_dict['Y[01]'] = Y_2Db
            elastic_constants_dict['G[10]'] = G_a
            elastic_constants_dict['G[01]'] = G_b
            elastic_constants_dict['v[10]'] = v_a
            elastic_constants_dict['v[01]'] = v_b
            elastic_constants_dict['B[10]'] = B_a
            elastic_constants_dict['B[01]'] = B_b
            elastic_constants_dict['L']   = layermodulus 
            elastic_constants_dict['Q']   = kl 
            elastic_constants_dict['M_1']   = lame_1
            elastic_constants_dict['M_2']   = lame_2
            elastic_constants_dict['S'] = gamma
            elastic_constants_dict['P[10]']   = B_a/G_a  
            elastic_constants_dict['P[01]']   = B_b/G_b  
            #elastic_constants_dict['R_f'] = Rf
            elastic_constants_dict['K_Clarke']   = K_Clarke
            elastic_constants_dict['K_Cahill']   = K_Cahill
            elastic_constants_dict['k_Cahill_lowT'] = thermalcalculator.cahill_thermalconduct(((V_la+V_lb)/2.)*1E3, 300, (T_Da+T_Db)/2.,n,vol)
            elastic_constants_dict['k_Slack']   = K_Slack
            elastic_constants_dict['k_Slack_lowT']   = K_Slack_int
            elastic_constants_dict['h_Cv']   = thermalcalculator.constant_volume_hc((T_Da+T_Db)/2., n,300)
            elastic_constants_dict['s_Entropy']   = thermalcalculator.entropy((T_Da+T_Db)/2., n, 300)
            elastic_constants_dict['s_DebyeEntropy']   = thermalcalculator.debye_entropy((T_Da+T_Db)/2., n, 300)
            elastic_constants_dict['C[100]']   = lc_1
            elastic_constants_dict['C[110]']   = lc_2
            elastic_constants_dict['C[010]']   = lc_3
            elastic_constants_dict['C_avg']   = (lc_1+lc_2+lc_3)/3.
            elastic_constants_dict['D']   = ductility_test((B_a/G_a + B_b/G_b)/2.) 
            # Build the elastic tensor for plotting
 

            cell = pos.get_cell()
            #length = linalg.norm(cell[2])
            c11 = elastic_tensor[0, 0] 
            c22 = elastic_tensor[1, 1]
            c12 = elastic_tensor[0, 1] 
            c66 = elastic_tensor[2, 2] 
            c16 = elastic_tensor[0, 2] 
            c26 = elastic_tensor[2, 1] 
            elastic_tensor_Q2d = np.zeros((6, 6))

            #Embedding method
            elastic_tensor_Q2d = np.array([
                [c11, c12, 0, 0, 0, c16],
                [c12, c22, 0, 0, 0, c26],
                [0, 0, c66, 0, 0, 0],
                [0, 0, 0, c11, c12, 0],
                [0, 0, 0, c12, c22, 0],
                [c16, c26, 0, 0, 0, c66]
            ])
                        
            np.savetxt("elastic_tensor.dat", elastic_tensor, fmt='%.4f', header='Elastic tensor in Voigt notation')

            rho_str = f"{rho_2D:.8f}"  # Format to 8 decimal places
            data_to_save = np.array([rho_str, dimensional])
            np.savetxt("massdensity_dim.dat", [data_to_save], fmt='%s', header='Mass density in Kg/m^2, Dimension')
            
            
            print_elastic_tensor(elastic_tensor=elastic_tensor_Q2d, dim=dimensional, density=rho_2D) 


            if plotparameters[0]:
                analysis_instance = MaterialAnalysis(elastic_tensor,rho_2D, plot=True,plotly=plotparameters[1])
                analysis_instance.plot_orientation_dependent_2D()
                analysis_instance.plot_contour_2D()
                #analysis_instance.plot_contour_polar_2D()
                analysis_instance.heatmap_2D()
                analysis_instance.plot_directional_poisson_2D()
                analysis_instance.plot_combined_directional_2D()
                #analysis_instance.plot_youngmodulus_directional_2D()
                run_christoffel_simulation(elastic_tensor, rho_2D, dimensional,latt_system)
                
                
            if "print" not in [item.lower() for item in elateparameters]:
                for elateparam in elateparameters:
                    if elateparam.lower() != "none" and elateparam.lower() != "print":
                        print("The followings are being plotted with Elate ", elateparameters)
                        visualize_elastic_anisotropy(elastic_tensor=elastic_tensor_Q2d, density=rho_2D, plot="2D", elastic_calc=elateparam)


        except BaseException:
            pass
    elif dimensional == '1D':
        try:
            cell = pos.get_cell()
            # We have1D structure is oriented along cell[2]
            length = linalg.norm(cell[2])*1E-10 # m
            area = pi*linalg.norm(cross(cell[0], cell[1]))*1E-20 #m^2
            rho_1D = M * 1E-3 / Na / length  / 1E-10 # kg/m
            Diameter = 2.*sqrt(area/pi)
            I_inertia = pi/64 * Diameter**4
            mass_density = M * 1E-3 / Na / volume / 1E-30 # kg/m^3
            
            if latt_system == 'Nanotube': 
                
                c33 = elastic_constants_dict['c33']
                c23 = elastic_constants_dict['c23']
                Tension = c33*area*1E+9 #N
                Y = c33
                v = -c23 / c33
                K = (c33-2*c23)/3.
                G = (c33-c23)/2.
                V = 1E-3* sqrt(abs(c33*1E+9) / mass_density)
                T_D = (h / k) * (3. * n / 4. / pi * Na * mass_density / M / 1E-3)**(1. / 3.) * V * 1E+3 
                Rf_eb =  sqrt( (Tension/rho_1D) + Y*I_inertia *1E+9/rho_1D/length**4)/2./pi/1E+9 # From Euler-Bernoulli beam theory: F_f = sqrt(T/mu + E*I_inertia/(mu*L^4))/2/pi
                #Rf_membrane = 2.4048/2./pi*Diameter/2 * sqrt(c33*1E+9/mass_density/3.5/1E-10) # More accuracte formula - See Ann. Phys. (Berlin), 1–18 (2014) SWCNT Diameter = thickness
                #Rf_plate = 10.21*3.5*1E-10 /4./pi/(Diameter/2.)**2*sqrt( Y*1E+9/3./mass_density/(1.- v**2) )

                dir_1 = np.array([1, 0, 0]) 
                dir_2 = np.array([1, 1, 0]) 
                dir_3 = np.array([0, 1, 0]) 
                dir_4 = np.array([1, 1, 1])  
                              
                analysis_instance_tmp = MaterialAnalysis(elastic_tensor,mass_density)
                #K_Clarke = 0.87 * k * (mean_mass) **(-2/3) * (c33*1E+9)**(1/2) * mass_density**(1./6.)
                #K_Cahill =  k/2.48 * (mass_density*Na/M/1E-3)**(2/3) *V*1E3


                thermalcalculator = ThermalConductivityCalculator(dimensional)


                K_Clarke = thermalcalculator.clarke_model(n,c33*1E+9,mass_density, total_mass)  
                K_Cahill =  thermalcalculator.cahill_thermal_conductivity(n,volume*1E-30,V*1E+3) 
                    
                
            
            
            
                            
                lc_1 = analysis_instance_tmp.linear_compressibility(dir_1) 
                
                lc_2 = analysis_instance_tmp.linear_compressibility(dir_2)
                lc_3 = analysis_instance_tmp.linear_compressibility(dir_3)
                lc_4 = analysis_instance_tmp.linear_compressibility(dir_4)
                
                elastic_constants_dict['Y'] = Y
                #elastic_constants_dict['B'] = K # Bulk modulus
                elastic_constants_dict['G'] = G #Shear modulus
                #elastic_constants_dict['v'] = v
                elastic_constants_dict['V'] = V
                elastic_constants_dict['T_D'] = T_D
                #elastic_constants_dict['P']   = K/G   
                elastic_constants_dict['K_Clarke']   = K_Clarke
                elastic_constants_dict['K_Cahill']   = K_Cahill
                elastic_constants_dict['C[100]']   = lc_1*1E3
                elastic_constants_dict['C[110]']   = lc_2*1E3
                elastic_constants_dict['C[010]']   = lc_3*1E3
                elastic_constants_dict['C[111]']   = lc_4*1E3
                elastic_constants_dict['C_avg']   =  (lc_1+lc_2+lc_3+lc_4)*1E3/4.
                
                #elastic_constants_dict['D']   = ductility_test(K/G)  
                #elastic_constants_dict['RF_Euler-Bernoulli'] = Rf_eb
                #elastic_constants_dict['RF_membrane'] = Rf_membrane
                #elastic_constants_dict['RF_membrane+plate'] = sqrt(Rf_membrane**2 + Rf_plate**2)/1E+9

                np.savetxt("elastic_tensor.dat", elastic_tensor, fmt='%.4f', header='Elastic tensor in Voigt notation')


                rho_str = f"{rho_1D:.8f}"  # Format to 8 decimal places
                data_to_save = np.array([rho_str, dimensional])
            
                np.savetxt("massdensity_dim.dat", [data_to_save], fmt='%s', header='Mass density in Kg/m^3, Dimension')
            
            
                hardnessvalues = calculate_hardness(None, None, None, None, None,dim=dimensional)
        except BaseException:
            pass
    return elastic_constants_dict,hardnessvalues
    
                                                                              
                                                                         
def ductility_test(ratio):
    if(ratio > 1.75):
        return "Material is ductile"
    else:
        return "Material is brittle"



def calculate_hardness(*args,dim='3D'):
    """
    Calculate hardness using various methods.

    Parameters:
    -----------
    bulk_modulus : float
        The bulk modulus value.
    shear_modulus : float
        The shear modulus value.
    youngs_modulus : float
        Young's modulus value.
    poisson_ratio : float
        Poisson's ratio value.

    vol : float
        Volume.
    Returns:
    -------
    tuple
        The hardness values calculated by 6 different methods:
        [H1a, H1b, and H1c] Correlation between hardness and elastic moduli of the covalent crystals. Jiang, et al. (2011).
        [H2] Computational alchemy: the search for new superhard materials. Teter (1998).
        [H3] Mechanical and electronic properties of B12-based ternary crystals of orthorhombic phase. Jiang et al. (2010).
        [H4] Theoretical investigation on the transition-metal borides with Ta3B4-type structure: A class of hard and refractory materials. Miao et al. (2011).
        [H5] Modeling hardness of polycrystalline materials and bulk metallic glasses. Chen et al. (2011).
        [H6] A model of hardnesss and fracture toughness, Journal of Appl. Phys. 126, 125109 (2019).
        [K1] Simpe and accuracte model of fracture toughness of solids, J. Appl. Phys. 125, 065105 (2019)
        [K2] A model of hardnesss and fracture toughness, J. Appl. Phys. 126, 125109 (2019).
    """
    bulk_modulus, shear_modulus, youngs_modulus, poisson_ratio, vol = args
    B = bulk_modulus
    G = shear_modulus
    Y = youngs_modulus
    v = poisson_ratio
    k = G / B
    chi_v = (1-8.5*v + 19.5*v**2)/(1.-7.5*v +12.2*v**2 +19.6*v**3)
    k_v = (1.-13.7*v + 48.6*v**2)/(1.-15.2*v + 70.2*v**2 - 81.5*v**3)
    H1a = (1 / 6.78) * G
    H1b = (1 / 16.48) * Y
    H1c = 0.0963 * B
    H2 = (1 / 15.76) * Y
    H3 = (0.1769 * G) - 2.899
    H4 = ((1 - 2 * v) * B) / (6 * (1 + v))
    H5 = 2 * ((k * k * G) ** 0.585) - 3
    H6 = 2 * ((  (9*Y*(1.-2.*v)**2)/(8*(1+v)**3))**0.585) - 3     
    H7 = 0.096 * chi_v * Y  
    if dim == "3D":    
        K1 = vol**(1./6.)*G*(B/G)**(1./2)
        K2 = 8840**(-1./2.) * vol**(1./6.) *(k_v * Y ) **(3./2.)
        K3 = 2**(-1./2.) * vol**(1./6.) *(k_v * Y ) **(3./2.)
    elif dim == "2D":    
        K1 = vol**(1./6)*G*(B/G)**(1./2)
        K2 = 8840**(-1./2.) * vol**(1./6.) *(k_v * Y ) **(3./2.)
        K3 = 2**(-1./2.) * vol**(1./6.) *(k_v * Y ) **(3./2.)
    elif dim == "1D":
        return "Hardness for 1D tubular structures not yet implemented"
    return H1a, H1b, H1c, H2, H3, H4, H5, H6, H7, K1, K2, K3




class ThermalConductivityCalculator:
    def __init__(self,dimension="3D"):
        self.k = 1.380649e-23  # Boltzmann constant, J/K
        self.hbar = 1.0545718e-34  # Reduced Planck's constant, J·s
        self.Na   = 6.023E+23
        self.dimension = dimension

    def callaway_integrand(self, x, t_ph):
        return ( x ** 3 * exp(x)) /t_ph/ expm1(x)**2 if  self.dimension =="2D" else  ( x ** 4 * exp(x)) /t_ph/ expm1(x)**2 



#    def callaway_integrand_gamma(self, x):
#        return ( x ** 3 * exp(x))/ expm1(x)**2 if  self.dimension =="2D" else  ( x ** 4 * exp(x)) / expm1(x)**2 


    def slack_integrand(self, x, t_ph):
        return x  * t_ph  if  self.dimension =="2D" else  x ** 2  * t_ph 

    def cahill_integrand(self, x):
        """
        Integrand function for Cahill thermal conductivity.

        Args:
            x: (float) (hbar * omega) / (k * T)

        Returns:
            (float) Value of the integrand at x
        """
        return (x **2 * exp(x)) / (expm1(x) ** 2) if self.dimension =="2D" else (x ** 3 * exp(x)) / (expm1(x) ** 2)
         
    def cahill_integrand_summation(self, v_i, T, theta):
        """
        Calculate the summation term for the Cahill thermal conductivity integrand model.

        Args:
            v_i: (float) sound velocity for the acoustic mode i (in m/s)
            T: (float) absolute temperature (in K)
            theta: (float) Debye temperature (in K)

        Returns:
            (float) Summation term for one acoustic mode i
        """
        integrand = lambda x: self.cahill_integrand(x)
        integral_result, error = quad(integrand, 0, theta / T)
        
        #integral_result, _ = quad(self.cahill_integrand, 0, theta / T)
        return v_i * (T / theta) ** 2 * integral_result
                 
         
       
    def slack_low_temp(self, v_m, theta, T=300, t_c=1E+12):
        """
        Calculate Slack thermal conductivity.

        Args:
            v_m (float): Sound velocity.
            theta (float): Debye temperature.
            T (float): Temperature
            t_c (float): Scattering rate.

        Returns:
            float: Slack thermal conductivity (in SI units, i.e. W/(m-K))
        """
    
        integrand = lambda x: self.callaway_integrand(x, t_c)
        integral_result, error = quad(integrand, 0, theta / T)
        return (self.k / (2 * pi ** 2 * v_m)) * (self.k * T / self.hbar) ** 3 * integral_result



#    def slack_low_temp_gamma(self, v_m, theta, Vol,gamma, T=300):
#        """
#        Calculate Slack thermal conductivity.

#        Args:
#            v_m (float): Sound velocity.
#            theta (float): Debye temperature.
#            T (float): Temperature
#            gamma (float): G parameter
#            Vol (float): volume

#        Returns:
#            float: Slack thermal conductivity (in SI units, i.e. W/(m-K))
#        """
    
#        integrand = lambda x: self.callaway_integrand_gamma(x)
#        integral_result, error = quad(integrand, 0, theta / T)
#        return (self.k/ (2 * pi * v_m*gamma**2)) * (theta / T) ** 2 * integral_result if self.dimension=="2D" else self.k/(2 * pi**2 * v_m*gamma**2)*(self.k*T/self.hbar)**2  * integral_result
        

    def slack_high_temp(self, v_m, theta, T=300, t_c=1E+12):
        """
        Calculate Slack thermal conductivity.

        Args:
            v_m (float): Sound velocity.
            theta (float): Debye temperature.
            T (float): Temperature
            t_c (float): Scattering rate.

        Returns:
            float: Slack thermal conductivity (in SI units, i.e. W/(m-K))
        """
    
        integrand = lambda x: self.slack_integrand(x, t_c)
        integral_result, error = quad(integrand, 0, theta / T)
        return (self.k / (2 * pi ** 2 * v_m)) * (self.k * T / self.hbar) ** 3 * integral_result


    def cahill_thermalconduct(self, velocities, T, theta,n,V):
        """
        Calculate Cahill thermal conductivity using the integrand model.

        Args:
            velocities: (list of float) sound velocities for each acoustic mode (in m/s)
            T: (float) absolute temperature (in K)
            theta: (float) Debye temperature (in K)
            Ref: 10.1038/nature13184

        Returns:
            (float) Cahill thermal conductivity (in W/(m*K))
        """
        #if self.dimension == "2D":
        #    n_d = (4 * pi * (n * self.Na / V)) ** (1. / 2.)
        #else:
        #    n_d = ((n*self.Na /V)) ** (2. / 3.)
        cahill_sum = self.cahill_integrand_summation(velocities, T, theta)
        
        return (pi / 6) ** (1. / 3.) * self.k * ((n/V))  * cahill_sum  if self.dimension == "2D" else  (pi / 6) ** (1. / 3.) * self.k * ((n /V)) ** (2. / 3.)  * cahill_sum
         
                
    def debye_model(self, M, E, rho):
        """
        Calculate Debye thermal conductivity based on the dimensionality of the material.

        Args:
            M (float): Molecular mass (in atomic mass units).
            E (float): Young's modulus (in Pascals).
            rho (float): Density (in kg/m^3 for 3D or kg/m^2 for 2D).
            dim (str): Dimensionality of the material ('2D' or '3D', or '1D').

        Returns:
            float: Debye thermal conductivity (in SI units, i.e., W/(m·K)).
        """
        return 2.489e-11 * (M ** (-1/3)) * (E ** 0.5) * (rho ** (-1/2)) if self.dimension == "2D" else 2.489e-11 * (M ** (-1/3)) * (E ** 0.5) * (rho ** (-1/6))

    def constant_volume_hc(self, t_d, n,temp):
        """
        Calculate the molar heat capacity at constant volume.
        
        Args:
            N (int): Number of atoms in the system = n*Na.
            t_d (float): Debye temperature.
            temp (float): Current temperature.
        
        Returns:
            float: Molar heat capacity at constant volume in J/(mol·K).
        """
        
        def integrand(x):
            exm1 = expm1(x)
            if exm1 == 0:
                return 0
            else:
                return (x**3 * exp(x)) / exm1**2 if self.dimension == "2D" else  (x**4 * exp(x)) / exm1**2

        t_ratio = temp / t_d
        upper_limit = min(t_d / temp, 10)
        
        #integrand = lambda x: (x**4 * exp(x) / (exp(x) - 1)**2)
        if self.dimension == "2D":
            integral_result, _ = quad(integrand, 0, upper_limit)
            t_ratio = t_ratio**2
            const = 3
        else:
            integral_result, _ = quad(integrand, 0, upper_limit)
            t_ratio = t_ratio**3
            const = 9
        
        c_v = const * self.Na*n*self.k * t_ratio * integral_result #quad(integrand, 0, 1/t_ratio)[0]
        return c_v


    def debye_entropy(self, t_d, n, temp):
         
        # Debye integral for entropy => \[ S(T) = 3Nk\left(-\ln\left[1 - e^{-\frac{\Theta_D}{T}}\right] + 4\frac{T^3}{\Theta_D^3}\int_{0}^{\frac{\Theta_D}{T}}\frac{x^3}{e^x - 1}dx\right) \]
        # J/(mol·K)
        def integrand(x):
            return (x**2) / expm1(x) if self.dimension == "2D" else (x**3) / expm1(x)

                    
        if self.dimension == "2D":
            integral_result, _ = quad(integrand, 0, t_d/temp)
            localterm = 2 * (temp/t_d)**3
            const = 1.0
        else:
            integral_result, _ = quad(integrand, 0, t_d/temp)
            localterm = 4 * (temp/t_d)**3
            const = 3.0

        # Calculating entropy
        log_arg = 1. - exp(-t_d/temp)
        if log_arg <= 0:
            log_arg = 1e-10  # A small positive number to avoid log(0) or log(negative)
        s = const * n * self.Na * self.k * ( localterm * integral_result - log(log_arg) )
        return s

    def entropy(self,Theta_D, n, T):
        def cv_over_t(T_prime):
            if T_prime <= 0:
                return 0
            else:
                return self.constant_volume_hc(Theta_D, n,T_prime)/ T_prime   
        
        # Integrate C_v/T' from a small positive value near 0 to T to avoid division by zero
        S, _ = quad(cv_over_t, 1e-3, T)#, limit=100) 
        return S
            
    
    
    def slack_simple_model(self, M, theta, v_a, gamma, n, T): 
        """
        Calculate the simple Slack thermal conductivity

        References
            - DOI: 10.1007/0-387-25100-6_2 (Title: "High Lattice Thermal Conductivity Solids")
            - DOI: 10.1002/adfm.201600718 (Title: "Minimum Thermal Conductivity in Weak Topological Insulators with
            Bismuth-Based Stack Structure")

        Args:
            A_0: (float) constant, depends on Gruneisen parameter
            M: (float) average atomic mass
            theta: (float) Debye temperature (K)
            v_a: (float) (v_a)**3 is the volume per atom (Angstroms)
            gamma: (float) Gruneisen parameter
            n: (int) number of atoms in primitive cell
            T: (float) absolute temperature (K)

        Returns: (float) Slack thermal conductivity (in SI units, i.e. W/(m.K))

        """
        if self.dimension == '2D':
            powr = 1.
        else:
            powr = 2./3.
        A_0 = 2.43*1.0E-8/(1.0 -0.514/gamma +0.228/gamma**2 ) # See Phys. Rev. 137, A128 (1965)   
        return (A_0 * M * theta ** 3 * v_a) / (gamma * n ** (powr) * T)


    def cahill_thermal_conductivity(self, n, V,v_l, *v_ts):
        """
        Calculate Cahill thermal conductivity.
        Args:
            n: (int) number of atoms in the unit cell
            V: (float) unit cell volume (in SI units, i.e., m^3)
            v_l: (float) longitudinal sound velocity (in SI units, i.e., m/s)
            *v_ts: (float) variable number of transverse sound velocities (in SI units, i.e., m/s)
                   Accepts either 1 or 2 velocities for 2D or 3D materials, respectively.

        Returns: (float) Cahill thermal conductivity (in W/(m·K))
        """
        # Sum the longitudinal and all transverse sound velocities
        #total_velocity = v_l + sum(v_ts)
        if self.dimension == "1D":
            total_velocity = v_l  # Only longitudinal velocity is considered for 1D
        else:
            total_velocity = v_l + sum(v_ts) 
            
        if self.dimension =="2D":
            factor = (n / V)
        else:
            factor = (n / V) ** (2. / 3.)
        thermal_conductivity = 1./2. * ((pi / 6) ** (1. / 3.)) * self.k * factor * total_velocity
        # k/2.48 * (rho_2D*Na/M/1E-3) *(V_ma)*1E3
        return thermal_conductivity
        

    def clarke_model(self, n, E, rho, m):
        """
        Calculate Clarke thermal conductivity for both 3D and 2D materials.

        Args:
            n: (int) number of atoms in primitive cell
            E: (float) Young's modulus (in SI units, i.e., Kgm(s)^(-2))
            rho: (float) density (in SI units, i.e., Kg/m^3 for 3D or kg/m^2 for 2D)
            m: (float) total mass per unit cell (in SI units, i.e., Kg)
            dim: (str) dimensionality of the material ('3D' or '2D')

        Returns:
            (float) Clarke thermal conductivity (in SI units, i.e., W/(m·K))
        """
        if self.dimension == '2D':
            return 0.87 * self.k/m * (E ** (1. / 2.)) * (rho ** (1. / 2.))
        else:
            return 0.87 * self.k /(m ** (2. / 3.)) * (E ** (1. / 2.)) * (rho ** (1. / 6.))
            


    def internalenergy(self, t_d, n,temp):
        """
        Calculate the Internal energy.
        
        Args:
            N (int): Number of atoms in the system = n*Na.
            t_d (float): Debye temperature.
            temp (float): Current temperature.
        
        Returns:
            float: Internal energy in J.
        """
        
        def integrand(x):
            exm1 = expm1(x)
            if exm1 == 0:
                return 0
            else:
                return x**3 / exm1


        def integrand_2D(x):
            exm1 = expm1(x)
            if exm1 == 0:
                return 0
            else:
                return x**2 / exm1
                
            
        t_ratio = temp / t_d
        upper_limit = min(t_d / temp, 10)
        
 
        if self.dimension == "2D":
            integral_result, _ = quad(integrand_2D, 0, upper_limit)
            t_ratio = t_ratio**2
            const = 1
        else:
            integral_result, _ = quad(integrand, 0, upper_limit)
            t_ratio = t_ratio**3
            const = 9
        
        U = const * self.Na*n*self.k * t_ratio * integral_result
        return U



    def enthalpy(self, t_d, n,temp):
      return (self.internalenergy(t_d, n,temp) - temp*self.entropy(t_d, n, temp))*n/self.Na/1.602E-19

                            
    def compute_thermal_conductivity_over_temperature_range(self, v_m, theta, M, v_a, gamma, n, vol, start, interval, end, t_c=1E+12):
        temperatures = np.arange(start, end + interval, interval)
        thermal_conductivities = []
        thermal_conductivities_slack_simple = []
        themal_cal_cahill = []
        heat_capacity = []
        entropies =  []
        debye_entropies = []
        free_energies = []

        for T in temperatures:
            conductivity = self.slack_low_temp(v_m, theta, T, t_c)
            conductivity_slack_simple = self.slack_simple_model(M, theta, v_a, gamma, n, T)
            conductivity_cahill = self.cahill_thermalconduct(v_m, T, theta,n,vol)
            capacity = self.constant_volume_hc(theta, n,T)
            entropy = self.entropy(theta, n, T)
            debye_entropy = self.debye_entropy(theta, n, T)
            free_energy = self.enthalpy(theta,n,T)
            heat_capacity.append(capacity)
            thermal_conductivities.append(conductivity)
            thermal_conductivities_slack_simple.append(conductivity_slack_simple)
            entropies.append(entropy)
            debye_entropies.append(debye_entropy)
            free_energies.append(free_energy)
            themal_cal_cahill.append(conductivity_cahill)

        return temperatures, thermal_conductivities, thermal_conductivities_slack_simple, themal_cal_cahill, heat_capacity, entropies, debye_entropies,free_energies
 
       
         
        

if __name__ == '__main__':
    elastic_constants_dict = {'B_vrh': 34.99, 'G_vrh': 17.20}
    cwd = '.'
    elastic_constants_dict = sound_velocity(elastic_constants_dict, cwd)
    print(elastic_constants_dict)
