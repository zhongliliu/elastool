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
import numpy as np
from material_analysis import MaterialAnalysis
from calc_elastic_elate import visualize_elastic_anisotropy,print_elastic_tensor
from elate import ELATE


def sound_velocity(elastic_constants_dict, elastic_tensor, cwd, dimensional,latt_system,plotparameters,elateparameters):
    pos = vasp.read_vasp('%s/OPT/CONTCAR' % cwd)
    # The Planck's constant in m^2 Kg s^-1
    h = 6.626E-34
    # The reduced Planck's constant hbar in m^2 Kg s^-1
    hbar = 1.0546E-34
    # The Boltzmann constant in m^2 Kg s^-2 K-1
    k = 1.381E-23
    # The Avogadro's constant
    Na = 6.02E+23
    AMU_to_grams = 1.66053906660e-27  # Conversion factor from AMU to kg
    # The total volume of the system cell 
    volume = pos.get_volume()
    # The total mass of all atoms in the material, in AMU
    M = sum(pos.get_masses())
    # The number of atoms in the system (material)
    n = pos.get_global_number_of_atoms()
    mean_mass = M * AMU_to_grams/ n  # in kg
    #mean_mass = molar_mass/n/Na/1000 # in kg
    hardnessvalues = None
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
            K_Clarke = 0.87 * k * (mean_mass) **(-2/3) * (elastic_constants_dict['E']*1E+9)**(1/2) * rho**(1./6.)
            K_Cahill =  k/2.48 * (rho*Na/M/1E-3)**(2/3) *(2.*V_s +V_b)*1E3
            
            
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
                        
            elastic_constants_dict['V_s'] = V_s
            elastic_constants_dict['V_b'] = V_b
            elastic_constants_dict['V_p'] = V_p
            elastic_constants_dict['V_m'] = V_m
            elastic_constants_dict['T_D'] = T_D
            elastic_constants_dict['Q']   = kl 
            elastic_constants_dict['M_1']   = lame_1
            elastic_constants_dict['M_2']   = lame_2
            elastic_constants_dict['P']   = B/G  
            elastic_constants_dict['K_Clarke']   = K_Clarke
            elastic_constants_dict['K_Cahill']   = K_Cahill
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
            # Note you can use (C_ii-C_ij)/2        cell = pos.get_cell()
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
            # 1E-3*sqrt(Y_2D*(1-v_b)/rho_2D/(1+v_b)/(1-2*v_b)) 
            V_lb = 1E-3 * sqrt(abs(B_b + G_b) / rho_2D)
            V_sb = 1E-3 * sqrt(abs(G_b) / rho_2D)
            V_ma = (1.0 / 2.0 * (1.0 / V_sa**2.0 + 1 / V_la**2.0))**(-1.0 / 2.0) #(1.0 / 3.0 * (2.0 / V_sa**3.0 + 1 / V_la**3.0))**(-1.0 / 3.0)
            T_Da = (hbar / k) * (4. * pi * n / area / 1E-20)**(1. / 2.) * V_ma * 1E+3

            V_mb = (1.0 / 2.0 * (1.0 / V_sb**2.0 + 1 / V_lb**2.0))**(-1.0 / 2.0)
            T_Db = (hbar / k) * (4. * pi * n / area / 1E-20)**(1. / 2.) * V_mb * 1E+3
            layermodulus = 0.25*(c11+c22+2.*c12)
            hardnessvalues = calculate_hardness((B_a+B_b)/2., (G_a+G_b)/2., (Y_2Da+Y_2Db)/2., (v_a+v_b)/2., vol,dim=dimensional)
            # Recast Clarke formula for minimal thermal conductivity for 2D materials
            #kappa_2D = C * kb * (M**0.5) * (E_2D**0.5) * (sigma**1.5)
            K_Clarke = 0.435 * k/mean_mass * ((Y_2Da+Y_2Db)/2)**(1/2) * rho_2D**(1./2.)
            K_Cahill =  k/2.48 * (rho_2D*Na/M/1E-3) *(V_ma)*1E3 
            Rf = sqrt( (Y_2Da + Y_2Db)/2./rho_2D/height**2)/2./pi/1E+9  #Resonance_F = sqrt(E/rho/height**2)/2/pi, h here is the thickness
            v_avg = (v_a + v_b)/2; E_avg = (Y_2Da + Y_2Db)/2.
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
            elastic_constants_dict['V_m[10]'] = V_ma
            elastic_constants_dict['V_m[01]'] = V_mb
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
            elastic_constants_dict['P[10]']   = B_a/G_a  
            elastic_constants_dict['P[01]']   = B_b/G_b  
            #elastic_constants_dict['R_f'] = Rf
            elastic_constants_dict['K_Clarke']   = K_Clarke
            elastic_constants_dict['K_Cahill']   = K_Cahill
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


            elastic_tensor_Q2d = np.array([
                [c11, c12, 0, 0, 0, c16],
                [c12, c22, 0, 0, 0, c26],
                [0, 0, c66, 0, 0, 0],
                [0, 0, 0, c11, c12, 0],
                [0, 0, 0, c12, c22, 0],
                [c16, c26, 0, 0, 0, c66]
            ])
            
            

#            elastic_tensor_Q2d = np.array([
#                [c11, c12, c12, 0, 0, c16],
#                [c12, c11, c12, 0, 0, c26],
#                [c12, c12, c11, 0, 0, 0],
#                [0, 0, 0, c66, 0, 0],
#                [0, 0, 0, 0, c66, 0],
#                [c16, c26, 0, 0, 0, c66]
#            ])
            
            
                        
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
                K_Clarke = 0.87 * k * (mean_mass) **(-2/3) * (c33*1E+9)**(1/2) * mass_density**(1./6.)
                K_Cahill =  k/2.48 * (mass_density*Na/M/1E-3)**(2/3) *V*1E3
                
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




if __name__ == '__main__':
    elastic_constants_dict = {'B_vrh': 34.99, 'G_vrh': 17.20}
    cwd = '.'
    elastic_constants_dict = sound_velocity(elastic_constants_dict, cwd)
    print(elastic_constants_dict)
