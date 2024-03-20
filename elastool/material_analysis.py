import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from numpy import sin,cos,pi,tan,vstack, linalg
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from matplotlib.colors import PowerNorm
import matplotlib.colors as mcolors
from mpl_toolkits.mplot3d import Axes3D


class MaterialAnalysis:

    def __init__(self, elastic_tensor, density, plot=False, plotly=False):
        if not (isinstance(elastic_tensor, np.ndarray) and 
                (elastic_tensor.shape == (3, 3) or elastic_tensor.shape == (6, 6))):
            raise ValueError("elastic_tensor must be a 3x3 matrix (for 2D) or a 6x6 matrix (for 3D).")

        
        self.Cdim = elastic_tensor
        self.plot = plot
        self.rho_mass = density
        det = np.linalg.det(self.Cdim)
        if det != 0:
        #    # Invert the matrix
           self.Cs = linalg.inv(elastic_tensor)
        else:
            # If the matrix is singular, use pseudo-inverse
           self.Cs = np.linalg.pinv(elastic_tensor)
                
        
        self.plotly = plotly
        #if elastic_tensor.shape == (6, 6):
        self.C_tensor = self.voigt_to_tensor(self.Cdim)  # Convert to 3x3x3x3 tensor for 3D materials
        self.Cs_tensor = self.voigt_to_tensor(self.Cs)

    @property
    def v_2D(self):
        """
        return v_2D=C12/C22
        """
        return self.Cdim[0][1]/self.Cdim[1][1]

    @property
    def d1(self):
        """
        return d1=C11/C22+1-(C11*C22-C12**2)/C22/C66;
        """
        return self.Cdim[0][0]/self.Cdim[1][1]+1 - \
               (self.Cdim[0][0]*self.Cdim[1][1]-self.Cdim[0][1]**2)/ \
               self.Cdim[1][1]/self.Cdim[2][2]

    @property
    def d2(self):
        """
        return d2=-(2*C12/C22-(C11*C22-C12**2)/C22/C66);
        """
        return -1*(2*self.Cdim[0][1]/self.Cdim[1][1]-\
                (self.Cdim[0][0]*self.Cdim[1][1]-self.Cdim[0][1]**2)/ \
                 self.Cdim[1][1]/self.Cdim[2][2])


    def print_Cs(self):
        print("Compliance matrix (Cs):")
        print(self.Cs)


    @property
    def d3(self):
        """
        return d3 =C11/C22
        """
        return self.Cdim[0][0]/self.Cdim[1][1]

    @property
    def Y_2D(self):
        """
        return Y_2D = C11*C22-C12**2)/C22
        """
        return (self.Cdim[0][0]*self.Cdim[1][1]-self.Cdim[0][1]**2)/self.Cdim[1][1]

    def compute_E_2D_polar(self, theta):
        denominator = cos(theta)**4 + self.d2 * cos(theta)**2 * sin(theta)**2 + self.d3 * sin(theta)**4
        E_theta = self.Y_2D / denominator
        return E_theta

    def compute_V_2D_polar(self, theta):
        numerator = self.v_2D * cos(theta)**4 - self.d1 * cos(theta)**2 * sin(theta)**2 + self.v_2D * sin(theta)**4
        denominator = cos(theta)**4 + self.d2 * cos(theta)**2 * sin(theta)**2 + self.d3 * sin(theta)**4
        V_theta = numerator / denominator
        return V_theta



    def compute_G_2D_polar(self, theta):
        S11, S12, S16, S22, S26, S66 = self.Cs[0, 0], self.Cs[0, 1], self.Cs[0, 2], self.Cs[1, 1], self.Cs[1, 2], self.Cs[2, 2]
        
        cos_theta = cos(theta)
        sin_theta = sin(theta)
        cos2_theta = cos_theta**2
        sin2_theta = sin_theta**2
        cos3_theta = cos_theta**3
        sin3_theta = sin_theta**3
        cos4_theta = cos_theta**4
        sin4_theta = sin_theta**4

        # Calculate 1/(4*G(theta))
        G_inv_theta = (S11 + S22 - 2 * S12) * cos2_theta * sin2_theta + \
                      0.25 * S66 * (cos4_theta + sin4_theta - 2 * sin2_theta * cos2_theta) + \
                      S16 * (sin3_theta * cos_theta - cos3_theta * sin_theta) + \
                      S26 * (cos3_theta * sin_theta - sin3_theta * cos_theta)

        return 1 / (4 * G_inv_theta)

        
    def compute_K_2D_polar(self, theta):
        K_theta = (self.Cdim[0][0] + self.Cdim[1][1] - 2*self.Cdim[0][1] + 
                   (self.Cdim[0][0] - self.Cdim[0][1] - self.Cdim[1][1]) * np.cos(2*theta)) / 2
        return K_theta



    def uvw_to_l(self):
        """
        Converts the [uvw] direction to the six-component vector l for cubic symmetry.
    
        Parameters:
        - u, v, w: Crystallographic direction [uvw]
    
        Returns:
        - l: six-component direction vector
        """
        npoints = 360
        theta = np.linspace(0, np.pi, npoints)
        phi = np.linspace(0, 2 * np.pi, npoints)

        theta, phi = np.meshgrid(theta, phi)

        u = np.sin(theta) * np.cos(phi)
        v = np.sin(theta) * np.sin(phi)
        w = np.cos(theta)
    
        # Normalize the [uvw] direction
        #norm = np.sqrt(u**2 + v**2 + w**2)
        #u /= norm
        #v /= norm
        #w /= norm

        # Convert [uvw] to the six-component vector l
        l = np.array([
            u**2,
            v**2,
            w**2,
            2*u*v,
            2*v*w,
            2*u*w
        ]).transpose(1, 2, 0)
    
        return l


    def adjust_colormap_brightness(self,cmap_name, scale_factor):
        cmap = plt.get_cmap(cmap_name)
        colors = cmap(np.arange(cmap.N))


    def compute_E_general(self,l):
        """
        Compute the orientation-dependent modulus E_uvw based on the compliance matrix S and direction [uvw].
    
        Parameters:
        - S: 6x6 compliance matrix
        - u, v, w: Crystallographic direction [uvw]
    
        Returns:
        - E_uvw: orientation-dependent modulus
        """
        #l = self.uvw_to_l()
        E_inv = sum(self.Cs[i, j] * l[i] * l[j] for i in range(6) for j in range(6))
        return 1.0 / E_inv



#    def compute_G_general(self,l):
#        """
#        Compute the orientation-dependent shear modulus G_uvw based on the compliance matrix S and direction l.
#    
#        Returns:
#        - G_uvw: orientation-dependent shear modulus
#        """
#        #l = self.uvw_to_l()
#        G_inv = sum(l[i] * l[j] * l[k] * l[m] * self.Cs[i, j, k, m] for i in range(6) for j in range(6) for k in range(6) for m in range(6))
#        print(f"For l: {l}, G_inv is: {1./G_inv}")
#        return 1.0 / G_inv


    def compute_G_general(self, l):
        # Mapping from Voigt notation to 4D tensor indices
        voigt_to_tensor_mapping = {
            0: (0, 0),
            1: (1, 1),
            2: (2, 2),
            3: (1, 2),
            4: (0, 2),
            5: (0, 1)
        }

        G_inv = 0
        for i in range(6):
            for j in range(6):
                tensor_i = voigt_to_tensor_mapping[i]
                tensor_j = voigt_to_tensor_mapping[j]
                G_inv += l[tensor_i[0]] * l[tensor_i[1]] * l[tensor_j[0]] * l[tensor_j[1]] * self.Cs[i, j]

        #print(f"For l: {l}, G_inv is: {G_inv}")
    
        if G_inv == 0:
          #  print("G_inv is zero!")
            return float('inf')  # or handle this case appropriately

        return 1.0 / G_inv

    def compute_V_general(self,l):
        """
        Compute the orientation-dependent Poisson's ratio nu_uvw based on E_uvw and G_uvw.
    
        Returns:
        - nu_uvw: orientation-dependent Poisson's ratio
        """
        E_uvw = self.compute_E_general(l)
        G_uvw = self.compute_G_general(l)
        return E_uvw / (2 * G_uvw) - 1


    def compute_K_general(self,l):
        """
        Compute the orientation-dependent bulk modulus K_uvw based on the compliance matrix S and direction l.
    
        Returns:
        - K_uvw: orientation-dependent bulk modulus
        """
       # l = self.uvw_to_l()
        K_inv = sum(l[i] * l[j] * self.Cs[i, j] for i in range(6) for j in range(6))
        return 1.0 / K_inv


    def voigt_to_tensor(self,Cscm):
        """
        Convert the Voigt matrix of elastic constants to a 3x3x3x3 tensor.
        Handles both 6x6 matrices (for 3D materials) and 3x3 matrices (for 2D materials).
        """

        C_tensor = np.zeros((3, 3, 3, 3))

        if Cscm.shape == (6, 6):  # 3D material
            voigt_notation = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]
            for i in range(6):
                for j in range(6):
                    indices = (voigt_notation[i][0], voigt_notation[i][1], voigt_notation[j][0], voigt_notation[j][1])
                    C_tensor[indices] = Cscm[i, j]
                    if i != j:
                        # Symmetric indices
                        C_tensor[voigt_notation[i][1], voigt_notation[i][0], voigt_notation[j][0], voigt_notation[j][1]] = Cscm[i, j]
                        C_tensor[voigt_notation[i][0], voigt_notation[i][1], voigt_notation[j][1], voigt_notation[j][0]] = Cscm[i, j]
                        C_tensor[voigt_notation[i][1], voigt_notation[i][0], voigt_notation[j][1], voigt_notation[j][0]] = Cscm[i, j]

        elif Cscm.shape == (3, 3):  # 2D material
            for i in range(3):
                for j in range(3):
                    C_tensor[i, j, i, j] = Cscm[i, j]
                    if i != j:
                        C_tensor[i, i, j, j] = Cscm[i, j]
                        C_tensor[j, j, i, i] = Cscm[i, j]
                        C_tensor[i, j, j, i] = Cscm[i, j]
                        C_tensor[j, i, i, j] = Cscm[i, j]

        else:
            raise ValueError("Elastic constants matrix must be either 3x3 or 6x6.")

        return C_tensor


    def invert_elastic_tensor(self):
        """
        Invert a 3x3x3x3 elastic tensor to get the compliance tensor.
        """
        C = self.C_tensor
        # Reshape the 3x3x3x3 tensor into a 9x9 matrix for inversion
        C_matrix = np.reshape(C, (9, 9))
        S_matrix = np.linalg.inv(C_matrix)

        # Reshape back into a 3x3x3x3 tensor
        S = np.reshape(S_matrix, (3, 3, 3, 3))
        return S
        
    def christoffel_matrix(self, n):
        """
        Compute the Christoffel matrix for a given direction n.
        n: Direction vector.
        """
        #if self.Cdim.shape == (6, 6):
        C = self.C_tensor *1E+9 # in N/m^2
        #else:
        #    C = self.Cdim  # Directly use the 3x3 matrix for 2D materials

        Gamma = np.zeros((3, 3))
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for l in range(3):
                        Gamma[i, j] += C[i, j, k, l] * n[k] * n[l]
        return Gamma

    def group_velocity(self, n):
        """
        Calculate the group velocity for a given direction n.
        n: Direction vector.
        """
        Gamma = self.christoffel_matrix(n)
        eigenvalues = np.linalg.eigvalsh(Gamma)
        sorted_eigenvalues = np.sort(eigenvalues)
        return np.sqrt(sorted_eigenvalues / self.rho_mass)/1000 # in km/s
        


    def christoffel_matrix_2D(self, n):
        """
        Compute the Christoffel matrix for a given direction n and 3x3 elastic constants C for 2D materials.
        C: Elastic constants matrix (3x3 for 2D materials).
        n: Direction vector (in-plane).
        """
        C = self.C_tensor
        Gamma = np.zeros((2, 2))
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    for l in range(2):
                        Gamma[i, j] += C[i, k, j, l] * n[k] * n[l]
        return Gamma

    def group_velocity_2D(self, n):
        Gamma = self.christoffel_matrix_2D(n)
        eigenvalues = np.linalg.eigvalsh(Gamma)
        sorted_eigenvalues = np.sort(eigenvalues)  # Sort the eigenvalues
        return np.sqrt(abs(sorted_eigenvalues) / self.rho_mass)/1000 # in km/s
        
                
    def compute_moduli_contour_2D(self, theta, phi):
        epsilon = 1e-8  # Small regularizing constant

        C11, C12, C16, _, C22, C26, _, _, C66 = self.Cdim.ravel()
        
        x = np.sin(theta) * np.cos(phi)
        y = np.sin(theta) * np.sin(phi)

        # Modified stiffness constants
        C11_prime = C11 * x**4 + C22 * y**4 + 2 * (C12 + 2 * C66) * x**2 * y**2 + 4 * (C16 * x**3 * y + C26 * x * y**3)
        C22_prime = C11 * y**4 + C22 * x**4 + 2 * (C12 + 2 * C66) * x**2 * y**2 + 4 * (C16 * x * y**3 + C26 * x**3 * y)
        C12_prime = (C11 + C22 - 4 * C66) * x**2 * y**2 + C12 * (x**4 + y**4) + 2 * (C16 * x * (x**2 - y**2) + C26 * y * (y**2 - x**2))
        C66_prime = (C11 + C22 - 2 * C12) * x**2 * y**2 + C66 * (x**4 + y**4) + 2 * (C16 - C26) * x * y * (x**2 - y**2)

#        zero_block = np.zeros_like(C11_prime)
#        C_prime = np.block([[C11_prime, C12_prime, zero_block],
#                            [C12_prime, C22_prime, zero_block],
#                            [zero_block, zero_block, C66_prime]])
                 
        # Calculate compliance values
        S11_prime = 1 / (C11_prime+ epsilon)
        S22_prime = 1 / (C22_prime + epsilon)
        S66_prime = 1 / (C66_prime + epsilon)

        # Compute S12_prime, guarding against division by zero
        denominator = C11_prime * C22_prime - C12_prime**2
        S12_prime = np.zeros_like(denominator)
        non_zero_mask = np.abs(denominator) >= epsilon
        S12_prime[non_zero_mask] = -C12_prime[non_zero_mask] / denominator[non_zero_mask]
        #S12_prime = np.where(np.abs(denominator) < epsilon, 0.0, -C12_prime / denominator)

        # Calculate moduli
        Ex = 1 / S11_prime
        Ey = 1 / S22_prime
        Gxy = 1 / S66_prime

        # Poisson's ratios
        nu_xy = -S12_prime / S11_prime
        nu_yx = -S12_prime / S22_prime

        # Stiffness constant
        K_prime = (C11_prime + C22_prime + 2 * C12_prime) / 3

        return Ex, Ey, Gxy, nu_xy, nu_yx, K_prime


    # Transform the stiffness matrix
    def transform_stiffness_3D(self, theta, phi):
        # Extract the 3x3 submatrix
        C_sub = self.Cdim[:3, :3]

        # Transformation matrix
        M = np.array([
            [np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)],
            [np.cos(theta) * np.cos(phi), np.cos(theta) * np.sin(phi), -np.sin(theta)],
            [-np.sin(phi), np.cos(phi), 0]
        ])
        return M @ C_sub @ M.T


    def safe_inverse(self,value, tol=1e-5):
        """Returns inverse if value is not too close to zero; otherwise returns infinity."""
        if np.isclose(value, 0, atol=tol):
            return np.inf
        else:
            return 1 / value



    # Compute directional moduli for 3D material
    def compute_moduli_directional_3D(self,Sij):
        #Sij = self.Cs

        # Young's moduli
        E1 = 1 / Sij[0, 0]
        E2 = 1 / Sij[1, 1]
        E3 = 1 / Sij[2, 2]

        # Shear moduli
        G23 = 1 / Sij[1, 2]
        G13 = 1 / Sij[0, 2]
        G12 = 1 / Sij[0, 1]

        # Poisson's ratios
        nu12 = -Sij[0, 1] * Sij[0, 0]
        nu21 = -Sij[1, 0] * Sij[1, 1]
        nu13 = -Sij[0, 2] * Sij[0, 0]
        nu31 = -Sij[2, 0] * Sij[2, 2]
        nu23 = -Sij[1, 2] * Sij[1, 1]
        nu32 = -Sij[2, 1] * Sij[2, 2]

        # Bulk modulus
        K = 1 / (Sij[0, 0] + 2 * Sij[0, 1])

        return E1, E2, E3, G23, G13, G12, nu12, nu21, nu13, nu31, nu23, nu32, K



    def transform_stiffness_2D(self, theta):
        # 2D transformation matrix for anisotropic plane strain
        C_sub = self.Cdim[:3, :3]
        M = np.array([
            [np.cos(theta), np.sin(theta), 0],
            [-np.sin(theta), np.cos(theta), 0],
            [0, 0, 1]
        ])
        return M @ C_sub @ M.T

    def compute_moduli_directional_2D(self, Sij):
        threshold = 1e-3
        close_to_zero_mask = np.abs(Sij) < threshold
        Sij[close_to_zero_mask] = 0
        # Young's moduli
        E1 = 1 / Sij[0, 0]
        E2 = 1 / Sij[1, 1]

        # Shear moduli
        G12 = 1 / self.safe_inverse(Sij[0, 1])
        G23 = 1 / self.safe_inverse(Sij[1, 2])
        G13 = 1 / self.safe_inverse(Sij[0, 2])

        # Poisson's ratios
        nu12 = -Sij[0, 1] * Sij[0, 0]
        nu21 = -Sij[1, 0] * Sij[1, 1]
        nu13 = -Sij[0, 2] * Sij[0, 0]
        nu31 = -Sij[2, 0] * Sij[2, 2]
        nu23 = -Sij[1, 2] * Sij[1, 1]
        nu32 = -Sij[2, 1] * Sij[2, 2]

        return E1, E2, G12, G23, G13, nu12, nu21, nu13, nu31, nu23, nu32





    def plot_orientation_dependent_2D(self, npoints=360, fname='EVGK_polar_2D', dpi=80):
        theta = np.linspace(0, 2*np.pi, npoints)

        E = self.compute_E_2D_polar(theta)
        V = self.compute_V_2D_polar(theta)
        G = self.compute_G_2D_polar(theta)
        K = self.compute_K_2D_polar(theta)

        # Save the data
        data = np.vstack((theta, E, V, G, K)).T
        np.savetxt(f"{fname}.dat", data, fmt='%10.4f %10.4f %10.4f %10.4f %10.4f', header='Theta E V G K')


        # Set global style parameters
        plt.rcParams.update({'font.size': 16, 'font.family': 'sans-serif'})

        sns.set_style("whitegrid")

        try:
            fig, axs = plt.subplots(2, 2, figsize=(12, 12), subplot_kw=dict(projection='polar'))

            # Define line properties
            line_properties = {'lw': 1, 'ls': '--', 'alpha': 0.6}

            # Young's Modulus
            axs[0, 0].plot(theta, E, color="tab:orange", marker='o', label="$E$", **line_properties)
            axs[0, 0].set_title("Young's Modulus", fontsize=16)

            # Poisson's Ratio
            axs[0, 1].plot(theta, V, color="tab:green", marker='o', label="$\nu$", **line_properties)
            axs[0, 1].set_title("Poisson's Ratio", fontsize=16)

            # Shear Modulus
            axs[1, 0].plot(theta, G, color="tab:blue", marker='o', label="$G$", **line_properties)
            axs[1, 0].set_title("Shear Modulus", fontsize=16)

            # Stiffness Constant
            axs[1, 1].plot(theta, K, color="tab:red", marker='o', label="$K$", **line_properties)
            axs[1, 1].set_title("Stiffness Constant", fontsize=16)

            plt.tight_layout()
            plt.savefig(f"{fname}.png", format='png', dpi=dpi)
            plt.close(fig) 


            if self.plotly:
                fig_plotly = make_subplots(rows=2, cols=2, subplot_titles=("Young's Modulus", "Poisson's Ratio", "Shear Modulus", "Stiffness Constant"), specs=[[{'type': 'polar'}, {'type': 'polar'}], [{'type': 'polar'}, {'type': 'polar'}]])

                fig_plotly.update_layout(
                    title_text="Material visualization with ElasTool, developed by C.E. Ekuma. If you use the software in your research, please cite Sci. Rep. 12, 3776 (2022) ",
                    font=dict(size=12, family='sans-serif'),
                    annotations=[
                        dict(xref='paper', yref='paper' ) ]  )

                # Add traces
                fig_plotly.add_trace(go.Scatterpolar(r=E, theta=theta*180/np.pi, name='E', marker=dict(symbol='circle')), 1, 1)
                fig_plotly.add_trace(go.Scatterpolar(r=V, theta=theta*180/np.pi, name='V', marker=dict(symbol='circle')), 1, 2)
                fig_plotly.add_trace(go.Scatterpolar(r=G, theta=theta*180/np.pi, name='G', marker=dict(symbol='circle')), 2, 1)
                fig_plotly.add_trace(go.Scatterpolar(r=K, theta=theta*180/np.pi, name='K', marker=dict(symbol='circle')), 2, 2)


                # Update layout
                fig_plotly.write_image(f"{fname}plotly.png")  # Save the image
                fig_plotly.show()


        except Exception as e:
            print(f"Error while plotting: {e}")



    def plot_orientation_dependent_3D(self, npoints=360, fname='EVGK_polar_3D',dpi = 80):
        theta = np.linspace(0, np.pi, npoints)
        phi = np.linspace(0, 2 * np.pi, npoints)

        theta, phi = np.meshgrid(theta, phi)

        u = np.sin(theta) * np.cos(phi)
        v = np.sin(theta) * np.sin(phi)
        w = np.cos(theta)
        l = np.array([
            u**2,
            v**2,
            w**2,
            2*u*v,
            2*v*w,
            2*u*w
        ]).transpose(1, 2, 0)

        rho_e = np.empty_like(theta)
        rho_v = np.empty_like(theta)
        rho_g = np.empty_like(theta)
        rho_k = np.empty_like(theta)

 
        for i in range(npoints):
            for j in range(npoints):
                l_point = l[i, j, :]  # Extract the 6-component vector at this point
                rho_e[i, j] = self.compute_E_general(l_point)
                rho_v[i, j] = self.compute_V_general(l_point)
                rho_g[i, j] = self.compute_G_general(l_point)
                rho_k[i, j] = self.compute_K_general(l_point)



        x_e, y_e, z_e = rho_e * np.sin(theta) * np.cos(phi), rho_e * np.sin(theta) * np.sin(phi), rho_e * np.cos(theta)
        x_v, y_v, z_v = rho_v * np.sin(theta) * np.cos(phi), rho_v * np.sin(theta) * np.sin(phi), rho_v * np.cos(theta)
        x_g, y_g, z_g = rho_g * np.sin(theta) * np.cos(phi), rho_g * np.sin(theta) * np.sin(phi), rho_g * np.cos(theta)
        x_k, y_k, z_k = rho_k * np.sin(theta) * np.cos(phi), rho_k * np.sin(theta) * np.sin(phi), rho_k * np.cos(theta)

        sns.set_style("whitegrid")
        plt.rcParams.update({'font.size': 16, 'font.family': 'sans-serif'})
        try:    
            fig = plt.figure(figsize=(12, 10))

            ax_e = fig.add_subplot(2, 2, 1, projection='3d')
            ax_e.plot_surface(x_e, y_e, z_e, cmap='viridis')
            ax_e.set_title("Young's Modulus",fontsize=16)

            ax_e.set_xlabel('X')
            ax_e.set_ylabel('Y')
            ax_e.set_zlabel('Z')

            ax_v = fig.add_subplot(2, 2, 2, projection='3d')
            ax_v.plot_surface(x_v, y_v, z_v, cmap='viridis')
            ax_v.set_title("Poisson's Ratio",fontsize=16)

            ax_v.set_xlabel('X')
            ax_v.set_ylabel('Y')
            ax_v.set_zlabel('Z')

            ax_g = fig.add_subplot(2, 2, 3, projection='3d')
            ax_g.plot_surface(x_g, y_g, z_g, cmap='viridis')
            ax_g.set_title("Shear Modulus",fontsize=16)

            ax_g.set_xlabel('X')
            ax_g.set_ylabel('Y')
            ax_g.set_zlabel('Z')

            ax_k = fig.add_subplot(2, 2, 4, projection='3d')
            ax_k.plot_surface(x_k, y_k, z_k, cmap='viridis')
            ax_k.set_title("Bulk Modulus",fontsize=16)

            ax_k.set_xlabel('X')
            ax_k.set_ylabel('Y')
            ax_k.set_zlabel('Z')

            plt.tight_layout()
            plt.savefig(f"{fname}.png", format='png', dpi=dpi)
            plt.close(fig) 
            #plt.show()

            if self.plotly:


                fig_plotly = make_subplots(
                    rows=2, cols=2,
                    specs=[[{'type': 'surface'}, {'type': 'surface'}],
                    [{'type': 'surface'}, {'type': 'surface'}]],
                    subplot_titles=('Young\'s Modulus', 'Poisson\'s Ratio', 'Shear Modulus', 'Bulk Modulus')
                   )

                # Young's Modulus plot
                fig_plotly.add_trace(
                   go.Surface(z=z_e, x=x_e, y=y_e, colorscale='Viridis', showscale=False),
                   row=1, col=1
                )

                # Poisson's Ratio plot
                fig_plotly.add_trace(
                    go.Surface(z=z_v, x=x_v, y=y_v, colorscale='Viridis', showscale=False),
                    row=1, col=2
                )

                # Shear Modulus plot
                fig_plotly.add_trace(
                    go.Surface(z=z_g, x=x_g, y=y_g, colorscale='Viridis', showscale=False),
                    row=2, col=1
                )

                # Bulk Modulus plot
                fig_plotly.add_trace(
                    go.Surface(z=z_k, x=x_k, y=y_k, colorscale='Viridis', showscale=False),
                    row=2, col=2
                )

                # Update layouts and axes titles
                fig_plotly.update_layout(scene=dict(xaxis_title='X', yaxis_title='Y', zaxis_title='Z'),
                           scene2=dict(xaxis_title='X', yaxis_title='Y', zaxis_title='Z'),
                           scene3=dict(xaxis_title='X', yaxis_title='Y', zaxis_title='Z'),
                           scene4=dict(xaxis_title='X', yaxis_title='Y', zaxis_title='Z'),
                           margin=dict(l=0, r=0, b=40, t=40),
                           title_text="Material visualization with ElasTool, developed by C.E. Ekuma. If you use the software in your research, please cite Comput. Phys. Commun. 270, 108180, (2022) ")

                fig_plotly.write_image(f"{fname}plotly.png")
                fig_plotly.show()

        except Exception as e:
            print(f"Error while plotting: {e}")


#Work on this
    def plot_contour_2D(self,npoints=200, fname='EVGK_contour_2D',dpi = 200):

        norm = PowerNorm(0.75)  # gamma value, which you can adjust


        theta = np.linspace(0, np.pi, npoints)
        phi = np.linspace(0, 2*np.pi, npoints)
        Theta, Phi = np.meshgrid(theta, phi)

        Ex, Ey, Gxy, nu_xy, nu_yx, K = self.compute_moduli_contour_2D(Theta, Phi)


        try:
            fig, axes = plt.subplots(3, 2, figsize=(10, 10)) #,subplot_kw=dict(projection='polar'))

            for ax in axes.ravel():
                ax.set_facecolor('white')
           # for i in range(3):
           #     for j in range(2):
           #         axes[i, j].set_theta_direction(-1)  # Clockwise
           #         axes[i, j].set_theta_offset(np.pi)  # Starting from the top

            # Young's Modulus Ex
            brighten_color = self.adjust_colormap_brightness('rainbow', 0.05)
            colormap_obj = plt.get_cmap(brighten_color)
            levels = 50
            contour1 = axes[0, 0].contourf(Theta, Phi, Ex, cmap=colormap_obj, levels=levels,norm=norm)
            fig.colorbar(contour1, ax=axes[0, 0], label="Young's Modulus Ex")
            axes[0, 0].set_title("Young's Modulus Ex")
            axes[0, 0].set_xlabel('θ (radians)')
            axes[0, 0].set_ylabel('φ (radians)')

            # Young's Modulus Ey
            contour2 = axes[0, 1].contourf(Theta, Phi, Ey, cmap=colormap_obj, levels=levels)
            fig.colorbar(contour2, ax=axes[0, 1], label="Young's Modulus Ey")
            axes[0, 1].set_title("Young's Modulus Ey")
            axes[0, 1].set_xlabel('θ (radians)')
            axes[0, 1].set_ylabel('φ (radians)')

            # Shear Modulus Gxy
            contour3 = axes[1, 0].contourf(Theta, Phi, Gxy, cmap=colormap_obj, levels=levels)
            fig.colorbar(contour3, ax=axes[1, 0], label="Shear Modulus Gxy")
            axes[1, 0].set_title("Shear Modulus Gxy")
            axes[1, 0].set_xlabel('θ (radians)')
            axes[1, 0].set_ylabel('φ (radians)')

            #levels_nu = np.linspace(nu_xy.min(), nu_xy.max(), levels)

            # Poisson's Ratio νxy
            contour4 = axes[1, 1].contourf(Theta, Phi, nu_xy, cmap=colormap_obj, levels=levels)
            fig.colorbar(contour4, ax=axes[1, 1], label="Poisson's Ratio νxy")
            axes[1, 1].set_title("Poisson's Ratio νxy")
            axes[1, 1].set_xlabel('θ (radians)')
            axes[1, 1].set_ylabel('φ (radians)')

            # Bulk Modulus K
            contour5 = axes[2, 0].contourf(Theta, Phi, K, cmap=colormap_obj, levels=levels)
            fig.colorbar(contour5, ax=axes[2, 0], label="Stiffness Constant K")
            axes[2, 0].set_title("Stiffness Constant K")
            axes[2, 0].set_xlabel('θ (radians)')
            axes[2, 0].set_ylabel('φ (radians)')
 
            # Poisson's Ratio νyx
            contour6 = axes[2, 1].contourf(Theta, Phi, nu_yx, cmap=colormap_obj, levels=levels)
            fig.colorbar(contour6, ax=axes[2, 1], label="Poisson's Ratio νyx")
            axes[2, 1].set_title("Poisson's Ratio νyx")
            axes[2, 1].set_xlabel('θ (radians)')
            axes[2, 1].set_ylabel('φ (radians)')

            plt.tight_layout(pad=2.0)
            #plt.show()

            plt.savefig(f"{fname}.png", format='png', dpi=dpi,facecolor=fig.get_facecolor())
            plt.close(fig) 

        except Exception as e:
            print(f"Error while plotting 2D contour: {e}")



    def plot_contour_polar_3D(self,npoints=200, fname='EVGK_polar_contour_3D',dpi = 200):

        theta = np.linspace(0, np.pi, npoints)
        phi = np.linspace(0, 2*np.pi, npoints)
        Theta, Phi = np.meshgrid(theta, phi)

        # Calculate the transformed stiffness matrix and moduli
        C_prime = np.array([self.transform_stiffness_3D(t, p) for t, p in zip(np.ravel(Theta), np.ravel(Phi))])
        moduli = np.array([self.compute_moduli_directional_3D(np.linalg.inv(cp)) for cp in C_prime])



        try:
            # Visualization
            fig, axes = plt.subplots(5, 3, figsize=(15, 20), subplot_kw={'projection': 'polar'})

      
            titles = [
                ["Young's Modulus E1", "Young's Modulus E2", "Young's Modulus E3"],
                ["Shear Modulus G23", "Shear Modulus G13", "Shear Modulus G12"],
                ["Poisson's Ratio ν12", "Poisson's Ratio ν21", "Poisson's Ratio ν13"],
                ["Poisson's Ratio ν31", "Poisson's Ratio ν23", "Poisson's Ratio ν32"],
                ["Bulk Modulus K", "", ""]
            ]
            E1, E2, E3, G23, G13, G12, nu12, nu21, nu13, nu31, nu23, nu32, K = [x.reshape(Theta.shape) for x in moduli.T]
            data = [E1, E2, E3, G23, G13, G12, nu12, nu21, nu13, nu31, nu23, nu32, K]
            levels = 10
            for i in range(5):
                for j in range(3):
                    if titles[i][j]:  # Only plot if there's a title
                        c = axes[i][j].contourf(Phi, Theta, data[i*3 + j], cmap='coolwarm', levels=levels)
                        axes[i][j].set_title(titles[i][j])
                        fig.colorbar(c, ax=axes[i][j])
                        contours = axes[i][j].contour(Phi, Theta, data[i*3 + j], levels=2, colors='black')
                        plt.clabel(contours, inline=True, fontsize=8)
                    else:
                        axes[i][j].set_visible(False)  # Hide unused subplots


            plt.tight_layout()
            #plt.show()

            plt.savefig(f"{fname}.png", format='png', dpi=dpi)
            plt.close(fig) 
            
        except Exception as e:
            print(f"Error while plotting 3D contour: {e}")



    def heatmap_3D(self, resolution=200, file_name='groupvelocity_heatmap_3D'):
        """
        Generate a 3D heat map of group velocities for each mode and save as a high resolution PNG file.
        """
        # Create a grid of directions
        import plotly.express as px
        theta = np.linspace(0, np.pi, resolution)
        phi = np.linspace(0, 2 * np.pi, resolution)
        theta, phi = np.meshgrid(theta, phi)
        n_x = np.sin(theta) * np.cos(phi)
        n_y = np.sin(theta) * np.sin(phi)
        n_z = np.cos(theta)

        # Calculate group velocities
        v_g = np.array([self.group_velocity(np.array([nx, ny, nz])) for nx, ny, nz in zip(n_x.flatten(), n_y.flatten(), n_z.flatten())])
        v_g = v_g.reshape(theta.shape + (3,))

                
        # Plot the heat maps for each mode
        try:
            fig, axes = plt.subplots(1, 3, figsize=(15, 5))
            titles = ['Transverse Mode 1', 'Transverse Mode 2', 'Longitudinal Mode']
            for i in range(3):
                im = axes[i].imshow(v_g[:, :, i], extent=[0, 2 * np.pi, 0, np.pi], origin='lower', aspect='auto')
                axes[i].set_title(titles[i])
                axes[i].set_xlabel('φ (radians)')
                axes[i].set_ylabel('θ (radians)')
                axes[i].grid(False)
                cbar = fig.colorbar(im, ax=axes[i])
                cbar.set_label('Sound Velocity (km/s)')
            plt.tight_layout()

            plt.savefig(f"{file_name}.png", format='png', dpi=300)
            #plt.show()
            plt.close(fig)
            
            if self.plotly:
              axis_font_size = 20
              tick_font_size = 20
              titles = ['Transverse Mode 1', 'Transverse Mode 2', 'Longitudinal Mode']
              fig = make_subplots(rows=1, cols=3, subplot_titles=titles,horizontal_spacing=0.1)
              colorbar_positions = [0.262, 0.63, 0.998]
              #colorbar_positions = [0.284, 0.641, 0.998]
              for i in range(3):
                 
                  fig.add_trace(go.Heatmap(
                      z=v_g[:, :, i],
                      x=np.linspace(0, 2*np.pi, resolution),
                      y=np.linspace(0, np.pi, resolution),
                      colorscale='RdBu_r',
                      colorbar=dict(
                          title='V (km/s)',
                          x= colorbar_positions[i]
                      )
                  ), row=1, col=i+1)

                  # Update subplot axis titles
                  fig.update_xaxes(title_text=r'$\Phi$ (radians)', row=1, col=i+1, title_font=dict(size=axis_font_size),tickfont=dict(size=tick_font_size))
                  fig.update_yaxes(title_text=r'$\Theta$ (radians)', row=1, col=i+1, title_font=dict(size=axis_font_size),tickfont=dict(size=tick_font_size))

              # Update layout for the entire figure
              fig.update_layout(
                  title_text="Sound velocity modes with ElasTool, developed by C.E. Ekuma. If you use the software in your research, please cite Comput. Phys. Commun. 270, 108180, (2022)",
                  height=500, width=1650
              )

              fig.write_image(f"{file_name}plotly.png")
              fig.show()

                # Save the plot as a static image
                #fig.write_image(file_name)
        except Exception as e:
            print(f"Error while plotting 3D heatmap: {e}")        
        


    def heatmap_2D(self, resolution=200, file_name='groupvelocity_contour_2D'):
        """
        Generate a 2D heat map of group velocities for each mode and save as a high resolution PNG file.
        For 2D materials, only in-plane directions are considered.
        """
        # Create a grid of directions in the XY plane
        theta = np.linspace(0, 2 * np.pi, resolution)
        n_x = np.cos(theta)
        n_y = np.sin(theta)
        # Calculate group velocities in the XY plane
        v_g = np.array([self.group_velocity_2D(np.array([nx, ny,0])) for nx, ny in zip(n_x, n_y)])
        v_g = v_g.reshape((resolution, 2))
        # Plot the heat maps for each mode
        
        try:
            fig, axes = plt.subplots(1, 2, figsize=(12, 6), subplot_kw={'projection': 'polar'})
            titles = ['Transverse', 'Longitudinal']
            for i in range(2):
                axes[i].plot(theta, v_g[:, i], color='grey', alpha=0.5)
                scatter = axes[i].scatter(theta, v_g[:, i], c=v_g[:, i], cmap='viridis', s=10)
                fig.colorbar(scatter, ax=axes[i], label='V (km/s)', pad=0.1)

                axes[i].set_title(titles[i])
                axes[i].set_theta_zero_location('N')
                axes[i].set_theta_direction(-1)
                axes[i].set_rlabel_position(90)
                #axes[i].set_xlabel(r'$\Theta$ (radians)')
                #axes[i].set_ylabel('Sound Velocity (km/s)',labelpad=30)
                #axes[i].grid(False)

            plt.tight_layout()
            plt.savefig(f"{file_name}.png", format='png', dpi=300)
            plt.close(fig)
            
            if self.plotly:
              fig = make_subplots(rows=1, cols=2, specs=[[{'type': 'polar'}, {'type': 'polar'}]], subplot_titles=['Transverse', 'Longitudinal'])

              # Add plots to each subplot
              for i in range(2):
                  fig.add_trace(go.Scatterpolar(
                      r=v_g[:, i], 
                      theta=np.rad2deg(theta), 
                      mode='lines+markers',
                      marker=dict(color=v_g[:, i], size=5, colorscale='Viridis', colorbar=dict(title='V (km/s)', x=(0.5 if i == 0 else 1))),
                      #name=['Transverse', 'Longitudinal'][i]
                  ), row=1, col=i+1)
              tick_font_size = 18
              # Update the layout for polar plots
              fig.update_layout(
                  polar=dict(
                      angularaxis=dict(
                          rotation=90,
                          direction='clockwise',
                          thetaunit='radians',
                          tickfont=dict(size=tick_font_size)
                      )
                  ),
                  polar2=dict(
                      angularaxis=dict(
                          rotation=90,
                          direction='clockwise',
                          thetaunit='radians',
                          tickfont=dict(size=tick_font_size)
                      )
                  ),
                  title_text="Sound velocity modes with ElasTool, developed by C.E. Ekuma. If you use the software in your research, please cite Sci. Rep. 12, 3776 (2022)",
                  height=800,  # Adjust as necessary
                  width=1600   # Adjust as necessary
              )
              fig.update_layout(showlegend=False)
              # Show the plot
              fig.write_image(f"{file_name}plotly.png")
              fig.show() 
              
            
        except Exception as e:
            print(f"Error while plotting 2D heatmap: {e}")
            

    def heatmapN_2D(self, resolution=200, file_name='groupvelocity_heatmap_2D.png'):
        """
        Generate a 2D heat map of group velocities for each mode and save as a high resolution PNG file.
        For 2D materials, only in-plane directions are considered.
        """
        # Create a grid of directions in the XY plane
        theta = np.linspace(0, 2 * np.pi, resolution)
        phi = np.linspace(0, 2*np.pi, resolution)
        theta, phi = np.meshgrid(theta, phi)
        n_x = np.sin(phi) * np.cos(theta)
        n_y = np.sin(phi) * np.sin(theta)


        # Calculate group velocities in the XY plane
        v_g = np.array([[group_velocity_2D(np.array([nx, ny, 0])) for nx, ny in zip(n_x[i], n_y[i])] for i in range(len(n_x))])
        

        try:
            fig, axes = plt.subplots(1, 2, figsize=(12, 6), subplot_kw={'projection': 'polar'})
            titles = ['Transverse Mode', 'Longitudinal Mode']

            for i in range(2):
                # Filled contour plot
                contour = axes[i].contourf(theta, phi, v_g[:, :, i], 20, cmap='viridis')
                fig.colorbar(contour, ax=axes[i], label='Sound Velocity (km/s)', pad=0.1)
                
                # Overlay regular line plot
                axes[i].plot(theta[0, :], v_g[:, 0, i], color='red', linewidth=2)

                # Set plot attributes
                axes[i].set_title(titles[i])
                axes[i].set_theta_zero_location('N')
                axes[i].set_theta_direction(-1)
                axes[i].set_rlabel_position(90)
                axes[i].set_xticks(np.linspace(0, 2*np.pi, 4))
                axes[i].set_xticklabels(['0', 'π/2', 'π', ''])

            plt.tight_layout()
            plt.savefig(file_name, dpi=300)  # Save the figure as a high-resolution PNG file
            plt.show()
        except Exception as e:
            print(f"Error while plotting 2D2 heatmap: {e}")
            

    def plot_contour_polar_2D(self, npoints=200, fname='EVGK_polar_contour_2D', dpi=200):
        theta_values = np.linspace(0, 2 * np.pi, npoints)
        radius = np.array([0, 1])

        # Calculate the transformed stiffness matrix and moduli
        C_prime = np.array([self.transform_stiffness_2D(t) for t in np.ravel(theta_values)])
        
        moduli = np.array([self.compute_moduli_directional_2D(np.linalg.pinv(cp)) for cp in C_prime])




        # Visualization
        titles = [
            "Young's Modulus E1", "Young's Modulus E2", "Shear Modulus G12",
            "Shear Modulus G13","Shear Modulus G23", "Poisson's Ratio ν12",
            "Poisson's Ratio ν21", "Poisson's Ratio ν13","Poisson's Ratio ν23"
        ]
        fig, axes = plt.subplots(2, 4, figsize=(15, 8), subplot_kw={'projection': 'polar'})

        theta, radius = np.meshgrid(theta_values, radius)

        levels = 50
        for i, ax in enumerate(axes.ravel()):
            if i < len(titles):
                z_data = np.tile(moduli[:, i], (2,1))

                c = ax.contourf(theta, radius, z_data, cmap='rainbow', levels=levels)
                ax.set_title(titles[i])

                try:
                    fig.colorbar(c, ax=ax)
                except ValueError as e:
                    print(f"Failed to create colorbar for plot {i}. Error: {e}")


                #fig.colorbar(c, ax=ax)
                contours = ax.contour(theta, radius, z_data, levels=2, colors='black')
                plt.clabel(contours, inline=True, fontsize=8)
            else:
                ax.set_visible(False)  # Hide unused subplots

        plt.tight_layout()
        plt.savefig(f"{fname}.png", format='png', dpi=dpi)
        plt.close(fig)
        #plt.show()


#class ElasticParametersHeatMap:
#    def __init__(self, elastic_tensor,plot=False,plotly=False):
#        self.Cdim = elastic_tensor
#        self.plotly = plotly
#        self.plot = plot

    def compute_youngs_modulus(self, direction):
        S = np.linalg.inv(self.Cdim) 
        
        # Normalize the direction vector
        n = np.array(direction)
        n /= np.linalg.norm(n)

        # Convert direction vector to Voigt notation (strain/stress vector)
        voigt_direction = np.array([
            n[0]**2, n[1]**2, n[2]**2,
            2*n[1]*n[2], 2*n[0]*n[2], 2*n[0]*n[1]
        ])

        # Calculate Young's modulus in Voigt notation
        E_inv = np.dot(voigt_direction, np.dot(S, voigt_direction))
        E = 1 / E_inv if E_inv != 0 else float('inf')
        return E


    def compute_poissons_ratioold(self, direction_stress, direction_strain):
        S = np.linalg.inv(self.Cdim)
        
        # Normalize direction vectors
        n_stress = np.array(direction_stress) / np.linalg.norm(direction_stress)
        n_strain = np.array(direction_strain) / np.linalg.norm(direction_strain)

        # Convert direction vectors to Voigt notation
        voigt_stress = np.array([
            n_stress[0]**2, n_stress[1]**2, n_stress[2]**2,
            2*n_stress[1]*n_stress[2], 2*n_stress[0]*n_stress[2], 2*n_stress[0]*n_stress[1]
        ])
        voigt_strain = np.array([
            n_strain[0]**2, n_strain[1]**2, n_strain[2]**2,
            2*n_strain[1]*n_strain[2], 2*n_strain[0]*n_strain[2], 2*n_strain[0]*n_strain[1]
        ])

        # Calculate Poisson's ratio
        nu = -np.dot(voigt_strain, np.dot(S, voigt_stress)) / np.dot(voigt_stress, np.dot(S, voigt_stress))
        return nu


    def compute_poissons_ratio(self, direction_stress):
        S = np.linalg.inv(self.Cdim)
        
        # Fixed direction for strain along the x-axis
        direction_strain = [1, 0, 0]

        # Normalize direction vector for stress
        n_stress = np.array(direction_stress) / np.linalg.norm(direction_stress)

        # Convert direction vectors to Voigt notation
        voigt_stress = np.array([
            n_stress[0]**2, n_stress[1]**2, n_stress[2]**2,
            2*n_stress[1]*n_stress[2], 2*n_stress[0]*n_stress[2], 2*n_stress[0]*n_stress[1]
        ])
        voigt_strain = np.array([1, 0, 0, 0, 0, 0])  # Strain along x-axis in Voigt notation

        # Calculate Poisson's ratio
        nu = np.dot(voigt_strain, np.dot(S, voigt_stress)) / np.dot(voigt_stress, np.dot(S, voigt_stress))
        return nu


    def compute_shear_modulus(self, direction):
        S = np.linalg.inv(self.Cdim)
        
        # Normalize the direction vector
        n = np.array(direction)
        n /= np.linalg.norm(n)

        # Convert direction vector to Voigt notation
        voigt_direction = np.array([
            2*n[1]*n[2], 2*n[0]*n[2], 2*n[0]*n[1], # Shear components
            n[0]**2, n[1]**2, n[2]**2              # Normal components
        ])

        # Calculate Shear modulus
        G_inv = np.dot(voigt_direction, np.dot(S, voigt_direction))
        G = 1 / G_inv if G_inv != 0 else float('inf')
        return G


    def plot_moduli_heatmaps(self, n_points=100,fname='EVK_heatmap_3D',dpi = 100):  #fname='EVGK_polar_3D',dpi = 80)
        theta = np.linspace(0, np.pi, n_points)
        phi = np.linspace(0, 2 * np.pi, n_points)
        theta_grid, phi_grid = np.meshgrid(theta, phi)

        youngs_modulus_map = np.zeros_like(theta_grid)
        #poisson_ratio = self.compute_poissons_ratio()
        poisson_ratio_map = np.zeros_like(theta_grid)  
        shear_modulus_map = np.zeros_like(theta_grid)

        for i in range(n_points):
            for j in range(n_points):
                direction = [np.sin(theta_grid[i, j]) * np.cos(phi_grid[i, j]),
                             np.sin(theta_grid[i, j]) * np.sin(phi_grid[i, j]),
                             np.cos(theta_grid[i, j])]
                youngs_modulus_map[i, j] = self.compute_youngs_modulus(direction)
                shear_modulus_map[i, j] = self.compute_shear_modulus(direction)
                poisson_ratio_map[i, j] = self.compute_poissons_ratio(direction)

        plt.rcParams.update({'font.size': 16, 'font.family': 'sans-serif'})
        fig = plt.figure(figsize=(18, 6))

        plt.subplot(1, 3, 1)
        plt.pcolormesh(phi_grid, theta_grid, youngs_modulus_map, shading='gouraud') # ['gouraud', 'nearest', 'flat', 'auto']. Setting shading='auto'.
        plt.colorbar(label='Young\'s Modulus (GPa)')
        #plt.title('Young\'s Modulus')
        plt.xlabel('Azimuthal Angle φ (radians)')
        plt.ylabel('Polar Angle θ (radians)')

        plt.subplot(1, 3, 2)
        plt.pcolormesh(phi_grid, theta_grid, poisson_ratio_map, shading='gouraud')
        plt.colorbar(label='Poisson\'s Ratio along [100]')
        #plt.title('Poisson\'s Ratio')
        plt.xlabel('Azimuthal Angle φ (radians)')
        plt.ylabel('Polar Angle θ (radians)')

        plt.subplot(1, 3, 3)
        plt.pcolormesh(phi_grid, theta_grid, shear_modulus_map, shading='gouraud')
        plt.colorbar(label='Shear Modulus (GPa)')
        #plt.title('Shear Modulus')
        plt.xlabel('Azimuthal Angle φ (radians)')
        plt.ylabel('Polar Angle θ (radians)')

        plt.tight_layout()
        plt.savefig(f"{fname}.png", format='png', dpi=dpi)
        plt.close(fig) 
        #plt.show()
        if self.plotly:
            import plotly.graph_objects as go
            from plotly.subplots import make_subplots


            fig = make_subplots(rows=1, cols=3, subplot_titles=['Young\'s Modulus', 'Poisson\'s Ratio', 'Shear Modulus'],
                                specs=[[{'type': 'xy'}, {'type': 'xy'}, {'type': 'xy'}]],
                                horizontal_spacing=0.1)
            colorbar_positions = [0.262, 0.63, 0.998]
            # Define plot titles and labels
            titles = ['Young\'s Modulus', 'Poisson\'s Ratio', 'Shear Modulus']
            xlabels = ['φ (radians)'] * 3
            ylabels = ['θ (radians)'] * 3
            colorbars = ['E (GPa)', 'ν', 'G (GPa)']
            axis_font_size = 20
            tick_font_size = 20 

            # Add heatmaps to each subplot
            for i, data in enumerate([youngs_modulus_map, poisson_ratio_map, shear_modulus_map]):
                fig.add_trace(go.Heatmap(
                    x=phi_grid.flatten(), 
                    y=theta_grid.flatten(), 
                    z=data.flatten(), 
                    #colorbar=dict(title=colorbars[i], len=0.3, y=0.5),
                    colorbar=dict( title=colorbars[i], x= colorbar_positions[i] ),
                    colorscale='rainbow'
                ), row=1, col=i+1)

                # Update axes labels for each subplot
                fig.update_xaxes(title_text=xlabels[i], row=1, col=i+1, title_font=dict(size=axis_font_size),tickfont=dict(size=tick_font_size))
                fig.update_yaxes(title_text=ylabels[i], row=1, col=i+1, title_font=dict(size=axis_font_size),tickfont=dict(size=tick_font_size))

            # Update layout for the entire figure
            fig.update_layout(
                title_text="Contour plot visualization of elastic parameters with ElasTool, developed by C.E. Ekuma. If you use the software in your research, please cite Comput. Phys. Commun. 270, 108180, (2022)",
                height=500, width=1650
            )

            fig.write_image(f"{fname}plotly.png")
            fig.show()

            

    def linear_compressibility(self,direction):
        """
        Calculate the linear compressibility for an anisotropic material
        S is the 6x6/3x3 compliance matrix, direction is a 3D/2D unit vector.
        """
        #direction = np.array([1, 0, 0])  # Replace with your direction vector
        dirr = direction / np.linalg.norm(direction)
        #S_tensor = self.Cs_tensor  #convert_6x6/3x3_to_3x3x3x3(S)
        beta = np.einsum('ijkl,i,j,k,l', self.Cs_tensor, dirr, dirr, dirr, dirr)
        return beta



    def linear_compressibilityold(self):
        # Extract the compliance tensor components
        compliance_tensor = self.Cs_tensor
        C11 = compliance_tensor[0, 0, 0, 0]
        C22 = compliance_tensor[1, 1, 1, 1]
        C33 = compliance_tensor[2, 2, 2, 2]
        C12 = compliance_tensor[0, 0, 1, 1]
        C13 = compliance_tensor[0, 0, 2, 2]
        C23 = compliance_tensor[1, 1, 2, 2]

        # Calculate the trace of the compliance tensor
        trace = C11 + C22 + C33 + 2 * (C12 + C13 + C23)

        # Calculate the linear compressibility
        linear_compressibility = -trace

        return linear_compressibility


# Other ways. Both are valid, but not adapted to suit this class
    def calc_linear_compressibility_3D(self, S, theta, phi):
        """
        Calculate linear compressibility for 3D materials.

        Parameters:
        S (numpy.ndarray): Compliance matrix for 3D material.
        theta (float): Angle theta in radians.
        phi (float): Angle phi in radians.

        Returns:
        float: Linear compressibility.
        """
        # Calculate components of the wave vector in 3D
        k1 = np.sin(theta) * np.cos(phi)
        k2 = np.sin(theta) * np.sin(phi)
        k3 = np.cos(theta)

        # Calculate k products for 3D
        k11 = k1 * k1
        k12 = k1 * k2
        k13 = k1 * k3
        k22 = k2 * k2
        k23 = k2 * k3
        k33 = k3 * k3

        # Calculate linear compressibility for 3D
        aa = ((S[0, 0] + S[0, 1] + S[0, 2]) * k11 +
              (S[0, 5] + S[1, 5] + S[2, 5]) * k12 +
              (S[0, 4] + S[1, 4] + S[2, 4]) * k13 +
              (S[0, 1] + S[1, 1] + S[1, 2]) * k22 +
              (S[0, 3] + S[1, 3] + S[2, 3]) * k23 +
              (S[0, 2] + S[1, 2] + S[2, 2]) * k33)
        return aa

    def calc_linear_compressibility_2D(self, S, theta):
        """
        Calculate linear compressibility for 2D materials.

        Parameters:
        S (numpy.ndarray): Compliance matrix for 2D material.
        theta (float): Angle theta in radians.

        Returns:
        float: Linear compressibility.
        """
        # Calculate components of the wave vector in 2D
        k1 = np.cos(theta)
        k2 = np.sin(theta)

        # Calculate k products for 2D
        k11 = k1 * k1
        k12 = k1 * k2
        k22 = k2 * k2

        # Calculate linear compressibility for 2D
        aa = ((S[0, 0] + S[0, 1]) * k11 +
              (S[0, 2] + S[1, 2]) * k12 +
              (S[0, 1] + S[1, 1]) * k22)
        return aa
        

    def plot_linear_compressibility_3D(self, fname='linear_compressibility', dpi=200):

        phi, theta = np.mgrid[0:2*np.pi:150j, 0:np.pi:100j]
        x = np.sin(theta) * np.cos(phi)
        y = np.sin(theta) * np.sin(phi)
        z = np.cos(theta)

        compressibility_3D = np.zeros_like(x)
        compressibility_2D = np.zeros_like(x)
        for i in range(x.shape[0]):
            for j in range(x.shape[1]):
                direction = np.array([x[i, j], y[i, j], z[i, j]])
                compressibility_3D[i, j] = self.linear_compressibility(direction)

                direction_c = np.array([x[i, j], y[i, j], np.sqrt(1 - x[i, j]**2 - y[i, j]**2)])
                compressibility_2D[i, j] = self.linear_compressibility(direction_c)

        fsize = 16
        fig = plt.figure(figsize=(16, 8))
        
        norm = plt.Normalize(vmin=min(compressibility_3D.min(), compressibility_2D.min()), 
                         vmax=max(compressibility_3D.max(), compressibility_2D.max()))

        # Now apply the colormap directly to the normalized data
        facecolors = plt.cm.viridis(norm(compressibility_3D)) #(compressibility_normalized)


        ax1 = fig.add_subplot(121, projection='3d')
        surface = ax1.plot_surface(
            x, y, z, rstride=1, cstride=1, facecolors= facecolors, # plt.cm.viridis(norm(compressibility_3D)),
            linewidth=1, antialiased=True, shade=False
        )

    
        ax1.set_title('Linear Compressibility 3D Plot', fontsize=fsize)
        ax1.set_xlabel('X-axis', fontsize=fsize)
        ax1.set_ylabel('Y-axis', fontsize=fsize)
        ax1.set_zlabel('Z-axis', fontsize=fsize)
        ax1.tick_params(axis='x', labelsize=fsize)
        ax1.tick_params(axis='y', labelsize=fsize)
        ax1.tick_params(axis='z', labelsize=fsize)
        cbar1 = fig.colorbar(surface, ax=ax1, shrink=0.5, aspect=10, pad=0.1)
        cbar1.set_label(r'LC (GPa$^{-1}$)', fontsize=fsize)
        

        # 2D contour plot
        ax2 = fig.add_subplot(122)
        contourf_plot = ax2.contourf(x, y, compressibility_2D, levels=30, cmap=plt.cm.viridis, alpha=0.8, norm=norm)
        contour_plot = ax2.contour(x, y, compressibility_2D, levels=10, linewidths=0.5, norm=norm)
        cbar2 = fig.colorbar(contourf_plot, ax=ax2, orientation='vertical', aspect=10, shrink=0.75, pad=0.01)
        cbar2.set_label(r'LC (GPa$^{-1}$)', fontsize=fsize)
        ax2.clabel(contour_plot, inline=True, fontsize=8, fmt='%1.1e')        
        ax2.set_xlabel('X-axis', fontsize=fsize)
        ax2.set_ylabel('Y-axis', fontsize=fsize)
        ax2.set_title('Linear Compressibility Contour Plot', fontsize=fsize)
        ax2.set_aspect('equal', 'box')
        ax2.grid(False)

        plt.tight_layout()
        plt.savefig(f"{fname}.png", format='png', dpi=dpi)
        plt.close(fig)

          
    def calc_poisson_directional(self, direction, khi, theta):
        # Normalize the direction vector
        v1_dir, v2_dir, v3_dir = direction / np.linalg.norm(direction)

        # Calculate the v vectors based on the given khi and theta
        v1 = np.cos(theta) * np.cos(khi) * v1_dir - np.sin(khi) * v2_dir
        v2 = np.cos(theta) * np.sin(khi) * v1_dir + np.cos(khi) * v2_dir
        v3 = -np.sin(theta) * v1_dir

        v11 = v1 * v1
        v12 = v1 * v2
        v13 = v1 * v3
        v22 = v2 * v2
        v23 = v2 * v3
        v33 = v3 * v3
        S = self.Cs

        # Calculate a11 and a12 using the compliance matrix S and the v values
        a11 = S[0, 0] * v11 + S[1, 1] * v22 + S[2, 2] * v33 + 2 * (S[0, 1] * v12 + S[0, 2] * v13 + S[1, 2] * v23)
        a12 = S[0, 0] * v11 * v1_dir + (S[0, 1] * v22 + S[1, 0] * v11) * v2_dir + (S[0, 2] * v33 + S[2, 0] * v11) * v3_dir + \
              S[1, 1] * v22 * v2_dir + (S[1, 2] * v33 + S[2, 1] * v22) * v3_dir + \
              S[2, 2] * v33 * v3_dir + 2 * (S[0, 3] * v23 * v1_dir + S[0, 4] * v13 * v1_dir + S[0, 5] * v12 * v1_dir + \
              S[1, 3] * v23 * v2_dir + S[1, 4] * v13 * v2_dir + S[1, 5] * v12 * v2_dir + \
              S[2, 3] * v23 * v3_dir + S[2, 4] * v13 * v3_dir + S[2, 5] * v12 * v3_dir + \
              S[3, 3] * v23 + S[3, 4] * v13 + S[3, 5] * v12 + \
              S[4, 4] * v13 + S[4, 5] * v12 + \
              S[5, 5] * v12)

        ratio = -a12 / a11
  
        return ratio    


    def calc_poisson_directional_2D(self, direction, khi, theta):
        # Normalize the direction vector
        S = self.Cs
        
        v1_dir, v2_dir = direction / np.linalg.norm(direction)

        # Calculate the v vectors based on the given khi and theta for 2D
        v1 = np.cos(theta) * np.cos(khi) * v1_dir - np.sin(khi) * v2_dir
        v2 = np.cos(theta) * np.sin(khi) * v1_dir + np.cos(khi) * v2_dir

        v11 = v1 * v1
        v12 = v1 * v2
        v22 = v2 * v2

        # Calculate a11 and a12 using the compliance matrix S and the v values for 2D
        a11 = S[0, 0] * v11 + S[1, 1] * v22 + 2 * S[0, 1] * v12
        a12 = S[0, 0] * v11 * v1_dir + (S[0, 1] * v22 + S[1, 0] * v11) * v2_dir + \
              S[1, 1] * v22 * v2_dir + 2 * (S[0, 2] * v12 * v1_dir + S[1, 2] * v12 * v2_dir + S[2, 2] * v12)

        ratio = -a12 / a11

        return ratio

    def calc_young_from_direction_2D(self,direction, khi, theta):
        # Normalize the direction vector
        d1, d2 = direction / np.linalg.norm(direction)
        S = self.Cdim

        # Calculate the k vectors based on the given khi and theta
        k1 = np.sin(theta) * np.cos(khi)
        k2 = np.sin(theta) * np.sin(khi)

        # Define the components of the k vectors
        k11 = k1 * k1
        k12 = k1 * k2
        k22 = k2 * k2

        # Calculate the Young's modulus using the 2D elastic tensor S and k vectors
        young_modulus = (
            d1*d1*(k11*S[0][0] + k12*S[0][1] + k12*S[0][2]) +
            d2*d2*(k22*S[1][1] + k12*S[1][2]) +
            2*d1*d2*(k11*S[0][1] + k12*S[1][1] + k12*S[1][2])
        )

        return young_modulus
        
    def calc_shear_3D(self, direction, khi, theta):
        # Normalize the direction vector
        d1, d2, d3 = direction / np.linalg.norm(direction)
        S = self.Cs
        # Calculate the k vectors based on the given khi and theta
        k1 = np.sin(theta) * np.cos(khi)
        k2 = np.sin(theta) * np.sin(khi)
        k3 = -np.cos(theta)

        k11 = k1 * k1
        k12 = k1 * k2
        k13 = k1 * k3
        k22 = k2 * k2
        k23 = k2 * k3
        k33 = k3 * k3

        # Calculate the v vectors based on the given khi and theta
        v1 = np.cos(theta) * np.cos(khi) * d1 - np.sin(khi) * d2
        v2 = np.cos(theta) * np.sin(khi) * d1 + np.cos(khi) * d2
        v3 = -np.sin(theta) * d1

        v11 = v1 * v1
        v12 = v1 * v2
        v13 = v1 * v3
        v22 = v2 * v2
        v23 = v2 * v3
        v33 = v3 * v3

        a66 = (k11 * v11 * S[0, 0] + 2 * k12 * v12 * S[0, 1] + 2 * k13 * v13 * S[0, 2] +
                k22 * v22 * S[1, 1] + 2 * k23 * v23 * S[1, 2] +
                k33 * v33 * S[2, 2])

        # Additional terms
        a66 += ((k12 * v13 + k13 * v12) * S[0, 3] + (k11 * v13 + k13 * v11) * S[0, 4] +
                (k11 * v12 + k12 * v11) * S[0, 5] + (k22 * v23 + k23 * v22) * S[1, 3] +
                (k12 * v23 + k23 * v12) * S[1, 4] + (k22 * v12 + k12 * v22) * S[1, 5] +
                (k33 * v23 + k23 * v33) * S[2, 3] + (k33 * v13 + k13 * v33) * S[2, 4] +
                (k13 * v23 + k23 * v13) * S[2, 5] + 1/4 * (k22 * v33 + 2 * k23 * v23 + k33 * v22) * S[3, 3] +
                1/2 * (k12 * v33 + k23 * v13 + k13 * v23 + k33 * v12) * S[3, 4] +
                1/2 * (k12 * v23 + k22 * v13 + k13 * v22 + k23 * v12) * S[3, 5] +
                1/4 * (k11 * v33 + 2 * k13 * v13 + k33 * v11) * S[4, 4] +
                1/2 * (k11 * v23 + k12 * v13 + k13 * v12 + k23 * v11) * S[4, 5] +
                1/4 * (k11 * v22 + 2 * k12 * v12 + k22 * v11) * S[5, 5])

        modulus = 1 / (4 * a66)
        return modulus
 
    def calc_young_from_direction(self, direction, khi, theta):
        # Normalize the direction vector
        d1, d2, d3 = direction / np.linalg.norm(direction)
        S = self.Cdim

        # Calculate the k vectors based on the given khi and theta
        k1 = np.sin(theta) * np.cos(khi)
        k2 = np.sin(theta) * np.sin(khi)
        k3 = -np.cos(theta)

        # Define the components of the k vectors
        k11 = k1 * k1
        k12 = k1 * k2
        k13 = k1 * k3
        k22 = k2 * k2
        k23 = k2 * k3
        k33 = k3 * k3

        # Calculate the Young's modulus using the elastic tensor S, k vectors, and direction vector components
        young_modulus = (
            d1*d1*(k11*S[0][0] + k12*S[0][1] + k13*S[0][2] + k12*S[0][3] + k13*S[0][4] + k12*S[0][5]) +
            d2*d2*(k22*S[1][1] + k23*S[1][2] + k22*S[1][3] + k12*S[1][4] + k12*S[1][5]) +
            d3*d3*(k33*S[2][2] + k23*S[2][3] + k13*S[2][4] + k13*S[2][5]) +
            2*d1*d2*(k11*S[0][1] + k12*S[1][1] + k13*S[0][4] + k12*S[1][4]) +
            2*d1*d3*(k11*S[0][2] + k13*S[2][2] + k12*S[0][5] + k13*S[2][5]) +
            2*d2*d3*(k22*S[1][2] + k23*S[2][2] + k12*S[1][5] + k23*S[2][5])
        )

        return young_modulus
                
    def plot_directional_poisson_3D(self,fname='poisson_ratio_contour_directional_3D', dpi=150): 
      # Define directions [100], [110], and [111]
      direction_100 = np.array([1, 0, 0])
      direction_110 = np.array([-1, 1, 0]) / np.sqrt(2)
      direction_111 = np.array([1, 1, 1]) / np.sqrt(3)

      # Define the range of khi and theta values
      khi_vals = np.linspace(0, np.pi, 50)
      theta_vals = np.linspace(0, np.pi, 50)

      # Initialize grids to store Poisson ratios
      ratio_grid_100 = np.zeros((len(khi_vals), len(theta_vals)))
      ratio_grid_110 = np.zeros((len(khi_vals), len(theta_vals)))
      ratio_grid_111 = np.zeros((len(khi_vals), len(theta_vals)))

      # Calculate the Poisson ratios for each direction
      for i, khi in enumerate(khi_vals):
          for j, theta in enumerate(theta_vals):
              ratio_grid_100[i, j] = self.calc_poisson_directional(direction_100, khi, theta)
              ratio_grid_110[i, j] = self.calc_poisson_directional(direction_110, khi, theta)
              ratio_grid_111[i, j] = self.calc_poisson_directional(direction_111, khi, theta)

      fig = plt.figure(figsize=(18, 6))
      fsize = 16
      # Plot for [100] direction
      plt.subplot(1, 3, 1)
      X, Y = np.meshgrid(khi_vals, theta_vals)
      plt.contourf(X, Y, ratio_grid_100.T, levels=20, cmap='rainbow_r')
      plt.colorbar(label='Poisson Ratio')
      contour_lines1 = plt.contour(X, Y, ratio_grid_100.T, levels=20, colors='black')
      plt.clabel(contour_lines1, inline=True, fontsize=8)
      plt.xlabel(r'$\Phi$ (radians)',fontsize=fsize)
      plt.ylabel(r'$\Theta$ (radians)',fontsize=fsize)
      plt.xticks(fontsize=fsize)
      plt.yticks(fontsize=fsize)
      plt.title('Poisson Ratio for [100] Direction')

      # Plot for [110] direction
      plt.subplot(1, 3, 2)
      plt.contourf(X, Y, ratio_grid_110.T, levels=20, cmap='rainbow_r')
      plt.colorbar(label='Poisson Ratio')
      contour_lines2 = plt.contour(X, Y, ratio_grid_110.T, levels=20, colors='black')
      plt.clabel(contour_lines2, inline=True, fontsize=8)
      plt.xlabel(r'$\Phi$ (radians)',fontsize=fsize)
      plt.ylabel(r'$\Theta$ (radians)',fontsize=fsize)
      plt.xticks(fontsize=fsize)
      plt.yticks(fontsize=fsize)
      plt.title(r'Poisson Ratio for [$\bar{1}$10] Direction')

      # Plot for [111] direction
      plt.subplot(1, 3, 3)
      plt.contourf(X, Y, ratio_grid_111.T, levels=20, cmap='rainbow_r')
      plt.colorbar(label='Poisson Ratio')
      contour_lines3 = plt.contour(X, Y, ratio_grid_111.T, levels=20, colors='black')
      plt.clabel(contour_lines3, inline=True, fontsize=8)
      plt.xlabel(r'$\Phi$ (radians)',fontsize=fsize)
      plt.ylabel(r'$\Theta$ (radians)',fontsize=fsize)
      plt.xticks(fontsize=fsize)
      plt.yticks(fontsize=fsize)      
      plt.title('Poisson Ratio for [111] Direction')

      plt.tight_layout()
      #plt.show()       
      plt.tight_layout()
      plt.savefig(f"{fname}.png", format='png', dpi=dpi)
      plt.close(fig)
      
      
      
    def plot_directional_poisson_2D(self, fname='poisson_ratio_contour_directional_2D', dpi=200):
        direction_10 = np.array([1, 0])
        direction_1n1 = np.array([0, -1])
        direction_11 = np.array([1, 1])/np.sqrt(2)

        # Define the range of khi and theta values
        khi_vals = np.linspace(0, np.pi, 50)
        theta_vals = np.linspace(0, np.pi, 50)
        ratio_grid_10 = np.zeros((len(khi_vals), len(theta_vals)))
        ratio_grid_0n1 = np.zeros((len(khi_vals), len(theta_vals)))
        ratio_grid_11 = np.zeros((len(khi_vals), len(theta_vals)))

        # Calculate the Poisson ratio for each combination of khi and theta
        for i, khi in enumerate(khi_vals):
            for j, theta in enumerate(theta_vals):
                ratio_10 = self.calc_poisson_directional_2D(direction_10, khi, theta)
                ratio_0n1 = self.calc_poisson_directional_2D(direction_1n1, khi, theta)
                ratio_11 = self.calc_poisson_directional_2D(direction_11, khi, theta)
                ratio_grid_10[i, j] = ratio_10
                ratio_grid_0n1[i, j] = ratio_0n1
                ratio_grid_11[i, j]  = ratio_11

        fig = plt.figure(figsize=(18, 6))
        fsize=16
        # Plot for [100] direction
        plt.subplot(1, 3, 1)
        X, Y = np.meshgrid(khi_vals, theta_vals)
        plt.contourf(X, Y, ratio_grid_10.T, levels=20, cmap='rainbow_r')
        plt.colorbar(label='Poisson Ratio')
        contour_lines1 = plt.contour(X, Y, ratio_grid_10.T, levels=10, colors='black')
        plt.clabel(contour_lines1, inline=True, fontsize=8)
        plt.xlabel(r'$\Phi$ (radians)',fontsize=fsize)
        plt.ylabel(r'$\Theta$ (radians)',fontsize=fsize)
        plt.xticks(fontsize=fsize)
        plt.yticks(fontsize=fsize)
        plt.title('Poisson Ratio for [100] Direction')

        # Plot for [-11] direction
        plt.subplot(1, 3, 2)
        plt.contourf(X, Y, ratio_grid_0n1.T, levels=20, cmap='rainbow_r')
        plt.colorbar(label='Poisson Ratio')
        contour_lines2 = plt.contour(X, Y, ratio_grid_0n1.T, levels=10, colors='black')
        plt.clabel(contour_lines2, inline=True, fontsize=8)
        plt.xlabel(r'$\Phi$ (radians)')
        plt.ylabel(r'$\Theta$ (radians)')
        plt.title(r'Poisson Ratio for [0$\bar{1}$0] Direction')

        # Plot for [111] direction
        plt.subplot(1, 3, 3)
        plt.contourf(X, Y, ratio_grid_11.T, levels=20, cmap='rainbow_r')
        plt.colorbar(label='Poisson Ratio')
        contour_lines3 = plt.contour(X, Y, ratio_grid_11.T, levels=10, colors='black')
        plt.clabel(contour_lines3, inline=True, fontsize=8)
        plt.xlabel(r'$\Phi$ (radians)')
        plt.ylabel(r'$\Theta$ (radians)')
        plt.title('Poisson Ratio for [110] Direction')

        plt.tight_layout()
        plt.savefig(f"{fname}.png", format='png', dpi=dpi)
        plt.close(fig)
                
        
    def plot_poisson_3Dprojection_3D(self,fname='poisson_ratio_3Dprojection_3D', dpi=150):
        direction_100 = np.array([1, 0, 0])
        direction_110 = np.array([-1, 1, 0]) / np.sqrt(2)
        direction_111 = np.array([1, 1, 1]) / np.sqrt(3)
        # Define the range of khi and theta values
        khi_vals = np.linspace(0, np.pi, 50)
        theta_vals = np.linspace(0, np.pi, 50)

        # Initialize grids to store Poisson ratios
        ratio_grid_100 = np.zeros((len(khi_vals), len(theta_vals)))
        ratio_grid_110 = np.zeros((len(khi_vals), len(theta_vals)))
        ratio_grid_111 = np.zeros((len(khi_vals), len(theta_vals)))

        # Calculate the Poisson ratios for each direction
        for i, khi in enumerate(khi_vals):
            for j, theta in enumerate(theta_vals):
                ratio_grid_100[i, j] = self.calc_poisson_directional(direction_100, khi, theta)
                ratio_grid_110[i, j] = self.calc_poisson_directional(direction_110, khi, theta)
                ratio_grid_111[i, j] = self.calc_poisson_directional(direction_111, khi, theta)

        fig = plt.figure(figsize=(18, 6))
        fsize = 16
        try:
            # Plot for [100] direction
            X, Y = np.meshgrid(khi_vals, theta_vals)
            ax1 = plt.subplot(1, 3, 1, projection='3d')
            surf1 = ax1.plot_surface(X, Y, ratio_grid_100, cmap='viridis', edgecolor='none')
            ax1.contour(X, Y, ratio_grid_100, zdir='z', offset=ratio_grid_100.min(), levels=20, cmap='viridis')
            ax1.set_xlabel(r'$\Phi$ (radians)', fontsize=fsize)
            ax1.set_ylabel(r'$\Theta$ (radians)', fontsize=fsize)
            ax1.set_zlabel('Poisson Ratio', fontsize=fsize)
            ax1.set_title('3D Poisson Ratio for [100] Direction')

            # Plot for [-110] direction
            ax2 = plt.subplot(1, 3, 2, projection='3d')
            surf2 = ax2.plot_surface(X, Y, ratio_grid_110, cmap='viridis', edgecolor='none')
            ax2.contour(X, Y, ratio_grid_110, zdir='z', offset=ratio_grid_110.min(), levels=20, cmap='viridis')
            ax2.set_xlabel(r'$\Phi$ (radians)', fontsize=fsize)
            ax2.set_ylabel(r'$\Theta$ (radians)', fontsize=fsize)
            ax2.set_zlabel('Poisson Ratio', fontsize=fsize)
            ax2.set_title(r'3D Poisson Ratio for [$\bar{1}$10] Direction')

            # Plot for [111] direction
            ax3 = plt.subplot(1, 3, 3, projection='3d')
            surf3 = ax3.plot_surface(X, Y, ratio_grid_111, cmap='viridis', edgecolor='none')
            ax3.contour(X, Y, ratio_grid_111, zdir='z', offset=ratio_grid_111.min(), levels=20, cmap='viridis')
            ax3.set_xlabel(r'$\Phi$ (radians)', fontsize=fsize)
            ax3.set_ylabel(r'$\Theta$ (radians)', fontsize=fsize)
            ax3.set_zlabel('Poisson Ratio', fontsize=fsize)
            ax3.set_title('3D Poisson Ratio for [111] Direction')
          
            plt.tight_layout()
            plt.savefig(f"{fname}.png", format='png', dpi=dpi)
            plt.close(fig)    
        except ValueError as e:
            print(f"Failed to plot poisson_ratio_3Dprojection_3D Error: {e}")   
            

    def plot_shear_3Dprojection_3D(self,fname='shearmodulus_3Dprojection_3D', dpi=150):
        direction_100 = np.array([1, 0, 0])
        direction_110 = np.array([-1, 1, 0]) / np.sqrt(2)
        direction_111 = np.array([1, 1, 1]) / np.sqrt(3)
        # Define the range of khi and theta values
        khi_vals = np.linspace(0, np.pi, 50)
        theta_vals = np.linspace(0, np.pi, 50)

        # Initialize grids to store Shear Moduluss
        ratio_grid_100 = np.zeros((len(khi_vals), len(theta_vals)))
        ratio_grid_110 = np.zeros((len(khi_vals), len(theta_vals)))
        ratio_grid_111 = np.zeros((len(khi_vals), len(theta_vals)))

        # Calculate the Shear Moduluss for each direction
        for i, khi in enumerate(khi_vals):
            for j, theta in enumerate(theta_vals):
                ratio_grid_100[i, j] = self.calc_shear_3D(direction_100, khi, theta)
                ratio_grid_110[i, j] = self.calc_shear_3D(direction_110, khi, theta)
                ratio_grid_111[i, j] = self.calc_shear_3D(direction_111, khi, theta)

        fig = plt.figure(figsize=(20, 6))
        fsize = 16
        try:
            # Plot for [100] direction
            X, Y = np.meshgrid(khi_vals, theta_vals)
            ax1 = plt.subplot(1, 3, 1, projection='3d')
            surf1 = ax1.plot_surface(X, Y, ratio_grid_100, cmap='viridis', edgecolor='none')
            ax1.contour(X, Y, ratio_grid_100, zdir='z', offset=ratio_grid_100.min(), levels=20, cmap='viridis')
            ax1.set_xlabel(r'$\Phi$ (radians)', fontsize=fsize)
            ax1.set_ylabel(r'$\Theta$ (radians)', fontsize=fsize)
            ax1.set_zlabel('Shear Modulus', fontsize=fsize)
            ax1.set_title('3D Shear Modulus for [100] Direction')

            # Plot for [-110] direction
            ax2 = plt.subplot(1, 3, 2, projection='3d')
            surf2 = ax2.plot_surface(X, Y, ratio_grid_110, cmap='viridis', edgecolor='none')
            ax2.contour(X, Y, ratio_grid_110, zdir='z', offset=ratio_grid_110.min(), levels=20, cmap='viridis')
            ax2.set_xlabel(r'$\Phi$ (radians)', fontsize=fsize)
            ax2.set_ylabel(r'$\Theta$ (radians)', fontsize=fsize)
            ax2.set_zlabel('Shear Modulus', fontsize=fsize)
            ax2.set_title(r'3D Shear Modulus for [$\bar{1}$10] Direction')

            # Plot for [111] direction
            ax3 = plt.subplot(1, 3, 3, projection='3d')
            surf3 = ax3.plot_surface(X, Y, ratio_grid_111, cmap='viridis', edgecolor='none')
            ax3.contour(X, Y, ratio_grid_111, zdir='z', offset=ratio_grid_111.min(), levels=20, cmap='viridis')
            ax3.set_xlabel(r'$\Phi$ (radians)', fontsize=fsize)
            ax3.set_ylabel(r'$\Theta$ (radians)', fontsize=fsize)
            ax3.set_zlabel('Shear Modulus', fontsize=fsize)
            ax3.set_title('3D Shear Modulus for [111] Direction')
          
            plt.tight_layout()
            plt.savefig(f"{fname}.png", format='png', dpi=dpi)
            plt.close(fig)    
        except ValueError as e:
            print(f"Failed to plot shearmodulus_3Dprojection_3D Error: {e}")     


    def plot_E_3Dprojection_3D(self,fname='youngmodulus_3Dprojection_3D', dpi=150):
        
        direction_100 = np.array([1, 0, 0])
        direction_110 = np.array([-1, 1, 0]) / np.sqrt(2)
        direction_111 = np.array([1, 1, 1]) / np.sqrt(3)
        # Define the range of khi and theta values
        khi_vals = np.linspace(0, np.pi, 50)
        theta_vals = np.linspace(0, np.pi, 50)

        # Initialize grids to store Young Moduluss
        ratio_grid_100 = np.zeros((len(khi_vals), len(theta_vals)))
        ratio_grid_110 = np.zeros((len(khi_vals), len(theta_vals)))
        ratio_grid_111 = np.zeros((len(khi_vals), len(theta_vals)))
        
        # Calculate the Young Moduluss for each direction
        for i, khi in enumerate(khi_vals):
            for j, theta in enumerate(theta_vals):
                ratio_grid_100[i, j] = self.calc_young_from_direction(direction_100, khi, theta)
                ratio_grid_110[i, j] = self.calc_young_from_direction(direction_110, khi, theta)
                ratio_grid_111[i, j] = self.calc_young_from_direction(direction_111, khi, theta)

        fig = plt.figure(figsize=(19, 6))
        fsize = 16

        try:
            # Plot for [100] direction
            X, Y = np.meshgrid(khi_vals, theta_vals)
            ax1 = plt.subplot(1, 3, 1, projection='3d')
            surf1 = ax1.plot_surface(X, Y, ratio_grid_100, cmap='rainbow', edgecolor='none')
            ax1.contour(X, Y, ratio_grid_100, zdir='z', offset=ratio_grid_100.min(), levels=20, cmap='rainbow')
            ax1.set_xlabel(r'$\Phi$ (radians)', fontsize=fsize)
            ax1.set_ylabel(r'$\Theta$ (radians)', fontsize=fsize)
            ax1.set_zlabel(r'Young$\'s$ Modulus', fontsize=fsize)
            ax1.set_title(r'Young$\'s$ Modulus for [100] Direction')

            # Plot for [-110] direction
            ax2 = plt.subplot(1, 3, 2, projection='3d')
            surf2 = ax2.plot_surface(X, Y, ratio_grid_110, cmap='rainbow', edgecolor='none')
            ax2.contour(X, Y, ratio_grid_110, zdir='z', offset=ratio_grid_110.min(), levels=20, cmap='rainbow')
            ax2.set_xlabel(r'$\Phi$ (radians)', fontsize=fsize)
            ax2.set_ylabel(r'$\Theta$ (radians)', fontsize=fsize)
            ax2.set_zlabel(r'Young$\'s$ Modulus', fontsize=fsize)
            ax2.set_title(r'Young$\'s$ Modulus for [$\bar{1}$10] Direction')

            # Plot for [111] direction
            ax3 = plt.subplot(1, 3, 3, projection='3d')
            surf3 = ax3.plot_surface(X, Y, ratio_grid_111, cmap='rainbow', edgecolor='none')
            ax3.contour(X, Y, ratio_grid_111, zdir='z', offset=ratio_grid_111.min(), levels=20, cmap='rainbow')
            ax3.set_xlabel(r'$\Phi$ (radians)', fontsize=fsize)
            ax3.set_ylabel(r'$\Theta$ (radians)', fontsize=fsize)
            ax3.set_zlabel(r'Young$\'s$ Modulus', fontsize=fsize)
            ax3.set_title(r'Young$\'s$ Modulus for [111] Direction')
          
            plt.tight_layout()
            plt.savefig(f"{fname}.png", format='png', dpi=dpi)
            plt.close(fig)    
        except ValueError as e:
            print(f"Failed to plot Youngmodulus_3Dprojection_3D Error: {e}")     


    def plot_combined_directional_2D(self, fname='EV_ploar_directional_2D', dpi=200):
        direction_100 = np.array([1, 0])
        direction_010 = np.array([0, 1])
        direction_110 = np.array([-1, 1]) / np.sqrt(2)

        # Define the range of angles
        angles = np.linspace(0, 2 * np.pi, 100)

        # Initialize arrays for Poisson's ratio and Young's modulus
        poisson_moduli_100, poisson_moduli_010, poisson_moduli_110 = [], [], []
        young_moduli_100, young_moduli_010, young_moduli_110 = [], [], []

        for angle in angles:
            # Calculate Poisson's ratio for each direction
            poisson_ratio_100 = self.calc_poisson_directional_2D(direction_100, angle, angle)
            poisson_ratio_010 = self.calc_poisson_directional_2D(direction_010, angle, angle)
            poisson_ratio_110 = self.calc_poisson_directional_2D(direction_110, angle, angle)
            poisson_moduli_100.append(poisson_ratio_100)
            poisson_moduli_010.append(poisson_ratio_010)
            poisson_moduli_110.append(poisson_ratio_110)

            # Calculate Young's modulus for each direction
            young_modulus_100 = self.calc_young_from_direction_2D(direction_100, angle, angle)
            young_modulus_010 = self.calc_young_from_direction_2D(direction_010, angle, angle)
            young_modulus_110 = self.calc_young_from_direction_2D(direction_110, angle, angle)
            young_moduli_100.append(young_modulus_100)
            young_moduli_010.append(young_modulus_010)
            young_moduli_110.append(young_modulus_110)

        # Create subplots
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6), subplot_kw={'projection': 'polar'})

        # Plot Poisson's ratio
        ax1.plot(angles, poisson_moduli_100, label=r'$\nu$ [100]')
        ax1.plot(angles, poisson_moduli_010, label=r'$\nu$ [010]')
        ax1.plot(angles, poisson_moduli_110, label=r'$\nu$ [$\bar{1}$10]')
        ax1.set_title('Poisson Ratio in Different Directions')
        ax1.legend()

        # Plot Young's modulus
        ax2.plot(angles, young_moduli_100, label=r'$E$ [100]')
        ax2.plot(angles, young_moduli_010, label=r'$E$ [010]')
        ax2.plot(angles, young_moduli_110, label=r'$E$ [$\bar{1}$10]')
        ax2.set_title('Young\'s Modulus in Different Directions')
        ax2.legend()

        # Adjust layout and save
        plt.tight_layout()
        plt.savefig(f"{fname}.png", format='png', dpi=dpi)
        plt.show()
                                
# Create an instance of the class and plot the heatmaps
#analyzer = ElasticParametersHeatMap(elastic_tensor)
#analyzer.plot_moduli_heatmaps()



