# ElasTool

**ElasTool** (*Elastic Analysis and Simulation Toolkit for Optimized Observations and Learning*) is a powerful and innovative Python-based toolkit optimized for computing the second-order elastic constants (SOECs) and mechanical properties of crystal systems across different dimensions, including 3D structures, 2D materials, nanoribbons, coaxially rolled 2D-based van der Waals nanostructures, and 1D nanotubes. It harnesses first-principles Density Functional Theory (DFT) and ab initio molecular dynamics, making it ideal for analyzing both zero-temperature and finite-temperature conditions.


##  About ElasTool

ElasTool offers a flexible approach to determining elastic constants and mechanical properties of various materials at zero and finite temperatures and pressures, providing a broad scope of utility across different material conditions. It seamlessly integrates with the VASP electronic structure code. However, its architecture is open and can easily implement interfaces to other DFT packages. If you are interested in such extensions, please don't hesitate to contact the authors for further guidance and support with the ElasTool source code. The software utilizes three kinds of strain-matrix sets: High-Efficiency Strain-Matrix Sets (OHESS), Universal Linear-Independent Coupling Strains (ULICS), and All-Single-Element Strain-Matrix Sets (ASESS) [1], enabling automatic and efficient calculation of the SOECs.

##  Key Features

The ElasTool toolkit offers several features, including:
- Easy to use (installation and running)
- High efficiency
- Automated flow of the SOECs calculations
- Choice of three kinds of strain-matrix sets: OHESS, ASESS, and ULICS
- Capability to compute SOECs at zero temperature
- Capability to compute SOECs at high-temperature and/or high-pressure
- Ability to calculate the elastic and mechanical properties of 3D, 2D, coaxially rolled 2D-based van der Waals, and 1D nanotube materials.


## Expanded Capabilities

ElasTool offers a comprehensive range of features for analyzing various material properties:

- **Elastic Tensors**: Generates a detailed text file for in-depth post-processing.
- **Young's Modulus**: Evaluates the material's stiffness under uniaxial stress.
- **Poisson Ratio**: Determines the ratio of transverse strain to axial strain.
- **Shear Modulus**: Assesses the material's response to shear stress.
- **Bulk Modulus (3D/1D), Stiffness Constant (2D)**: Measures volume change under pressure and in-plane stiffness, respectively.
- **Pugh Modulus Ratio**:  Provides insights into the material's ductility.
- **Layer Modulus (2D Materials)**: Evaluates in-plane elasticity of layers.
- **Sound Velocities**: Longitudinal, transverse, and average sound velocities.
- **Debye Speed**: Estimates phonon propagation speeds.
- **Linear Compressibility**: Assesses the material's response to linear compressive stress.
- **Debye Temperature**: Evaluates the material's thermal properties.
- **Minimum Thermal Conductivity**: Utilizes Clarke and Cahill equations to estimate thermal conductivity limits.
- **Strain Energy Density**: Provides insights into the material's energy absorption capacity.
- **Hardness Estimation**: Employs various empirical equations to predict Vickers hardness.
- **Fracture Toughness Analysis**: Evaluates the material's resistance to crack propagation.
- **Elastic Anisotropy**: Examines directional variations in material properties.
- **Compliance Matrix Eigenvalues**: Investigates the material's mechanical response characteristics.
- **Solution of the Christoffel Equations***: 

The Christoffel equations have been integrated into our approach for calculating wave velocities, power factors, and enhancement factors. This integration is based on the methodology described in [Christoffel](https://www.sciencedirect.com/science/article/pii/S0010465516301795). The feature is turned on once `plotparameters` is set to `yes` in the ElasTool input file `elastol.in`. For 2D materials, we implemented an embedding method to visualize the properties in a three-dimensional plane. Specifically, the solution of the Christoffel equations generates the following output files: 

- `slow_secondary.dat`, `fast_secondary.dat`, `primary.dat` 
  These files provide acoustic data for various directions, determined by the `Numtheta` and `Numphi` grid settings. They offer insights into the slowest, intermediate, and fastest acoustic modes. Each file comprises 17 columns with the following data:
    1. Theta (rad)
    2. Phi (rad)
    3,4,5. Cartesian projection on the unit cube
    6. Phase velocity (km/s)
    7. Phase velocity compared to isotropic velocity (%)
    8,9,10. Phase polarization (unit vector, dimensionless)
    11. Group velocity (km/s)
    12. Group velocity compared to isotropic velocity (%)
    13,14,15. Cartesian coordinates of the Ray surface (km/s)
    16. Power flow angle (deg)
    17. Enhancement factor (dimensionless)

- `anisotropy.dat` 
  This file details the maximum and minimum velocities for the three acoustic modes and their directional occurrences. It also quantifies the material's overall anisotropy.

- `directions.dat` 
  Contains information on sound velocities in various directions.

- `sound.out` 
  Offers general data on the etensor tensor, bulk and shear modulus, and the isotropic sound velocities of the material.

### 3.4) Visualization
Upon enabling the plot feature, ElasTool generates Gnuplot scripts for automatic plotting of key features. These are saved as high-resolution PNG files in the `property_plots` folder. Due to the volume of files produced, enabling plot features in high-throughput settings is not recommended unless necessary. Alternatively, the post-processing mode (`run_mode=3`) can be used for plot generation for specific materials. For outputs related to the Christoffel equation solution, the following file naming conventions are important:

- `phase` ----- Phase velocity
- `group` ----- Group velocity (absolute value only)
- `ray` ------- Ray surface (directional group velocity)
- `pfangle` --- Power flow angle
- `enh` ------- Enhancement factor
- `pol` ------- Phase velocity polarization
- `relative` -- Relative to isotropic sound velocities (for phase and group)
- `sphere` ---- Projection onto the unit sphere
- `cube` ------ Projection onto the unit cube
- `eqar` ------ Equal area plane projection
- `stereo` ---- Stereographic plane projection, preserving local shape
- `radius` ---- For phase velocity, showing data on a sphere with radius scaled by absolute velocity

These features establish ElasTool as a versatile and indispensable toolkit for scientists and engineers specializing in materials science and engineering. Its comprehensive computational capabilities are ideal for advanced material science research, providing critical insights into the mechanical and thermal properties of various materials. Whether for academic research, industrial applications, or innovative material design, ElasTool serves as a pivotal resource in exploring and understanding the intricate behavior of 3D structures, 2D materials, and nanoscale systems.


## ElasTool: Advanced Data Visualization and Analysis Toolkit

ElasTool is a highly efficient computational tool designed for calculating and visualizing the elastic, mechanical, and related properties of materials in 1D, 2D, and 3D systems. It stands out with its advanced visualization capabilities, offering insightful and enriching user experiences. **The visualization capabilities in ElasTool** are designed to enhance user engagement, provide deeper insights into material properties, and facilitate the efficient presentation of complex data. Whether for academic research, material design, or engineering applications, these features make ElasTool a valuable asset in the field of material science.


### Advanced Visualization Capabilities
ElasTool is not just a computational toolkit; it's also a powerful visualization platform, enabling:
- Integration with external data: ElasTool can import elastic tensor matrices computed externally, supporting:
  - 6x6 matrix for 3D systems
  - 3x3 matrix for 2D systems
- Required files:
  - `massdensity_dim.dat` for material's dimension and mass density.
    ```
    # Mass density in Kg/m^2, Dimension
    0.00000224 2D
    ```
  - `elastic_tensor.dat` for the elastic tensor matrix.
    ```
    # Elastic tensor in Voigt notation for 2D material
    52.2849 28.6494 0.0000
    28.6494 36.5780 0.0000
    0.0000 0.0000 22.8516
    ```

### Post-Processing Capabilities
The postprocessing is for generating the plots for visualizing the elastic and mechanical properties of materials after a complete stress-strain calculation has been carried out. The `massdensity_dim.dat` and the `elastic_tensor.dat` files are only needed if you're using elastic tensor parameters from other electronic structure codes or besides the ones you have obtained. To utilize ElasTool's post-processing features:
1. Set up the `massdensity_dim.dat` and `elastic_tensor.dat` files as described.
2. Use the following terminal commands (case insensitive):
    ```
    elastool -pp
    elastool -Postprocess
    elastool -postprocess
    elastool -POSTPROCESS
    ```
3. For interactive web plots with Plotly, add `-plotly`:
    ```
    elastool -pp -plotly
    elastool -Postprocess -plotly
    elastool -postprocess -plotly
    elastool -POSTPROCESS -plotly
    ```

### Automated Plotting and Visualization
ElasTool excels in:
- **Angular Dependence Visualization:**  ElasTool automatically generates plots depicting the angular dependence of elastic and mechanical properties. These plots are saved as PNG files directly in the run directory, offering a convenient way to review and share results.
- **Elastic Parameters Heatmaps:** For materials with anisotropic properties, ElasTool plots heatmaps illustrating variations in key elastic moduli, such as Young's Modulus, Shear Modulus, and Poisson's Ratio.
- **Sound Velocity Heatmaps:** ElasTool uses the Christoffel matrix to visualize sound velocity in materials in both transverse and longitudinal directions, crucial for acoustic properties.
- **Spatial Dependence of Elastic Parameters:** Plots the spatial variation of key elastic parameters, providing a comprehensive view of material behavior.
- **Linear Compressibility Analysis:** Offers detailed insights into the linear compressibility of materials, showcasing directional dependencies and anisotropic behaviors. This feature enhances the understanding of how materials compress under linear stress, crucial for applications in material science and engineering.
- **Spatial Distribution of the Strain Energy Distribution Function (SEDF):** ElasTool provides a detailed visualization of the SEDF, offering insights into the material's response to various strain states. This feature allows users to explore how the energy density changes spatially within the material under different mechanical deformations, enhancing the understanding of the material's mechanical properties and behavior under stress.


### Interactive Interface
- **Plotly Integration:** Offers an interactive interface for an engaging and detailed analysis of elastic parameters using plotly.

### Flexibility in Usage
- **Runtime and Post-Processing Accessibility:** ElasTool's visualization tools are accessible both during runtime and post-processing, with options to enable or disable plotting as needed.


### Seamless Integration with Elate
You can also use ElasTool for post-processing visualization using the Elate web interface. After completing stress-strain calculations, users can seamlessly integrate ElasTool with Elate.
- To do this, execute the command `elastool -elate`, or `elastool -Elate`, or `elastool -ELATE` in your terminal.
- ElasTool facilitates a smooth workflow by automatically fetching the elastic tensor data from your work directory for Elate analysis, eliminating the need for manual data transfer.
- During this process, you will be prompted to select your default web browser for optimal visualization and interaction.

These features enhance the user experience by simplifying the analysis process and providing intuitive, accessible data visualizations directly from ElasTool's interface.


## Installation

ElasTool offers straightforward installation options suitable for various user preferences. These methods ensure a hassle-free setup, allowing you to commence your material science investigations with ElasTool promptly. Detailed instructions can be found in the INSTALL file, but here are the general methods:

1. **Using pip**:
   - Quickly install ElasTool with pip by executing: 
     ```
     pip install -U elastool
     ```

2. **From Source Code**:
   - Alternatively, download the source code with:
     ```
     git clone [git@github.com:zhongliliu/elastool.git]
     ```
   - Then, install ElasTool by navigating to the master directory and running:
     ```
     pip install .
     ```

3. **Installation via setup.py**:
   - As an alternative, ElasTool can be installed using the `setup.py` script:
     ```
     python setup.py install [--prefix=/path/to/install/]
     ```
   - The optional `--prefix` argument is useful for installations in environments like shared High-Performance Computing (HPC) systems, where administrative privileges might be restricted.
   - Please note that while this method remains supported, its usage is gradually declining in favor of more modern installation practices. It is recommended primarily for specific scenarios where standard installation methods like `pip` are not applicable.

     
## Usage and Running ElasTool

Learning to use ElasTool is made easy with the provided examples in the example folder. Here are the key steps for using ElasTool effectively:

1. **Create a Calculation Directory**:
   - Start by creating a directory for your calculations.
   - Run `elastool -0` to generate basic template input files (INCARs, KPOINTS, and elastool.in).

2. **Modify Input Files**:
   - Customize the generated files according to your project's requirements.

3. **Initialize the Job**:
   - Execute `elastool` to begin the calculation process.

4. **Understanding ElasTool Options**:
   - The main input file `elastool.in` includes descriptive text for each flag, making it user-friendly.
   - For additional help or to explore more features, use `elastool -h` or `elastool -help`.

This streamlined process allows users to quickly start and efficiently conduct their material analysis with ElasTool.


## Citing ElasTool

If you have used ElasTool in your research, please, consider citing the appropriate publications from the list below:

### Main ElasTool Implementation
- Please cite for ElasTool's primary implementation:
  - [ElasTool: An automated toolkit for elastic constants calculation](https://doi.org/10.1016/j.cpc.2021.108180) - Liu et al., 2022

@article{Liu2020elastool,
  title = {ElasTool: An automated toolkit for elastic constants calculation},
  journal = {Computer Physics Communications},
  volume = {270},
  pages = {108180},
  year = {2022},
  issn = {0010-4655},
  doi = {https://doi.org/10.1016/j.cpc.2021.108180},
  url = {https://www.sciencedirect.com/science/article/pii/S0010465521002927},
  author = {Zhong-Li Liu and C.E. Ekuma and Wei-Qi Li and Jian-Qun Yang and Xing-Ji Li}
}

  - [Efficient prediction of temperature-dependent elastic and mechanical properties of 2D materials](https://www.nature.com/articles/s41598-022-07819-8) - Kastuar et al., 2022
  
@article{Kastuar2022efficient,
  title={Efficient prediction of temperature-dependent elastic and mechanical properties of 2D materials},
  author={Kastuar, SM and Ekuma, CE and Liu, Z-L},
  journal={Scientific Reports},
  volume={12},
  number={1},
  pages={3776},
  year={2022},
  url = {https://www.nature.com/articles/s41598-022-07819-8},
  publisher={Nature Publishing Group UK London}
}

  - [ElasTool v3.0: Efficient computational and visualization toolkit for elastic and mechanical properties of materials](https://www.sciencedirect.com/science/article/abs/pii/S0010465524000845?via%3Dihub) - Ekuma and Liu 

@article{Ekuma2024,
  title = {ElasTool v3.0: Efficient computational and visualization toolkit for elastic and mechanical properties of materials},
  journal = {Computer Physics Communications},
  volume = {300},
  pages = {109161},
  year = {2024},
  doi = {10.1016/j.cpc.2024.109161},
  url = {https://www.sciencedirect.com/science/article/abs/pii/S0010465524000845?via%3Dihub},
  author = {Chinedu E. Ekuma and Zhong-Li Liu }
}

### Work Related to 2D Materials
- For work specifically on 2D materials, refer to:
  - [Efficient prediction of temperature-dependent elastic and mechanical properties of 2D materials](https://www.nature.com/articles/s41598-022-06650-1) - Kastuar et al., 2022

### Work Related to Tubular 2D-Based Nanostructures and Nanotubes, and Advanced Visualization
- For studies on tubular 2D-based nanostructures and nanotubes, please consider citing:
  - [ElasTool v3.0: Efficient computational and visualization toolkit for elastic and mechanical properties of materials](https://www.sciencedirect.com/science/article/abs/pii/S0010465524000845?via%3Dihub) - Ekuma and Liu

### OHESS Method and Strain-Stress Methods
- For the OHESS method and other strain-stress methods, please cite:
  - [Investigating elastic constants across diverse strain-matrix sets](https://doi.org/10.1016/j.commatsci.2023.112521) - Liu et al., 2023


### Related Works
- For related research, please consider citing:
  - [Calculations of single-crystal elastic constants made simple](https://doi.org/10.1016/j.cpc.2009.11.017) - Yu et al., 2010
  - [Mechanical properties and hardness of boron pnicogens BX](https://doi.org/10.1016/j.mtla.2020.100904) - Ekuma and Liu, 2020


## Contact Information

We welcome your interest in extending ElasTool's capabilities and are happy to assist with integrating it with other electronic structure codes. If you have queries about ElasTool, need help using it, or wish to share suggestions for its improvement, please reach out to us. Our team is dedicated to supporting your work and enhancing ElasTool's functionality.

Feel free to contact us via email:
- [cekuma1@gmail.com](mailto:cekuma1@gmail.com)
- [zl.liu@163.com](mailto:zl.liu@163.com)

Your feedback and questions are invaluable to us, and we look forward to hearing from you.


