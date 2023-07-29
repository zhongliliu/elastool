# ElasTool: An Automated Toolkit for Elastic Constant Calculations 

**ElasTool** is an innovative Python-based toolkit specifically designed for computing the second-order elastic constants (SOECs) of crystal systems in both two- and three-dimensional structures. The software uses three kinds of strain-matrix sets: High-Efficiency Strain-Matrix Sets (OHESS) [1], Universal Linear-Independent Coupling Strains (ULICS) [2], and All-Single-Element Strain-Matrix Sets (ASESS) [1]. This variety allows for automatic and efficient calculation of the SOECs. 

ElasTool offers a flexible approach to determining elastic constants and mechanical properties of various materials at zero and finite temperatures and pressures, providing a broad scope of utility across different material conditions. 

Currently, ElasTool integrates seamlessly with the VASP electronic structure code. However, its architecture is open to further expansions and can easily implement interfaces to other DFT packages. If you are interested in such extensions, don't hesitate to contact the authors for further guidance and support with the ElasTool source code.

## Installation 

For detailed instructions on how to install ElasTool, please refer to the `INSTALL` file in the repository. 

## Usage 

To learn how to run ElasTool, kindly refer to the examples provided in the `example` folder. 

## Relevant Literature

1. For more information about ElasTool, please see the following articles:

    * [ElasTool: An automated toolkit for elastic constants calculation](https://www.sciencedirect.com/science/article/abs/pii/S0010465521002927) - Z.-L. Liu, C.E. Ekuma, W.-Q. Li, J.-Q. Yang, and X.-J. Li, Computer Physics Communications *270*, 108180, 2022.
       
    * [Calculations of single-crystal elastic constants made simple](https://arxiv.org/abs/2002.06535) - R. Yu, J. Zhu, and H. Q. Ye. Comput. Phys. Commun., 181:671, 2010.
    
    * [Mechanical properties and hardness of boron pnicogens BX (X = N, P, As)](https://doi.org/10.1016/j.mtla.2020.100904) - C.E. Ekuma and Z. L. Liu. Materialia 14, 100904 (2020). 

3. For an example of high-pressure elastic constants calculation using ElasTool, refer to the article [here](https://arxiv.org/abs/2005.04331).

## Key Features

ElasTool offers an array of robust features that makes it a comprehensive toolkit for computing SOECs:

* User-friendly installation and operation process
* High-efficiency calculations 
* Automated workflow for SOECs calculation
* A choice between three kinds of strain-matrix sets: OHESS, ASESS, and ULICS
* Capability for both zero-temperature SOECs and high-temperature and/or high-pressure SOECs calculations.
