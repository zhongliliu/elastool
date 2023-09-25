# ElasTool

ElasTool is an innovative Python-based toolkit specifically designed for computing the second-order elastic constants (SOECs) and mechanical properties of crystal systems of 3D, 2D, coaxially rolled 2D-based van der Waals, and 1D nanotube structures. The software utilizes three kinds of strain-matrix sets: High-Efficiency Strain-Matrix Sets (OHESS) [1], Universal Linear-Independent Coupling Strains (ULICS) [2], and All-Single-Element Strain-Matrix Sets (ASESS) [1], enabling automatic and efficient calculation of the SOECs.

## 1. About ElasTool

ElasTool offers a flexible approach to determining elastic constants and mechanical properties of various materials at zero and finite temperatures and pressures, providing a broad scope of utility across different material conditions. It seamlessly integrates with the VASP electronic structure code. However, its architecture is open and can easily implement interfaces to other DFT packages. If you are interested in such extensions, please don't hesitate to contact the authors for further guidance and support with the ElasTool source code.

### Publications
- [1] Z.-L. Liu, C.E. Ekuma, W.-Q. Li, J.-Q. Yang, and X.-J. Li, *ElasTool: An automated toolkit for elastic constants calculation*, Computer Physics Communications 270, 108180, 2022.
- [2] R. Yu, J. Zhu, and H. Q. Ye, *Calculations of single-crystal elastic constants made simple*, Comput. Phys. Commun., 181:671, 2010.
- [3] C.E. Ekuma and Z. L. Liu, *Mechanical properties and hardness of boron pnicogens BX (X = N, P, As)*, Materialia 14, 100904 (2020).
- [4] Z.-L. Liu, Y.-D. Wei, X.-D. Xu, W.-Q. Li, G. Lv, J.-Q. Yang, X.-J. Li, and C. E. Ekuma, *Investigating elastic constants across diverse strain-matrix sets*, Comput. Mater. Sci., 2023.
- [5] S.M. Kastuar, C.E. Ekuma, Z.-L. Liu, *Efficient prediction of temperature-dependent elastic and mechanical properties of 2D materials*, Nature Scientific Reports 12, 3776 (2022)

## 2. Key Features

The ElasTool toolkit offers several features, including:
- Easy to use (installation and running)
- High efficiency
- Automated flow of the SOECs calculations
- Choice of three kinds of strain-matrix sets: OHESS, ASESS, and ULICS
- Capability to compute SOECs at zero temperature
- Capability to compute SOECs at high-temperature and/or high-pressure
- Ability to calculate the elastic and mechanical properties of 3D, 2D, coaxially rolled 2D-based van der Waals, and 1D nanotube materials.
## 3. Usage

To run ElasTool, please follow any of the examples given in the example folder. If you are interested in extending ElasTool to other electronic structure codes, please email the authors for assistance with the description of ElasTool's source code.

## 4. Contact Information

If you have any questions about using ElasTool or have suggestions for improving it, please don't hesitate to contact us:
- Email: che218@lehigh.edu or zl.liu@163.com 

## 5. Citing ElasTool

If you have used the ElasTool code in your research, please consider citing the following references:
- For the main ElasTool implementation please, cite [1].
- For related work please, cite [2], [3], [4].
- For work related to 2D materials please, cite [5] .

