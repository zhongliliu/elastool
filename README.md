# ElasTool

ElasTool is an innovative Python-based toolkit specifically designed for computing the second-order elastic constants (SOECs) of crystal systems in both two- and three-dimensional structures. The software uses three kinds of strain-matrix sets: 

- **High-Efficiency Strain-Matrix Sets (OHESS)** [[1]](#references)
- **Universal Linear-Independent Coupling Strains (ULICS)** [[2]](#references)
- **All-Single-Element Strain-Matrix Sets (ASESS)** [[1]](#references)

This variety of matrix sets allows for automatic and efficient calculation of the SOECs.

ElasTool offers a flexible approach to determining the elastic constants and mechanical properties of various materials at zero and finite temperatures and pressures, providing a broad scope of utility across different material conditions.

Currently, ElasTool integrates seamlessly with the VASP electronic structure code. However, its architecture is open to further expansion and can easily implement interfaces to other DFT packages. If you are interested in such extensions, don't hesitate to contact the authors for further guidance and support with the ElasTool source code.

To run Elastool, please follow any of the examples given in the example folder.

Presently, ElasTool interfaces with VASP electronic structure code. However, the interfaces to other DFT packages can also be easily implemented. If you're interested in extending ElasTool to other electronic structure codes, please email the authors if you need assistance with the description of the ElasTool source code.

---
## Key Features

The ElasTool toolkit has many features, including:

- Very easy to use (installation and running)
- High efficiency
- Automated flow of the SOECs calculations
- The choice of three kinds of strain-matrix sets: the OHESS, ASESS, and ULICS
- Zero-temperature SOECs
- High-temperature and/or high-pressure SOECs.

---

<a name="Contact Information"></a> 
## Contact Information
**Emails:** [che218@lehigh.edu](mailto:che218@lehigh.edu) or [zl.liu@163.com](mailto:zl.liu@163.com)
Please don't hesitate to contact us if you have any questions about using ElasTool or suggestions for improving ElasTool.

---
<a name="references"></a> 
## References
### If you have used the ElasTool code in your research, please consider citing the following references:

***For the main ElasTool*** implementation, please cite:
- [1] ElasTool: An automated toolkit for elastic constants calculation - Z.-L. Liu, C.E. Ekuma, W.-Q. Li, J.-Q. Yang, and X.-J. Li, Computer Physics Communications 270, 108180, 2022.
  
### Other related citations are:
- [2] Calculations of single-crystal elastic constants made simple - R. Yu, J. Zhu, and H. Q. Ye. Comput. Phys. Commun., 181:671, 2010.
- [3] Mechanical properties and hardness of boron pnicogens BX (X = N, P, As) - C.E. Ekuma and Z. L. Liu. Materialia 14, 100904 (2020).
- [4] Z. L. Liu. High-efficiency calculation of elastic constants enhanced by the optimized strain-matrix sets (arxiv:2002.00005). 2020.

### 2D materials
***If you have used Elastool*** for computing the elastic and mechanical properties of 2D materials, please consider citing the article(s):
- [5] Efficient prediction of temperature-dependent elastic and mechanical properties of 2D materials - S.M. Kastuar, C.E. Ekuma, Z.-L. Liu, Nature Scientific Reports **12**, 3776 (2022).


--- 
