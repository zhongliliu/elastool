# ElasTool

ElasTool is an innovative Python-based toolkit specifically designed for computing the second-order elastic constants (SOECs) of crystal systems in both two- and three-dimensional structures. The software uses three kinds of strain-matrix sets: 

- **High-Efficiency Strain-Matrix Sets (OHESS)** [[1]](#references)
- **Universal Linear-Independent Coupling Strains (ULICS)** [[2]](#references)
- **All-Single-Element Strain-Matrix Sets (ASESS)** [[1]](#references)

This variety allows for automatic and efficient calculation of the SOECs.

ElasTool offers a flexible approach to determining the elastic constants and mechanical properties of various materials at zero and finite temperatures and pressures, providing a broad scope of utility across different material conditions.

Currently, ElasTool integrates seamlessly with the VASP electronic structure code. However, its architecture is open to further expansions and can easily implement interfaces to other DFT packages. If you are interested in such extensions, don't hesitate to contact the authors for further guidance and support with the ElasTool source code.

To run Elastool, please follow any of the examples given in the example folder.

Presently, ElasTool interfaces to VASP electronic structure code. But the interfaces to other DFT packages can also be easily implemented. If you're interested in extending ElasTool to other electronic structure codes, please email the authors if you need assistance on the description of ElasTool source code.

## 1. About ElasTool

- [1] ElasTool: An automated toolkit for elastic constants calculation - Z.-L. Liu, C.E. Ekuma, W.-Q. Li, J.-Q. Yang, and X.-J. Li, Computer Physics Communications 270, 108180, 2022.
- [2] Calculations of single-crystal elastic constants made simple - R. Yu, J. Zhu, and H. Q. Ye. Comput. Phys. Commun., 181:671, 2010.
- [3] Mechanical properties and hardness of boron pnicogens BX (X = N, P, As) - C.E. Ekuma and Z. L. Liu. Materialia 14, 100904 (2020).
- [4] Z. L. Liu. High-efficiency calculation of elastic constants enhanced by the optimized strain-matrix sets (arxiv:2002.00005). 2020.

## 2. Key Features

The ElasTool toolkit has many features including:

- Very easy to use (installation and run)
- High efficiency
- Automated flow of the SOECs calculation
- The choice of three kinds of strain-matrix sets: the OHESS, ASESS, and ULICS
- Zero-temperature SOECs
- High-temperature and/or high-pressure SOECs.

---

<a name="references"></a> 
### References

For quick reference, here are the associated papers:

- [1] ElasTool: An automated toolkit for elastic constants calculation
- [2] Calculations of single-crystal elastic constants made simple
- [3] Mechanical properties and hardness of boron pnicogens
- [4] High-efficiency calculation of elastic constants

--- 
