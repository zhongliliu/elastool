<TeXmacs|1.99.13>

<style|<tuple|generic|old-dots|british|pagella-font|captions-above>>

<\body>
  <doc-data|<doc-title|<with|font-shape|small-caps|ElasTool> Users
  Manual>|<doc-author|<author-data|<author-name|>>>|<doc-author|<author-data|<author-name|An
  automated toolkit for elastic constants
  calculation>>>|<doc-author|<\author-data>
    \;
  </author-data|<\author-affiliation>
    \;
  </author-affiliation>|<\author-affiliation>
    VERSION: 1.0.2

    \;

    Zhong-Li Liu\ 

    E-mail: zl.liu@163.com\ 

    Copyright (2020) Luoyang Normal University.\ 

    This software is distributed under the GNU General Public License.

    \;

    \;

    \;

    \;

    \;

    \;
  </author-affiliation>>>>

  <new-page>

  <\table-of-contents|toc>
    <vspace*|1fn><with|font-series|bold|math-font-series|bold|1<space|2spc>Introduction>
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-1><vspace|0.5fn>

    <with|par-left|1tab|1.1<space|2spc>About
    <with|font-shape|small-caps|ElasTool>
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-2>>

    <with|par-left|1tab|1.2<space|2spc>Features
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-3>>

    <vspace*|1fn><with|font-series|bold|math-font-series|bold|2<space|2spc>Installation
    and run> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-4><vspace|0.5fn>

    <with|par-left|1tab|2.1<space|2spc>Installation requirements
    \ <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-5>>

    <with|par-left|1tab|2.2<space|2spc>Installation
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-6>>

    <with|par-left|1tab|2.3<space|2spc>Run
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-7>>

    <with|par-left|2tab|2.3.1<space|2spc>Interface to VASP
    \ <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-8>>

    <with|par-left|1tab|2.4<space|2spc>Input files
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-9>>

    <vspace*|1fn><with|font-series|bold|math-font-series|bold|3<space|2spc>Input
    parameters> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-10><vspace|0.5fn>

    <vspace*|1fn><with|font-series|bold|math-font-series|bold|4<space|2spc>Example
    of run> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-12><vspace|0.5fn>

    <with|par-left|1tab|4.1<space|2spc>Zero-temperature elastic constants
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-13>>

    <with|par-left|1tab|4.2<space|2spc>High-temperataure elastic constants
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-14>>

    <with|par-left|2tab|4.2.1<space|2spc>3D case
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-15>>

    <with|par-left|2tab|4.2.2<space|2spc>2D case
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-16>>

    <with|par-left|1tab|4.3<space|2spc>Output files
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-17>>

    <vspace*|1fn><with|font-series|bold|math-font-series|bold|5<space|2spc>How
    to cite> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-18><vspace|0.5fn>

    <vspace*|1fn><with|font-series|bold|math-font-series|bold|6<space|2spc>Acknowledgment>
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-19><vspace|0.5fn>

    <vspace*|1fn><with|font-series|bold|math-font-series|bold|Bibliography>
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-20><vspace|0.5fn>
  </table-of-contents>

  <new-page>

  \;

  <section|Introduction>

  <subsection|About <with|font-shape|small-caps|ElasTool>>\ 

  <with|font-shape|small-caps|ElasTool> is an automated toolkit for
  calculating the second-order elastic constants (SOECs) of any crystal
  systems belonging to two- and three-dimensional. It can utilize three kinds
  of strain-matrix sets, the high efficiency strain-matrix sets (OHESS)
  <cite|Liu2020>, the universal linear-independent coupling strains (ULICS)
  <cite|Yu2010> and the all single-element strain-matrix sets (ASESS)
  \ <cite|Liu2020> to calculate the SOECs automatically. In an automatic
  manner, <with|font-shape|small-caps|ElasTool> can deal with both zero- and
  high-temperature elastic constants.

  Presently, <with|font-shape|small-caps|ElasTool> interfaces to VASP package
  for calculating the accurate stresses of strained cystal. But the
  interfaces to other DFT packages can also be easily implemented.

  <subsection|Features>

  <with|font-shape|small-caps|ElasTool> \ has many features including

  <\itemize-dot>
    <item>Very easy to use (installation and run);

    <item>High efficiency;

    <item>Automated flow of the SOECs calculation;

    <item>Three kinds of strain-matrix sets, the OHESS, the ASESS, and the
    ULICS;

    <item>Zero-temperature SOECs;

    <item>High-temperature SOECs;

    <item>Interfaced to VASP.
  </itemize-dot>

  <section|Installation and run>

  <with|font-shape|small-caps|ElasTool> \ is based on Python and its
  installation is very easy. But before the installation of
  <with|font-shape|small-caps|ElasTool>, the necessary libraries should be
  installed first. The following packages are required:

  <subsection|Installation requirements >

  <with|font-shape|small-caps|ElasTool> depends on Python 3. The following
  packages are required:\ 

  Python3.5 or later.

  \<bullet\> NumPy.

  \<bullet\> Spglib.

  \<bullet\> ASE.

  \<bullet\> Pandas.

  \<bullet\> VASP.

  As for Spglib, please install spglib for the structrues' space group
  determination.\ 

  <subsection|Installation>

  In the Python 3 enviroment, the necessary libraries can be installed via
  <samp|pip> command, e.g., <samp|pip install numpy>. For the construction of
  the environment and libraries, the <verbatim|miniconda3> management
  platform of Python packages is highly recommended. After installing
  <verbatim|miniconda3>, the basic Python 3 language environment is
  cocnstructed and the other libraries can be installed via either
  <samp|conda> or <samp|pip> commands. For example, one can install numpy via
  <samp|conda install -c conda-forge numpy>.

  <subsection|Run>\ 

  To run <with|font-shape|small-caps|ElasTool>, one only needs to execute
  <samp|elastool> in the working directory.
  <with|font-shape|small-caps|ElasTool> \ will automatically prepare
  necessary files for calculating the stresses of crystal under deformation,
  and then call VASP to optimize initial crystal structure and calculate the
  stresses for each deformaton defined by OHESS, ASESS, or ULICS. Finally,
  <verbatim|ElasTool> analyzes stress-strain relationship according Hooke's
  law and calculate all the elastic constants.

  <subsubsection|Interface to VASP >

  To run <with|font-shape|small-caps|ElasTool> with VASP, one need provide
  <samp|<verbatim|elastic.in>>, <verbatim|INCARs>, and <verbatim|POTCAR-XX>
  files for each kind of atom, where <samp|<verbatim|XX>> is the short name
  of each atom. For example, you can use <verbatim|POTCAR-Mg>, and
  <verbatim|POTCAR-O> for MgO.

  <subsection|Input files>

  <with|font-shape|small-caps|ElasTool> needs one main input file for setting
  the calculation details of elastic constants. The main input file is named
  <samp|elatool.in>. The crystal structure file is provided either in POSCAR
  or cif format for reading in the structure information of the crystal. For
  VASP stress tensor calculations, INCARs, KPOINTS-static, KPOINTS-dynamic,
  and POTCAR-XX files are also necessary. In the INCARs file, several INCAR
  files of VASP are collected for optimization, the static calculation, and
  the molecular dynamics simulations for high-termperature elastic constants.
  The KPOINTS-static file is for the structure optimization and the static
  calculation of stress tensors. The KPOINTS-dynamic file is for the
  calculation of high-temperature elastic constants using molecular dynamics.
  XX in the POTCAR file name is the abbreviated name for an element of the
  crystal.

  <section|Input parameters><label|ip>

  There are totally 10 controlling parameters of
  <with|font-shape|small-caps|ElasTool>, as listed in
  Table.<reference|parameters>. The <samp|run_mode> sets the running mode for
  ElasTool, 1 for automatic run, 2 for pre-processing, and 3 for
  post-processing. If <samp|run_mode = 2 or 3>, one should ensure the
  structure has already been optimizated at fixed pressure or volume,
  <em|i.e.> both the CONTCAR and OUTCAR files are in <samp|./OPT> directory.
  In running mode 2, <with|font-shape|small-caps|ElasTool> will directly
  prepare all the necessary files for calculating stress tensors. After all
  the stress tensors calculations are finished, run mode 3 can analyze the
  output files and extract stress tensors, and then fit the first-order
  function to the stress-tensor data to obtain elastic constants. The
  <samp|dimensional> defines the dimensional of the system, 2D or 3D. If the
  system is 2D, <with|font-shape|small-caps|ElasTool> supposes the layered
  sheet in the <em|xy>-plane. The <samp|structure_file> specifies the
  original crystal structure file in POSCAR (.vasp) or cif (.cif) format. The
  <samp|if_conventional_cell> determines the usage of primitive cell (no) or
  conventional cell (yes). The <samp|method_stress_statistics> chooses the
  elastic constants calculation method, static or dynamic, <samp|static> for
  0 K elastic constants, <samp|dynamic> for high-temperature. The static
  method uses the static stress to compute elastic constants, while the
  dynamic method deduces elastic constants from the thermal stresses obtained
  by molecular dynamics simulations. The <samp|strains_matrix> defines the
  type of strain-matrix set, OHESS, ASESS, or ULICS. The <samp|strains_list>
  gives one or more strains for calculating stresses via the strain-matrix
  set of OHESS, ASESS, or ULICS. The <samp|repeat_num> controls how to build
  a supercell from the primitive or conventional cell defined by
  <samp|if_conventional_cell> for the dynamic method. The
  <samp|num_last_samples >is the number of last MD steps to average thermal
  stresses. The <samp|parallel_submit_command> is the parallel submitting
  command of <em|ab initio> code, e.g. VASP.

  <with|par-columns|1|<\big-table|<tabular|<tformat|<cwith|1|1|1|-1|cell-tborder|1ln>|<cwith|1|1|1|-1|cell-bborder|1ln>|<cwith|2|2|1|-1|cell-tborder|1ln>|<cwith|1|1|1|1|cell-lborder|0ln>|<cwith|11|11|1|-1|cell-tborder|0ln>|<cwith|10|10|1|-1|cell-bborder|0ln>|<cwith|11|11|1|-1|cell-bborder|1ln>|<cwith|11|11|1|1|cell-lborder|0ln>|<cwith|1|1|1|-1|cell-halign|c>|<table|<row|<cell|Parameters>|<cell|Values>>|<row|<cell|<samp|run_mode>>|<cell|<samp|1/2/3>>>|<row|<cell|<samp|dimensional>>|<cell|<samp|2D/3D>>>|<row|<cell|<samp|structure_file>>|<cell|<samp|file
  name ended with .vasp or .cif>>>|<row|<cell|<samp|if_conventional_cell>>|<cell|<samp|yes/no>>>|<row|<cell|<samp|method_stress_statistics>>|<cell|<samp|static/dynamic>>>|<row|<cell|<samp|strains_matrix>>|<cell|<samp|ohess/asess/ulics>>>|<row|<cell|<samp|strains_list>>|<cell|<samp|one
  or more numbers>>>|<row|<cell|<samp|repeat_num>>|<cell|<samp|3
  integers>>>|<row|<cell|<samp|num_last_samples>>|<cell|<samp|1
  integer>>>|<row|<cell|<samp|parallel_submit_command>>|<cell|<samp|DFT
  parallel run command>>>>>>>
    The controlling parameters and possible values of
    <with|font-shape|small-caps|ElasTool><label|parameters>
  </big-table>>

  <section|Example of run><label|exa><label|exa>

  The best way to learn <with|font-shape|small-caps|ElasTool> is to start
  from the examples. <with|font-shape|small-caps|ElasTool> can calculate
  zero-temperature and high-temperature elastic constants. The
  zero-temperature calculations can be conducted by the static stress
  computation. The high-temperature elastic constants can be derived by
  molecular dynamics simulations.

  <subsection|Zero-temperature elastic constants>

  We take the 0 K elastic constants calculation of diamond as the static
  example. The content of the input file <samp|elastool.in> is as follows.

  <\samp>
    run_mode = 1

    dimensional = 3D

    structure_file = diamond.cif

    if_conventional_cell = no

    method_stress_statistics = static

    strains_matrix = ohess

    strains_list = -0.06 -0.03 0.03 0.06

    #repeat_num = 1 1 1

    #num_last_samples = 1

    parallel_submit_command = mpirun -np 28 vasp544
  </samp>

  <subsection|High-temperataure elastic constants>

  <subsubsection|3D case>

  The high temperature elastic constants calculation of metal copper is the
  high-temperature example. We build a <math|3\<times\>3*\<times\>3>
  supercell from the conventional cell of face-centered-cubic of Cu and then
  perform long-time MD simulations defined in the INCAR-dynamic file. Because
  MD is very time consuming, there is only one strain of -0.06 is used. The
  last 500 MD steps are used to average thermal stresses.

  <\samp>
    run_mode = 1

    dimensional = 3D

    structure_file = CONTCAR.vasp

    if_conventional_cell = yes

    method_stress_statistics = dynamic

    strains_matrix = ohess

    strains_list = -0.06

    repeat_num = 3 3 3

    num_last_samples = 500

    parallel_submit_command = mpirun -np 28 vasp544
  </samp>

  <subsubsection|2D case>

  <\samp>
    run_mode = 1

    dimensional = 2D

    structure_file = CONTCAR.vasp

    if_conventional_cell = no

    method_stress_statistics = dynamic

    strains_matrix = ohess

    strains_list = -0.06

    repeat_num = 4 4 1

    num_last_samples = 500

    parallel_submit_command = mpirun -np 28 vasp544
  </samp>

  <subsection|Output files>

  The <samp|elastool.out> file is the unique output file of
  <with|font-shape|small-caps|ElasTool>. It includes the calculated elastic
  constants data, the elastic muduli, the sound velocity, the Debye
  temparture, the elastic anisotropy, and the stability analysis results of
  the crystal structure based on Born elastic creteria. The printed
  information on screen of the diamond example is as follows.

  <\samp>
    Reading controlling parameters from elastool.in...

    Calculating stresses using the OHESS strain matrices...

    strain = -0.060

    strain = -0.030

    strain = 0.030

    strain = 0.060

    Fitting the first-order function to the collected\ 

    stress-strain data according to Hooke's law...

    The finnal results are as follows:

    +===================+

    \|This is a 3D Cubic lattice. \ \ \ \ \ \ \ \ \ \ \ \ \ \ \|

    \|------------------------------------------------\|

    \|Mean Pressure = -0.08 GPa \ \ \ \ \ \ \ \ \ \ \ \ \|

    \|------------------------------------------------\|

    \|Elastic constants: \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \|

    \|C11 = 1055.04 GPa \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \|

    \|C12 = 136.56 GPa \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \|

    \|C44 = 567.76 GPa \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \|

    \|------------------------------------------------\|

    \|Elastic moduli: \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \|

    \|B_V = 442.72 GPa \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \|

    \|B_R = 442.72 GPa \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \|

    \|G_V = 524.36 GPa \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \|

    \|G_R = 518.73 GPa \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \|

    \|B_VRH = 442.72 GPa \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \|

    \|G_VRH = 521.54 GPa \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \|

    \|Young's modulus (E) = 1123.47 GPa \|

    \|Possion's ratio (V) = 0.0771 \ \ \ \ \ \ \ \ \ \ \ \|

    \|------------------------------------------------\|

    \|Sound velocity: \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \|

    \|V_S = 12.20 Km/s \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \|

    \|V_B = 11.24 Km/s \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \|

    \|V_P = 18.03 Km/s \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \|

    \|V_M = 13.31 Km/s \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \|

    \|------------------------------------------------\|

    \|Debye temperature: \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \|

    \|T_D = 1761.03 K \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \|

    \|------------------------------------------------\|

    \|Elastic anisotropy: \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \|

    \|A_U = 0.0542 \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \|

    \|A_C = 0.0054 \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \|

    \|------------------------------------------------\|

    \|Structure stability analysis... \ \ \ \ \ \ \ \ \ \ \ \ \|

    \|This structure is mechanically STABLE.\|

    +=====================+

    Results are also saved in the elastool.out file.

    Well done! GOOD LUCK!
  </samp>

  <section|How to cite>

  Please cite the following article when you use
  <with|font-shape|small-caps|ElasTool>:

  Z.<nbsp>L.<nbsp>Liu. <with|font-shape|small-caps|ElasTool>: An automated
  toolkit for elastic constants calculation (arxiv:2002.06535).
  <newblock>2020.<newblock>

  Z.<nbsp>L.<nbsp>Liu. <newblock>High-efficiency calculation of elastic
  constants enhanced by the optimized strain-matrix sets (arxiv:2002.00005).
  <newblock>2020.<newblock>

  <section|Acknowledgment>

  Sgplib is greatly appreciated for <with|font-shape|small-caps|ElasTool>
  development and support.\ 

  \;

  <\bibliography|bib|tm-plain|Refs>
    \;
  </bibliography>
</body>

<\initial>
  <\collection>
    <associate|font|pagella>
    <associate|font-base-size|11>
    <associate|font-family|rm>
    <associate|math-font|math-pagella>
    <associate|page-medium|paper>
    <associate|page-screen-margin|false>
  </collection>
</initial>

<\references>
  <\collection>
    <associate|auto-1|<tuple|1|3>>
    <associate|auto-10|<tuple|3|5>>
    <associate|auto-11|<tuple|1|5>>
    <associate|auto-12|<tuple|4|5>>
    <associate|auto-13|<tuple|4.1|6>>
    <associate|auto-14|<tuple|4.2|6>>
    <associate|auto-15|<tuple|4.2.1|6>>
    <associate|auto-16|<tuple|4.2.2|6>>
    <associate|auto-17|<tuple|4.3|7>>
    <associate|auto-18|<tuple|5|8>>
    <associate|auto-19|<tuple|6|9>>
    <associate|auto-2|<tuple|1.1|3>>
    <associate|auto-20|<tuple|6|9>>
    <associate|auto-3|<tuple|1.2|3>>
    <associate|auto-4|<tuple|2|3>>
    <associate|auto-5|<tuple|2.1|3>>
    <associate|auto-6|<tuple|2.2|4>>
    <associate|auto-7|<tuple|2.3|4>>
    <associate|auto-8|<tuple|2.3.1|4>>
    <associate|auto-9|<tuple|2.4|4>>
    <associate|exa|<tuple|4|5>>
    <associate|ip|<tuple|3|5>>
    <associate|parameters|<tuple|1|5>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|bib>
      Liu2020

      Yu2010

      Liu2020
    </associate>
    <\associate|table>
      <tuple|normal|<\surround|<hidden-binding|<tuple>|1>|>
        The controlling parameters and possible values of
        <with|font-shape|<quote|small-caps>|ElasTool>
      </surround>|<pageref|auto-11>>
    </associate>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|1<space|2spc>Introduction>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>

      <with|par-left|<quote|1tab>|1.1<space|2spc>About
      <with|font-shape|<quote|small-caps>|ElasTool>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2>>

      <with|par-left|<quote|1tab>|1.2<space|2spc>Features
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-3>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|2<space|2spc>Installation
      and run> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-4><vspace|0.5fn>

      <with|par-left|<quote|1tab>|2.1<space|2spc>Installation requirements
      \ <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-5>>

      <with|par-left|<quote|1tab>|2.2<space|2spc>Installation
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-6>>

      <with|par-left|<quote|1tab>|2.3<space|2spc>Run
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-7>>

      <with|par-left|<quote|2tab>|2.3.1<space|2spc>Interface to VASP
      \ <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-8>>

      <with|par-left|<quote|1tab>|2.4<space|2spc>Input files
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-9>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|3<space|2spc>Input
      parameters> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-10><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|4<space|2spc>Example
      of run> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-12><vspace|0.5fn>

      <with|par-left|<quote|1tab>|4.1<space|2spc>Zero-temperature elastic
      constants <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-13>>

      <with|par-left|<quote|1tab>|4.2<space|2spc>High-temperataure elastic
      constants <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-14>>

      <with|par-left|<quote|2tab>|4.2.1<space|2spc>3D case
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-15>>

      <with|par-left|<quote|2tab>|4.2.2<space|2spc>2D case
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-16>>

      <with|par-left|<quote|1tab>|4.3<space|2spc>Output files
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-17>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|5<space|2spc>How
      to cite> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-18><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|6<space|2spc>Acknowledgment>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-19><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|Bibliography>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-20><vspace|0.5fn>
    </associate>
  </collection>
</auxiliary>