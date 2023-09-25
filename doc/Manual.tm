<TeXmacs|2.1>

<style|<tuple|generic|old-dots|british|pagella-font|captions-above|old-lengths>>

<\body>
  <doc-data|<doc-title|<with|font-shape|small-caps|ElasTool> Users
  Manual>|<doc-author|<author-data|<author-name|>>>|<doc-author|<author-data|<author-name|An
  automated toolkit for elastic constants
  calculation>>>|<doc-author|<author-data|<\author-affiliation>
    \;
  </author-affiliation>|<\author-affiliation>
    VERSION: 1.1.0

    \;

    Zhong-Li Liu and Chinedu Ekuma

    E-mail: zl.liu@163.com/che218@lehigh.edu

    Copyright (2024) Luoyang Normal University & Lehigh University

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

    <vspace*|1fn><with|font-series|bold|math-font-series|bold|2<space|2spc>Installation
    and run> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-2><vspace|0.5fn>

    <with|par-left|1tab|2.1<space|2spc>The third-party libraries
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-3>>

    <with|par-left|1tab|2.2<space|2spc>Installation
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-4>>

    <with|par-left|1tab|2.3<space|2spc>Run
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-5>>

    <with|par-left|2tab|2.3.1<space|2spc>Interface to VASP
    \ <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-6>>

    <with|par-left|1tab|2.4<space|2spc>Input files
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-7>>

    <vspace*|1fn><with|font-series|bold|math-font-series|bold|3<space|2spc>Input
    parameters> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-8><vspace|0.5fn>

    <vspace*|1fn><with|font-series|bold|math-font-series|bold|4<space|2spc>Example>
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-10><vspace|0.5fn>

    <with|par-left|1tab|4.1<space|2spc>Zero-temperature elastic constants
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-11>>

    <with|par-left|1tab|4.2<space|2spc>High-temperataure elastic constants
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-12>>

    <with|par-left|2tab|4.2.1<space|2spc>3D case
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-13>>

    <with|par-left|2tab|4.2.2<space|2spc>2D case
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-14>>

    <with|par-left|1tab|4.3<space|2spc>Output files
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-15>>

    <vspace*|1fn><with|font-series|bold|math-font-series|bold|5<space|2spc>How
    to cite> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-16><vspace|0.5fn>

    <vspace*|1fn><with|font-series|bold|math-font-series|bold|6<space|2spc>Acknowledgment>
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-17><vspace|0.5fn>

    <vspace*|1fn><with|font-series|bold|math-font-series|bold|Bibliography>
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-18><vspace|0.5fn>
  </table-of-contents>

  <new-page>

  \;

  <section|Introduction>

  <space|1em><with|font-shape|small-caps|ElasTool> is an automated toolkit
  for calculating the second-order elastic constants (SOECs) of any crystal
  systems belonging to two- and three-dimensional. It can utilize three kinds
  of strain-matrix sets, the high efficiency strain-matrix sets (OHESS)
  <cite|Liu2020>, the universal linear-independent coupling strains (ULICS)
  <cite|Yu2010> and the all single-element strain-matrix sets (ASESS)
  \ <cite|Liu2020> to calculate the SOECs automatically. In an automatic
  manner, <with|font-shape|small-caps|ElasTool> can deal with both zero- and
  high-temperature elastic constants.

  <space|1em>Presently, <with|font-shape|small-caps|ElasTool> interfaces to
  VASP package for calculating the accurate stresses of strained cystal. But
  the interfaces to other DFT packages can also be easily implemented.

  <section|Installation and run>

  <space|1em><with|font-shape|small-caps|ElasTool> \ is based on Python and
  its installation is very easy. The necessary libraries can be installed
  automatically by only one command.

  <subsection|The third-party libraries>

  <space|1em><with|font-shape|small-caps|ElasTool> depends on Python 3
  (version 3.5 or later). The following libraries are required:\ 

  \<bullet\> NumPy

  \<bullet\> Spglib

  \<bullet\> ASE

  \<bullet\> Pandas

  <subsection|Installation>

  <space|1em>In the Python 3 enviroment, the install of
  <with|font-shape|small-caps|ElasTool> is very easy. One only needs to
  execute: <strong|python3 setup.py install --prefix=/path/to/install>. Then
  execute: <strong|export PATH=/path/to/install/bin:$PATH>, or write it in
  the ~/.bashrc file.

  <subsection|Run>\ 

  <space|1em>To run <with|font-shape|small-caps|ElasTool>, one only needs to
  execute <strong|elastool> in the working directory.
  <with|font-shape|small-caps|ElasTool> \ will automatically prepare
  necessary files for calculating the stresses of crystal under deformation,
  and then call VASP to optimize initial crystal structure and calculate the
  stresses for each deformaton defined by OHESS, ASESS, or ULICS. Finally, it
  analyzes stress-strain relationship according Hooke's law and calculate all
  the elastic constants.

  <subsubsection|Interface to VASP >

  <space|1em>To run <with|font-shape|small-caps|ElasTool> with VASP, one need
  provide <strong|elastic.in>, <strong|INCARs>, and <strong|POTCAR-XX> files
  for each kind of atom, where <strong|XX> is the short name of each atom.
  For example, you can use <strong|POTCAR-Mg>, and <strong|POTCAR-O> for MgO.

  <subsection|Input files>

  <space|1em><with|font-shape|small-caps|ElasTool> needs one main input file
  for setting the calculation details of elastic constants. The main input
  file is named <strong|elatool.in>. The crystal structure file is provided
  either in <strong|POSCAR> or cif format for the structure information of
  the crystal. For VASP stress tensor calculations, <strong|INCARs>,
  <strong|KPOINTS-static>, <strong|KPOINTS-dynamic>, and <strong|POTCAR-XX>
  files are also necessary. In the <strong|INCARs> file, several INCAR files
  of VASP are collected for optimization, the static calculation, and the
  molecular dynamics simulations for high-termperature elastic constants. The
  <strong|KPOINTS-static> file is for the structure optimization and the
  static calculation of stress tensors. The <strong|KPOINTS-dynamic> file is
  for the calculation of high-temperature elastic constants using molecular
  dynamics. <strong|XX> in the <strong|POTCAR> file name is the abbreviated
  name for an element of the crystal.

  <section|Input parameters><label|ip>

  <space|1em>There are totally 10 controlling parameters for
  <with|font-shape|small-caps|ElasTool>, as listed in
  Table.<reference|parameters>. The <strong|run_mode> sets the running mode
  for ElasTool, 1 for automatic run, 2 for pre-processing, and 3 for
  post-processing. If <strong|run_mode> = 2 or 3, one should ensure the
  structure has already been optimizated at fixed pressure or volume,
  <em|i.e.> both the <strong|CONTCAR> and <strong|OUTCAR> files are in
  .<strong|/OPT> directory. In running mode 2,
  <with|font-shape|small-caps|ElasTool> will directly prepare all the
  necessary files for calculating stress tensors. After all the stress
  tensors calculations are finished, run mode 3 can analyze the output files
  and extract stress tensors, and then fit the first-order function to the
  stress-tensor data to obtain elastic constants.

  <space|1em>The <strong|dimensional> defines the dimensional of the system,
  1D, 2D, or 3D. If the system is 2D, <with|font-shape|small-caps|ElasTool>
  supposes the layered sheet in the <em|xy>-plane.

  <space|1em>The <strong|structure_file> specifies the original crystal
  structure file in POSCAR (.vasp) or cif (.cif) format.

  <space|1em>The <strong|if_conventional_cell> determines the usage of
  primitive cell (no) or conventional cell (yes).

  <space|1em>The <strong|method_stress_statistics> chooses the elastic
  constants calculation method, static or dynamic, <strong|static> for 0 K
  elastic constants, and <strong|dynamic> for high-temperature. The static
  method uses the static stresses to compute elastic constants, while the
  dynamic method deduces elastic constants from the thermal stresses obtained
  by molecular dynamics simulations.

  <space|1em>The <strong|strains_matrix> defines the type of strain-matrix
  set, OHESS, ASESS, or ULICS.

  <space|1em>The <strong|strains_list> gives one or more strains for
  calculating stresses via the strain-matrix set of OHESS, ASESS, or ULICS.

  <space|1em>The <strong|repeat_num> controls how to build a supercell from
  the primitive or conventional cell defined by <strong|if_conventional_cell>
  for the dynamic method.

  <space|1em>The <strong|num_last_samples> is the number of last MD steps to
  average thermal stresses

  <space|1em>The <strong|parallel_submit_command> is the parallel submitting
  command of <em|ab initio> code, e.g. VASP.

  <with|par-columns|1|<\big-table|<tabular|<tformat|<cwith|1|1|1|-1|cell-tborder|1ln>|<cwith|1|1|1|-1|cell-bborder|1ln>|<cwith|2|2|1|-1|cell-tborder|1ln>|<cwith|1|1|1|1|cell-lborder|0ln>|<cwith|11|11|1|-1|cell-tborder|0ln>|<cwith|10|10|1|-1|cell-bborder|0ln>|<cwith|11|11|1|-1|cell-bborder|1ln>|<cwith|11|11|1|1|cell-lborder|0ln>|<cwith|1|1|1|-1|cell-halign|c>|<cwith|1|-1|1|-1|font-base-size|10>|<table|<row|<cell|Parameters>|<cell|Values>>|<row|<cell|<samp|run_mode>>|<cell|<samp|1/2/3>>>|<row|<cell|<samp|dimensional>>|<cell|<samp|2D/3D>>>|<row|<cell|<samp|structure_file>>|<cell|<samp|file
  name ended with .vasp or .cif>>>|<row|<cell|<samp|if_conventional_cell>>|<cell|<samp|yes/no>>>|<row|<cell|<samp|method_stress_statistics>>|<cell|<samp|static/dynamic>>>|<row|<cell|<samp|strains_matrix>>|<cell|<samp|ohess/asess/ulics>>>|<row|<cell|<samp|strains_list>>|<cell|<samp|one
  or more numbers>>>|<row|<cell|<samp|repeat_num>>|<cell|<samp|3
  integers>>>|<row|<cell|<samp|num_last_samples>>|<cell|<samp|1
  integer>>>|<row|<cell|<samp|parallel_submit_command>>|<cell|<samp|DFT
  parallel run command>>>>>>>
    The controlling parameters and possible values of
    <with|font-shape|small-caps|ElasTool><label|parameters>
  </big-table>>

  <section|Example><label|exa><label|exa>

  <space|1em><with|font-shape|small-caps|ElasTool> can calculate
  zero-temperature and high-temperature elastic constants. The
  zero-temperature calculations can be conducted by the static stress
  computation. The high-temperature elastic constants can be derived by
  molecular dynamics simulations.

  <subsection|Zero-temperature elastic constants>

  <space|1em>We take the 0 K elastic constants calculation of diamond as the
  static example. The content of the input file <strong|elastool.in> is as
  follows.

  <\samp>
    <\with|font-base-size|10>
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
    </with>
  </samp>

  <subsection|High-temperataure elastic constants>

  <subsubsection|3D case>

  <space|1em>The high temperature elastic constants calculation of metal
  copper is the high-temperature example. We build a
  <math|3\<times\>3*\<times\>3> supercell from the conventional cell of
  face-centered-cubic of Cu and then perform long-time MD simulations defined
  in the INCAR-dynamic file. Because MD is very time consuming, there is only
  one strain of -0.06 is used. The last 500 MD steps are used to average
  thermal stresses.

  <with|font-base-size|10|<\samp>
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
  </samp>>

  <subsubsection|2D case>

  <\samp>
    <\with|font-base-size|10>
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
    </with>
  </samp>

  <subsection|Output files>

  <space|1em>The <samp|elastool.out> file is the unique output file of
  <with|font-shape|small-caps|ElasTool>. It includes the calculated elastic
  constants data, the elastic muduli, the sound velocity, the Debye
  temparture, the elastic anisotropy, and the stability analysis results of
  the crystal structure based on Born elastic creteria. The printed
  information on screen of the diamond example is as follows.

  <with|font-base-size|10|<\samp>
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
  </samp>>

  <section|How to cite>

  <space|1em>Please cite the following article when you use
  <with|font-shape|small-caps|ElasTool>:

  Z.<nbsp>L.<nbsp>Liu. <with|font-shape|small-caps|ElasTool>: An automated
  toolkit for elastic constants calculation (arxiv:2002.06535).
  <newblock>2020.<newblock>

  Z.<nbsp>L.<nbsp>Liu. <newblock>High-efficiency calculation of elastic
  constants enhanced by the optimized strain-matrix sets (arxiv:2002.00005).
  <newblock>2020.<newblock>

  <section|Acknowledgment>

  <space|1em>Sgplib is greatly appreciated for
  <with|font-shape|small-caps|ElasTool> development and support.\ 

  \;

  <\bibliography|bib|tm-plain|Refs>
    <\bib-list|2>
      <bibitem*|1><label|bib-Liu2020>Z.<nbsp>L.<nbsp>Liu.
      <newblock>High-efficiency calculation of elastic constants enhanced by
      the optimized strain-matrix sets (arxiv:2002.00005).
      <newblock>2020.<newblock>

      <bibitem*|2><label|bib-Yu2010>R.<nbsp>Yu, J.<nbsp>Zhu<localize|, and
      >H.<nbsp>Q.<nbsp>Ye. <newblock>Calculations of single-crystal elastic
      constants made simple. <newblock><with|font-shape|italic|Comput. Phys.
      Commun.>, 181:671, 2010.<newblock>
    </bib-list>
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
    <associate|auto-10|<tuple|4|5>>
    <associate|auto-11|<tuple|4.1|5>>
    <associate|auto-12|<tuple|4.2|6>>
    <associate|auto-13|<tuple|4.2.1|6>>
    <associate|auto-14|<tuple|4.2.2|6>>
    <associate|auto-15|<tuple|4.3|6>>
    <associate|auto-16|<tuple|5|8>>
    <associate|auto-17|<tuple|6|8>>
    <associate|auto-18|<tuple|6|8>>
    <associate|auto-2|<tuple|2|3>>
    <associate|auto-3|<tuple|2.1|3>>
    <associate|auto-4|<tuple|2.2|3>>
    <associate|auto-5|<tuple|2.3|3>>
    <associate|auto-6|<tuple|2.3.1|4>>
    <associate|auto-7|<tuple|2.4|4>>
    <associate|auto-8|<tuple|3|4>>
    <associate|auto-9|<tuple|1|5>>
    <associate|bib-Liu2020|<tuple|1|8>>
    <associate|bib-Yu2010|<tuple|2|8>>
    <associate|exa|<tuple|4|5>>
    <associate|ip|<tuple|3|4>>
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
      </surround>|<pageref|auto-9>>
    </associate>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|1<space|2spc>Introduction>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|2<space|2spc>Installation
      and run> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2><vspace|0.5fn>

      <with|par-left|<quote|1tab>|2.1<space|2spc>The third-party libraries
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-3>>

      <with|par-left|<quote|1tab>|2.2<space|2spc>Installation
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-4>>

      <with|par-left|<quote|1tab>|2.3<space|2spc>Run
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-5>>

      <with|par-left|<quote|2tab>|2.3.1<space|2spc>Interface to VASP
      \ <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-6>>

      <with|par-left|<quote|1tab>|2.4<space|2spc>Input files
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-7>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|3<space|2spc>Input
      parameters> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-8><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|4<space|2spc>Example>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-10><vspace|0.5fn>

      <with|par-left|<quote|1tab>|4.1<space|2spc>Zero-temperature elastic
      constants <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-11>>

      <with|par-left|<quote|1tab>|4.2<space|2spc>High-temperataure elastic
      constants <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-12>>

      <with|par-left|<quote|2tab>|4.2.1<space|2spc>3D case
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-13>>

      <with|par-left|<quote|2tab>|4.2.2<space|2spc>2D case
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-14>>

      <with|par-left|<quote|1tab>|4.3<space|2spc>Output files
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-15>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|5<space|2spc>How
      to cite> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-16><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|6<space|2spc>Acknowledgment>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-17><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|Bibliography>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-18><vspace|0.5fn>
    </associate>
  </collection>
</auxiliary>
