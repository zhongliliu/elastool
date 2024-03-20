"""
  Elastool -- Elastic toolkit for zero and finite-temperature elastic constants and mechanical properties calculations

  Copyright (C) 2019-2024 by Chinedu Ekuma

  This program is free software; you can redistribute it and/or modify it under the
  terms of the GNU General Public License as published by the Free Software Foundation
  version 3 of the License.

  This program is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
  PARTICULAR PURPOSE.  See the GNU General Public License for more details.

  E-mail: cekuma1@gmail.com

"""

import os
import shutil
# Directory to store the Gnuplot scripts
#script_dir = "plot_christoffel"

# Create the directory if it doesn't exist
#os.makedirs(script_dir, exist_ok=True)

# Gnuplot script contents
phase_cube_script = """
#########################################################################
#                                                                       #
#   Elastool -- Elastic toolkit for zero and finite-temperature         #
#   elastic constants and mechanical properties calculations            #
#                                                                       #
#   Modified from the Christofel software                               #
#   Copyright (C) 2019-2024 by Chinedu Ekuma                            #
#                                                                       #
#   This program is free software; you can redistribute it and/or       #
#   modify it under the terms of the GNU General Public License         #
#   as published by the Free Software Foundation, version 3 of          #
#   the License.                                                        #
#                                                                       #
#   This program is distributed in the hope that it will be useful,     #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of      #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       #
#   GNU General Public License for more details.                        #
#                                                                       #
#   E-mail: cekuma1@gmail.com                                           #
#                                                                       #
#########################################################################
reset
set terminal pngcairo truecolor enhanced font 'Arial, 22' size 2250,750
set output "phase_velocity_cube.png"

set palette defined (-1 "blue", 0 "white", 1 "red")
set pm3d depthorder explicit
set view equal xyz
set xyplane 0
unset border
unset xtics
unset ytics
unset ztics
set cblabel "v_p (km/s)"
set cbtics 0.5  # Adjust the interval as needed

set view 60,120,0.95
set arrow from first -1,-1,+1 to first -1,-1,1.3 front
set arrow from first -1,+1,-1 to first -1,1.4,-1 front
set arrow from first +1,-1,-1 to first 1.5,-1,-1 front

set label 'X' at first 1.45,-0.8,-1 front 
set label 'Y' at first -0.65,1.3,-1 front
set label 'Z' at first -1.1,-0.9,1.2 front

set arrow from +1,+1,+1 to -1,+1,+1 nohead front
set arrow from +1,+1,+1 to +1,-1,+1 nohead front
set arrow from +1,+1,+1 to +1,+1,-1 nohead front

set arrow from -1,+1,+1 to -1,+1,-1 nohead front
set arrow from -1,+1,+1 to -1,-1,+1 nohead front

set arrow from +1,-1,+1 to -1,-1,+1 nohead front
set arrow from +1,-1,+1 to +1,-1,-1 nohead front

set arrow from +1,+1,-1 to -1,+1,-1 nohead front
set arrow from +1,+1,-1 to +1,-1,-1 nohead front

set multiplot layout 1,3
set title "Slow Secondary"
splot "slow_secondary.dat" u 3:4:5:6 w pm3d notitle, "slow_secondary.dat" u (-$3):(-$4):(-$5):6 w pm3d notitle;
set title "Fast Secondary"
splot "fast_secondary.dat" u 3:4:5:6 w pm3d notitle, "fast_secondary.dat" u (-$3):(-$4):(-$5):6 w pm3d notitle;
set title "Primary"
splot "primary.dat" u 3:4:5:6 w pm3d notitle, "primary.dat" u (-$3):(-$4):(-$5):6 w pm3d notitle;
unset multiplot
unset output

reset
"""

phase_eqar_script = """
#########################################################################
#                                                                       #
#   Elastool -- Elastic toolkit for zero and finite-temperature         #
#   elastic constants and mechanical properties calculations            #
#                                                                       #
#   Modified from the Christofel software                               #
#   Copyright (C) 2019-2024 by Chinedu Ekuma                            #
#                                                                       #
#   This program is free software; you can redistribute it and/or       #
#   modify it under the terms of the GNU General Public License         #
#   as published by the Free Software Foundation, version 3 of          #
#   the License.                                                        #
#                                                                       #
#   This program is distributed in the hope that it will be useful,     #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of      #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       #
#   GNU General Public License for more details.                        #
#                                                                       #
#   E-mail: cekuma1@gmail.com                                           #
#                                                                       #
#########################################################################
reset
set terminal pngcairo truecolor enhanced font 'Arial, 22' size 2250,750
set output "phase_velocity_eqar.png"

set palette defined (-1 "blue", 0 "white", 1 "red")
set view equal xyz
set xyplane 0
unset border
unset xtics
unset ytics
unset ztics
set cblabel "v_p (km/s)"
#set cbtics 0.5  # Adjust the interval as needed

set arrow from first 1,0 to first 1.2,0 back
set arrow from first 0,1 to first 0,1.25 back

set label 'X' at first 1.1,0.1 front 
set label 'Y' at first 0.1,1.15 front

set view 0,0,1.1

set multiplot layout 1,3
set title "Slow Secondary"
splot "slow_secondary.dat" u (sqrt(2)*sin(0.5*$1)*cos($2)):(sqrt(2)*sin(0.5*$1)*sin($2)):6:6 w pm3d notitle;
set title "Fast Secondary"
splot "fast_secondary.dat" u (sqrt(2)*sin(0.5*$1)*cos($2)):(sqrt(2)*sin(0.5*$1)*sin($2)):6:6 w pm3d notitle;
set title "Primary"
splot "primary.dat" u (sqrt(2)*sin(0.5*$1)*cos($2)):(sqrt(2)*sin(0.5*$1)*sin($2)):6:6 w pm3d notitle;
unset multiplot
unset output

reset
"""


phase_radius = """
#########################################################################
#                                                                       #
#   Elastool -- Elastic toolkit for zero and finite-temperature         #
#   elastic constants and mechanical properties calculations            #
#                                                                       #
#   Modified from the Christofel software                               #
#   Copyright (C) 2019-2024 by Chinedu Ekuma                            #
#                                                                       #
#   This program is free software; you can redistribute it and/or       #
#   modify it under the terms of the GNU General Public License         #
#   as published by the Free Software Foundation, version 3 of          #
#   the License.                                                        #
#                                                                       #
#   This program is distributed in the hope that it will be useful,     #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of      #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       #
#   GNU General Public License for more details.                        #
#                                                                       #
#   E-mail: cekuma1@gmail.com                                           #
#                                                                       #
#########################################################################
reset
set terminal pngcairo truecolor enhanced font 'Arial, 22' size 2250,750
set output "phase_velocity_radius.png"

set palette defined (-1 "blue", 0 "white", 1 "red")
set pm3d depthorder explicit
set view equal xyz
set xyplane 0
unset border
unset xtics
unset ytics
unset ztics
set cblabel "v_p (km/s)"
set cbtics 1  # Adjust the interval as needed

set view 60,120,1.3

set multiplot layout 1,3
set title "Slow Secondary"
splot "slow_secondary.dat" u ($6*sin($1)*cos($2)):($6*sin($1)*sin($2)):($6*cos($1)):6 w pm3d notitle, "slow_secondary.dat" u (-$6*sin($1)*cos($2)):(-$6*sin($1)*sin($2)):(-$6*cos($1)):6 w pm3d notitle;
set title "Fast Secondary"
splot "fast_secondary.dat" u ($6*sin($1)*cos($2)):($6*sin($1)*sin($2)):($6*cos($1)):6 w pm3d notitle, "fast_secondary.dat" u (-$6*sin($1)*cos($2)):(-$6*sin($1)*sin($2)):(-$6*cos($1)):6 w pm3d notitle;
set title "Primary"
splot "primary.dat" u ($6*sin($1)*cos($2)):($6*sin($1)*sin($2)):($6*cos($1)):6 w pm3d notitle, "primary.dat" u (-$6*sin($1)*cos($2)):(-$6*sin($1)*sin($2)):(-$6*cos($1)):6 w pm3d notitle;
unset multiplot
unset output

reset
"""



phase_rel_cube = """
#########################################################################
#                                                                       #
#   Elastool -- Elastic toolkit for zero and finite-temperature         #
#   elastic constants and mechanical properties calculations            #
#                                                                       #
#   Modified from the Christofel software                               #
#   Copyright (C) 2019-2024 by Chinedu Ekuma                            #
#                                                                       #
#   This program is free software; you can redistribute it and/or       #
#   modify it under the terms of the GNU General Public License         #
#   as published by the Free Software Foundation, version 3 of          #
#   the License.                                                        #
#                                                                       #
#   This program is distributed in the hope that it will be useful,     #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of      #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       #
#   GNU General Public License for more details.                        #
#                                                                       #
#   E-mail: cekuma1@gmail.com                                           #
#                                                                       #
#########################################################################
reset
set terminal pngcairo truecolor enhanced font 'Arial, 22' size 2250,750
set output "phase_velocity_relative_cube.png"

# Get min and max values from the data files
stats "slow_secondary.dat" u 7 nooutput
S_min = floor(STATS_min)
S_max = ceil(STATS_max)

stats "fast_secondary.dat" u 7 nooutput
if (floor(STATS_min) < S_min) S_min = floor(STATS_min)
if (ceil(STATS_max) > S_max) S_max = ceil(STATS_max)

stats "primary.dat" u 7 nooutput
P_min = floor(STATS_min)
P_max = ceil(STATS_max)

if (-S_min > S_max) S_max = -S_min
if (-P_min > P_max) P_max = -P_min

set palette defined (-1 "blue", 0 "white", 1 "red")
set pm3d depthorder explicit
set view equal xyz
set xyplane 0
unset border
unset xtics
unset ytics
unset ztics
set cblabel "v_p - v_{iso} (%)"

set view 60,120,0.95
set arrow from first -1,-1,+1 to first -1,-1,1.3 front
set arrow from first -1,+1,-1 to first -1,1.4,-1 front
set arrow from first +1,-1,-1 to first 1.5,-1,-1 front

set label 'X' at first 1.45,-0.8,-1 front 
set label 'Y' at first -0.65,1.3,-1 front
set label 'Z' at first -1.1,-0.9,1.2 front

set arrow from +1,+1,+1 to -1,+1,+1 nohead front
set arrow from +1,+1,+1 to +1,-1,+1 nohead front
set arrow from +1,+1,+1 to +1,+1,-1 nohead front

set arrow from -1,+1,+1 to -1,+1,-1 nohead front
set arrow from -1,+1,+1 to -1,-1,+1 nohead front

set arrow from +1,-1,+1 to -1,-1,+1 nohead front
set arrow from +1,-1,+1 to +1,-1,-1 nohead front

set arrow from +1,+1,-1 to -1,+1,-1 nohead front
set arrow from +1,+1,-1 to +1,-1,-1 nohead front

set multiplot layout 1,3

set cbrange[-S_max:S_max]
set title "Slow Secondary"
splot "slow_secondary.dat" u 3:4:5:7 w pm3d notitle, "slow_secondary.dat" u (-$3):(-$4):(-$5):7 w pm3d notitle

set title "Fast Secondary"
splot "fast_secondary.dat" u 3:4:5:7 w pm3d notitle, "fast_secondary.dat" u (-$3):(-$4):(-$5):7 w pm3d notitle

set cbrange[-P_max:P_max]
set title "Primary"
splot "primary.dat" u 3:4:5:7 w pm3d notitle, "primary.dat" u (-$3):(-$4):(-$5):7 w pm3d notitle

unset multiplot
unset output

reset
"""


phase_rel_eqar = """
#########################################################################
#                                                                       #
#   Elastool -- Elastic toolkit for zero and finite-temperature         #
#   elastic constants and mechanical properties calculations            #
#                                                                       #
#   Modified from the Christofel software                               #
#   Copyright (C) 2019-2024 by Chinedu Ekuma                            #
#                                                                       #
#   This program is free software; you can redistribute it and/or       #
#   modify it under the terms of the GNU General Public License         #
#   as published by the Free Software Foundation, version 3 of          #
#   the License.                                                        #
#                                                                       #
#   This program is distributed in the hope that it will be useful,     #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of      #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       #
#   GNU General Public License for more details.                        #
#                                                                       #
#   E-mail: cekuma1@gmail.com                                           #
#                                                                       #
#########################################################################
reset
set terminal pngcairo truecolor enhanced font 'Arial, 22' size 2250,750
set output "phase_velocity_relative_eqar.png"

# Get min and max values from the data files
stats "slow_secondary.dat" u 7 nooutput
S_min = floor(STATS_min)
S_max = ceil(STATS_max)

stats "fast_secondary.dat" u 7 nooutput
if (floor(STATS_min) < S_min) S_min = floor(STATS_min)
if (ceil(STATS_max) > S_max) S_max = ceil(STATS_max)

stats "primary.dat" u 7 nooutput
P_min = floor(STATS_min)
P_max = ceil(STATS_max)

if (-S_min > S_max) S_max = -S_min
if (-P_min > P_max) P_max = -P_min

set palette defined (-1 "blue", 0 "white", 1 "red")
set pm3d depthorder explicit
set view equal xyz
set xyplane 0
unset border
unset xtics
unset ytics
unset ztics
set cblabel "v_p - v_{iso} (%)"

set arrow from first 1,0 to first 1.2,0 back
set arrow from first 0,1 to first 0,1.25 back

set label 'X' at first 1.1,0.1 front 
set label 'Y' at first 0.1,1.15 front

set view 0,0,1.1

set multiplot layout 1,3

set title "Slow Secondary"
set cbrange [-S_max:S_max]
splot "slow_secondary.dat" u (sqrt(2)*sin(0.5*$1)*cos($2)):(sqrt(2)*sin(0.5*$1)*sin($2)):7:7 w pm3d notitle

set title "Fast Secondary"
splot "fast_secondary.dat" u (sqrt(2)*sin(0.5*$1)*cos($2)):(sqrt(2)*sin(0.5*$1)*sin($2)):7:7 w pm3d notitle

set title "Primary"
set cbrange [-P_max:P_max]
splot "primary.dat" u (sqrt(2)*sin(0.5*$1)*cos($2)):(sqrt(2)*sin(0.5*$1)*sin($2)):7:7 w pm3d notitle

unset multiplot
unset output

reset

"""


phase_rel_sphere = """
#########################################################################
#                                                                       #
#   Elastool -- Elastic toolkit for zero and finite-temperature         #
#   elastic constants and mechanical properties calculations            #
#                                                                       #
#   Modified from the Christofel software                               #
#   Copyright (C) 2019-2024 by Chinedu Ekuma                            #
#                                                                       #
#   This program is free software; you can redistribute it and/or       #
#   modify it under the terms of the GNU General Public License         #
#   as published by the Free Software Foundation, version 3 of          #
#   the License.                                                        #
#                                                                       #
#   This program is distributed in the hope that it will be useful,     #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of      #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       #
#   GNU General Public License for more details.                        #
#                                                                       #
#   E-mail: cekuma1@gmail.com                                           #
#                                                                       #
#########################################################################
reset
set terminal pngcairo truecolor enhanced font 'Arial, 22' size 2250,750
set output "phase_velocity_relative_sphere.png"

# Get min and max values from the data files
stats "slow_secondary.dat" u 7 nooutput
S_min = floor(STATS_min)
S_max = ceil(STATS_max)

stats "fast_secondary.dat" u 7 nooutput
if (floor(STATS_min) < S_min) S_min = floor(STATS_min)
if (ceil(STATS_max) > S_max) S_max = ceil(STATS_max)

stats "primary.dat" u 7 nooutput
P_min = floor(STATS_min)
P_max = ceil(STATS_max)

if (-S_min > S_max) S_max = -S_min
if (-P_min > P_max) P_max = -P_min

set palette defined (-1 "blue", 0 "white", 1 "red")
set pm3d depthorder explicit
set view equal xyz
set xyplane 0
unset border
unset xtics
unset ytics
unset ztics
set cblabel "v_p - v_{iso} (%)"

set view 60,120,1.3
set arrow from first 0,0,1 to first 0,0,1.35 front
set arrow from first 0,1,0 to first 0,1.4,0 front
set arrow from first 1,0,0 to first 1.4,0,0 front

set label 'X' at first 1.7,0,0 front 
set label 'Y' at first 0.2,1.35,0 front
set label 'Z' at first 0,0.07,1.3 front

set multiplot layout 1,3
set title "Slow Secondary"
set cbrange [-S_max:S_max]
splot "slow_secondary.dat" u (sin($1)*cos($2)):(sin($1)*sin($2)):(cos($1)):7 w pm3d notitle, "slow_secondary.dat" u (-sin($1)*cos($2)):(-sin($1)*sin($2)):(-cos($1)):7 w pm3d notitle;
set title "Fast Secondary"
splot "fast_secondary.dat" u (sin($1)*cos($2)):(sin($1)*sin($2)):(cos($1)):7 w pm3d notitle, "fast_secondary.dat" u (-sin($1)*cos($2)):(-sin($1)*sin($2)):(-cos($1)):7 w pm3d notitle;
set title "Primary"
set cbrange [-P_max:P_max]
splot "primary.dat" u (sin($1)*cos($2)):(sin($1)*sin($2)):(cos($1)):7 w pm3d notitle, "primary.dat" u (-sin($1)*cos($2)):(-sin($1)*sin($2)):(-cos($1)):7 w pm3d notitle;
unset multiplot
unset output

reset
"""


phase_rel_stereo = """
#########################################################################
#                                                                       #
#   Elastool -- Elastic toolkit for zero and finite-temperature         #
#   elastic constants and mechanical properties calculations            #
#                                                                       #
#   Modified from the Christofel software                               #
#   Copyright (C) 2019-2024 by Chinedu Ekuma                            #
#                                                                       #
#   This program is free software; you can redistribute it and/or       #
#   modify it under the terms of the GNU General Public License         #
#   as published by the Free Software Foundation, version 3 of          #
#   the License.                                                        #
#                                                                       #
#   This program is distributed in the hope that it will be useful,     #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of      #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       #
#   GNU General Public License for more details.                        #
#                                                                       #
#   E-mail: cekuma1@gmail.com                                           #
#                                                                       #
#########################################################################
reset
set terminal pngcairo truecolor enhanced font 'Arial, 22' size 2250,750
set output "phase_velocity_relative_stereo.png"

# Get min and max values from the data files
stats "slow_secondary.dat" u 7 nooutput
S_min = floor(STATS_min)
S_max = ceil(STATS_max)

stats "fast_secondary.dat" u 7 nooutput
if (floor(STATS_min) < S_min) S_min = floor(STATS_min)
if (ceil(STATS_max) > S_max) S_max = ceil(STATS_max)

stats "primary.dat" u 7 nooutput
P_min = floor(STATS_min)
P_max = ceil(STATS_max)

if (-S_min > S_max) S_max = -S_min
if (-P_min > P_max) P_max = -P_min


set palette defined (-1 "blue", 0 "white", 1 "red")
set pm3d depthorder explicit
set view equal xyz
set xyplane 0
unset border
unset xtics
unset ytics
unset ztics
set cblabel "v_p - v_{iso} (%)"

set arrow from first 1,0 to first 1.2,0 back
set arrow from first 0,1 to first 0,1.25 back

set label 'X' at first 1.1,0.1 front 
set label 'Y' at first 0.1,1.15 front

set view 0,0,1.1

set multiplot layout 1,3
set title "Slow Secondary"
set cbrange [-S_max:S_max]
splot "slow_secondary.dat" u (tan(0.5*$1)*cos($2)):(tan(0.5*$1)*sin($2)):7:7 w pm3d notitle;
set title "Fast Secondary"
splot "fast_secondary.dat" u (tan(0.5*$1)*cos($2)):(tan(0.5*$1)*sin($2)):7:7 w pm3d notitle;
set title "Primary"
set cbrange [-P_max:P_max]
splot "primary.dat" u (tan(0.5*$1)*cos($2)):(tan(0.5*$1)*sin($2)):7:7 w pm3d notitle;
unset multiplot
unset output

reset
"""


phase_sphere = """
#########################################################################
#                                                                       #
#   Elastool -- Elastic toolkit for zero and finite-temperature         #
#   elastic constants and mechanical properties calculations            #
#                                                                       #
#   Modified from the Christofel software                               #
#   Copyright (C) 2019-2024 by Chinedu Ekuma                            #
#                                                                       #
#   This program is free software; you can redistribute it and/or       #
#   modify it under the terms of the GNU General Public License         #
#   as published by the Free Software Foundation, version 3 of          #
#   the License.                                                        #
#                                                                       #
#   This program is distributed in the hope that it will be useful,     #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of      #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       #
#   GNU General Public License for more details.                        #
#                                                                       #
#   E-mail: cekuma1@gmail.com                                           #
#                                                                       #
#########################################################################
reset
set terminal pngcairo truecolor enhanced font 'Arial, 22' size 2250,750
set output "phase_velocity_sphere.png"

set palette defined (-1 "blue", 0 "white", 1 "red")
set pm3d depthorder explicit
set view equal xyz
set xyplane 0
unset border
unset xtics
unset ytics
unset ztics
set cblabel "v_p (km/s)"
set cbtics 0.5 
set view 60,120,1.3
set arrow from first 0,0,1 to first 0,0,1.35 front
set arrow from first 0,1,0 to first 0,1.4,0 front
set arrow from first 1,0,0 to first 1.4,0,0 front

set label 'X' at first 1.7,0,0 front 
set label 'Y' at first 0.2,1.35,0 front
set label 'Z' at first 0,0.07,1.3 front

set multiplot layout 1,3
set title "Slow Secondary"
splot "slow_secondary.dat" u (sin($1)*cos($2)):(sin($1)*sin($2)):(cos($1)):6 w pm3d notitle, "slow_secondary.dat" u (-sin($1)*cos($2)):(-sin($1)*sin($2)):(-cos($1)):6 w pm3d notitle;
set title "Fast Secondary"
splot "fast_secondary.dat" u (sin($1)*cos($2)):(sin($1)*sin($2)):(cos($1)):6 w pm3d notitle, "fast_secondary.dat" u (-sin($1)*cos($2)):(-sin($1)*sin($2)):(-cos($1)):6 w pm3d notitle;
set title "Primary"
splot "primary.dat" u (sin($1)*cos($2)):(sin($1)*sin($2)):(cos($1)):6 w pm3d notitle, "primary.dat" u (-sin($1)*cos($2)):(-sin($1)*sin($2)):(-cos($1)):6 w pm3d notitle;
unset multiplot
unset output

reset
"""


phase_stereo = """
#########################################################################
#                                                                       #
#   Elastool -- Elastic toolkit for zero and finite-temperature         #
#   elastic constants and mechanical properties calculations            #
#                                                                       #
#   Modified from the Christofel software                               #
#   Copyright (C) 2019-2024 by Chinedu Ekuma                            #
#                                                                       #
#   This program is free software; you can redistribute it and/or       #
#   modify it under the terms of the GNU General Public License         #
#   as published by the Free Software Foundation, version 3 of          #
#   the License.                                                        #
#                                                                       #
#   This program is distributed in the hope that it will be useful,     #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of      #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       #
#   GNU General Public License for more details.                        #
#                                                                       #
#   E-mail: cekuma1@gmail.com                                           #
#                                                                       #
#########################################################################
reset
set terminal pngcairo truecolor enhanced font 'Arial, 22' size 2250,750
set output "phase_velocity_stereo.png"

set palette defined (-1 "blue", 0 "white", 1 "red")
set view equal xyz
set xyplane 0
unset border
unset xtics
unset ytics
unset ztics
set cblabel "v_p (km/s)"
set cbtics 0.5 
set arrow from first 1,0 to first 1.2,0 back
set arrow from first 0,1 to first 0,1.25 back

set label 'X' at first 1.1,0.1 front 
set label 'Y' at first 0.1,1.15 front

set view 0,0,1.1

set multiplot layout 1,3
set title "Slow Secondary"
splot "slow_secondary.dat" u (tan(0.5*$1)*cos($2)):(tan(0.5*$1)*sin($2)):6:6 w pm3d notitle;
set title "Fast Secondary"
splot "fast_secondary.dat" u (tan(0.5*$1)*cos($2)):(tan(0.5*$1)*sin($2)):6:6 w pm3d notitle;
set title "Primary"
splot "primary.dat" u (tan(0.5*$1)*cos($2)):(tan(0.5*$1)*sin($2)):6:6 w pm3d notitle;
unset multiplot
unset output

reset
"""


enh_cube = """
#########################################################################
#                                                                       #
#   Elastool -- Elastic toolkit for zero and finite-temperature         #
#   elastic constants and mechanical properties calculations            #
#                                                                       #
#   Modified from the Christofel software                               #
#   Copyright (C) 2019-2024 by Chinedu Ekuma                            #
#                                                                       #
#   This program is free software; you can redistribute it and/or       #
#   modify it under the terms of the GNU General Public License         #
#   as published by the Free Software Foundation, version 3 of          #
#   the License.                                                        #
#                                                                       #
#   This program is distributed in the hope that it will be useful,     #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of      #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       #
#   GNU General Public License for more details.                        #
#                                                                       #
#   E-mail: cekuma1@gmail.com                                           #
#                                                                       #
#########################################################################
reset
set terminal pngcairo truecolor enhanced font 'Arial, 22' size 2250,750
set output "enhancement_factor_cube.png"

set palette defined (-1 "blue", 0 "white", 1 "red", 2 "black")
set pm3d depthorder explicit
set view equal xyz
set xyplane 0
unset border
unset xtics
unset ytics
unset ztics
set cblabel "log_{10}(A)"

set view 60,120,0.95
set arrow from first -1,-1,+1 to first -1,-1,1.3 front
set arrow from first -1,+1,-1 to first -1,1.4,-1 front
set arrow from first +1,-1,-1 to first 1.5,-1,-1 front

set label 'X' at first 1.45,-0.8,-1 front 
set label 'Y' at first -0.65,1.3,-1 front
set label 'Z' at first -1.1,-0.9,1.2 front

set arrow from +1,+1,+1 to -1,+1,+1 nohead front
set arrow from +1,+1,+1 to +1,-1,+1 nohead front
set arrow from +1,+1,+1 to +1,+1,-1 nohead front

set arrow from -1,+1,+1 to -1,+1,-1 nohead front
set arrow from -1,+1,+1 to -1,-1,+1 nohead front

set arrow from +1,-1,+1 to -1,-1,+1 nohead front
set arrow from +1,-1,+1 to +1,-1,-1 nohead front

set arrow from +1,+1,-1 to -1,+1,-1 nohead front
set arrow from +1,+1,-1 to +1,-1,-1 nohead front

set multiplot layout 1,3

# Adjust the color bar range as needed
set cbrange [-2.5:5]
set title "Slow Secondary"
splot "slow_secondary.dat" u 3:4:5:(log10($17)) w pm3d notitle, \
      "slow_secondary.dat" u (-$3):(-$4):(-$5):(log10($17)) w pm3d notitle

set title "Fast Secondary"
splot "fast_secondary.dat" u 3:4:5:(log10($17)) w pm3d notitle, \
      "fast_secondary.dat" u (-$3):(-$4):(-$5):(log10($17)) w pm3d notitle

# Adjust the color bar range for the primary plot if different
set cbrange [-0.6:1.2]
set title "Primary"
splot "primary.dat" u 3:4:5:(log10($17)) w pm3d notitle, \
      "primary.dat" u (-$3):(-$4):(-$5):(log10($17)) w pm3d notitle

unset multiplot
unset output

reset
"""


enh_eqar = """
#########################################################################
#                                                                       #
#   Elastool -- Elastic toolkit for zero and finite-temperature         #
#   elastic constants and mechanical properties calculations            #
#                                                                       #
#   Modified from the Christofel software                               #
#   Copyright (C) 2019-2024 by Chinedu Ekuma                            #
#                                                                       #
#   This program is free software; you can redistribute it and/or       #
#   modify it under the terms of the GNU General Public License         #
#   as published by the Free Software Foundation, version 3 of          #
#   the License.                                                        #
#                                                                       #
#   This program is distributed in the hope that it will be useful,     #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of      #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       #
#   GNU General Public License for more details.                        #
#                                                                       #
#   E-mail: cekuma1@gmail.com                                           #
#                                                                       #
#########################################################################
reset
set terminal pngcairo truecolor enhanced font 'Arial, 22' size 2250,750
set output "enhancement_factor_eqar.png"

set palette defined (-1 "blue", 0 "white", 1 "red", 2 "black")
set view equal xyz
set xyplane 0
unset border
unset xtics
unset ytics
unset ztics
set cblabel "log_{10}(A)"

set arrow from first 1,0 to first 1.2,0 back
set arrow from first 0,1 to first 0,1.25 back

set label 'X' at first 1.1,0.1 front 
set label 'Y' at first 0.1,1.15 front

set view 0,0,1.1

set multiplot layout 1,3

# Adjust the color bar range as needed
set cbrange [-2.5:5]
set title "Slow Secondary"
splot "slow_secondary.dat" u (sqrt(2)*sin(0.5*$1)*cos($2)):(sqrt(2)*sin(0.5*$1)*sin($2)):(log10($17)):(log10($17)) w pm3d notitle

set title "Fast Secondary"
splot "fast_secondary.dat" u (sqrt(2)*sin(0.5*$1)*cos($2)):(sqrt(2)*sin(0.5*$1)*sin($2)):(log10($17)):(log10($17)) w pm3d notitle

# Adjust the color bar range for the primary plot if different
set cbrange [-0.6:1.2]
set title "Primary"
splot "primary.dat" u (sqrt(2)*sin(0.5*$1)*cos($2)):(sqrt(2)*sin(0.5*$1)*sin($2)):(log10($17)):(log10($17)) w pm3d notitle

unset multiplot
unset output

reset
"""


enh_sphere = """
#########################################################################
#                                                                       #
#   Elastool -- Elastic toolkit for zero and finite-temperature         #
#   elastic constants and mechanical properties calculations            #
#                                                                       #
#   Modified from the Christofel software                               #
#   Copyright (C) 2019-2024 by Chinedu Ekuma                            #
#                                                                       #
#   This program is free software; you can redistribute it and/or       #
#   modify it under the terms of the GNU General Public License         #
#   as published by the Free Software Foundation, version 3 of          #
#   the License.                                                        #
#                                                                       #
#   This program is distributed in the hope that it will be useful,     #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of      #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       #
#   GNU General Public License for more details.                        #
#                                                                       #
#   E-mail: cekuma1@gmail.com                                           #
#                                                                       #
#########################################################################
reset
set terminal pngcairo truecolor enhanced font 'Arial, 22' size 2250,750
set output "enhancement_factor_sphere.png"

set palette defined (-1 "blue", 0 "white", 1 "red", 2 "black")
set pm3d depthorder explicit
set view equal xyz
set xyplane 0
unset border
unset xtics
unset ytics
unset ztics

set cblabel "log_{10}(A)"

set view 60,120,1.3
set arrow from first 0,0,1 to first 0,0,1.35 front
set arrow from first 0,1,0 to first 0,1.4,0 front
set arrow from first 1,0,0 to first 1.4,0,0 front

set label 'X' at first 1.7,0,0 front 
set label 'Y' at first 0.2,1.35,0 front
set label 'Z' at first 0,0.07,1.3 front

set multiplot layout 1,3

# Adjust the color bar range as needed
set cbrange [-2.5:5]
set title "Slow Secondary"
splot "slow_secondary.dat" u (sin($1)*cos($2)):(sin($1)*sin($2)):(cos($1)):(log10($17)) w pm3d notitle, \
      "slow_secondary.dat" u (-sin($1)*cos($2)):(-sin($1)*sin($2)):(-cos($1)):(log10($17)) w pm3d notitle

set title "Fast Secondary"
splot "fast_secondary.dat" u (sin($1)*cos($2)):(sin($1)*sin($2)):(cos($1)):(log10($17)) w pm3d notitle, \
      "fast_secondary.dat" u (-sin($1)*cos($2)):(-sin($1)*sin($2)):(-cos($1)):(log10($17)) w pm3d notitle

# Adjust the color bar range for the primary plot if different
set cbrange [-0.6:1.2]
set title "Primary"
splot "primary.dat" u (sin($1)*cos($2)):(sin($1)*sin($2)):(cos($1)):(log10($17)) w pm3d notitle, \
      "primary.dat" u (-sin($1)*cos($2)):(-sin($1)*sin($2)):(-cos($1)):(log10($17)) w pm3d notitle

unset multiplot
unset output

reset

"""


enh_stereo = """
#########################################################################
#                                                                       #
#   Elastool -- Elastic toolkit for zero and finite-temperature         #
#   elastic constants and mechanical properties calculations            #
#                                                                       #
#   Modified from the Christofel software                               #
#   Copyright (C) 2019-2024 by Chinedu Ekuma                            #
#                                                                       #
#   This program is free software; you can redistribute it and/or       #
#   modify it under the terms of the GNU General Public License         #
#   as published by the Free Software Foundation, version 3 of          #
#   the License.                                                        #
#                                                                       #
#   This program is distributed in the hope that it will be useful,     #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of      #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       #
#   GNU General Public License for more details.                        #
#                                                                       #
#   E-mail: cekuma1@gmail.com                                           #
#                                                                       #
#########################################################################
reset
set terminal pngcairo truecolor enhanced font 'Arial, 22' size 2250,750
set output "enhancement_factor_stereo.png"

set palette defined (-1 "blue", 0 "white", 1 "red", 2 "black")
set view equal xyz
set xyplane 0
unset border
unset xtics
unset ytics
unset ztics
set cblabel "log_{10}(A)"

set arrow from first 1,0 to first 1.2,0 back
set arrow from first 0,1 to first 0,1.25 back

set label 'X' at first 1.1,0.1 front 
set label 'Y' at first 0.1,1.15 front

set view 0,0,1.1

set multiplot layout 1,3

# Adjust the color bar range as needed
set cbrange [-2.5:5]
set title "Slow Secondary"
splot "slow_secondary.dat" u (tan(0.5*$1)*cos($2)):(tan(0.5*$1)*sin($2)):(log10($17)):(log10($17)) w pm3d notitle

set title "Fast Secondary"
splot "fast_secondary.dat" u (tan(0.5*$1)*cos($2)):(tan(0.5*$1)*sin($2)):(log10($17)):(log10($17)) w pm3d notitle

# Adjust the color bar range for the primary plot if different
set cbrange [-0.6:1.2]
set title "Primary"
splot "primary.dat" u (tan(0.5*$1)*cos($2)):(tan(0.5*$1)*sin($2)):(log10($17)):(log10($17)) w pm3d notitle

unset multiplot
unset output

reset
"""


slowness = """
#########################################################################
#                                                                       #
#   Elastool -- Elastic toolkit for zero and finite-temperature         #
#   elastic constants and mechanical properties calculations            #
#                                                                       #
#   Modified from the Christofel software                               #
#   Copyright (C) 2019-2024 by Chinedu Ekuma                            #
#                                                                       #
#   This program is free software; you can redistribute it and/or       #
#   modify it under the terms of the GNU General Public License         #
#   as published by the Free Software Foundation, version 3 of          #
#   the License.                                                        #
#                                                                       #
#   This program is distributed in the hope that it will be useful,     #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of      #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       #
#   GNU General Public License for more details.                        #
#                                                                       #
#   E-mail: cekuma1@gmail.com                                           #
#                                                                       #
#########################################################################
reset
set terminal pngcairo truecolor enhanced font 'Arial, 22' size 2200,690
set output "slowness.png"

set palette defined (-1 "red", 0 "white", 1 "blue")
set pm3d depthorder explicit
set view equal xyz
set xyplane 0
set hidden3d front
unset border
unset xtics
unset ytics
unset ztics
set cblabel "Slowness (s/km)"
set format cb "%.2f"
set view 60,120,1.5

# Function to calculate cbtics interval
calculate_cbtics_interval(min, max) = abs(1.0/max - 1.0/min) / 5.0
#calculate_cbtics_interval(min, max) = abs(max - min) / 6.0

set multiplot layout 1,3
#set cbtics 0.02

# Plot 1: Slow Secondary
stats "slow_secondary.dat" u 6 nooutput
S_MIN = STATS_min
S_MAX = STATS_max
set cbtics calculate_cbtics_interval(S_MIN, S_MAX)


set title "Slow Secondary"
splot "slow_secondary.dat" u ((1.0/$6)*sin($1)*cos($2)):((1.0/$6)*sin($1)*sin($2)):((1.0/$6)*cos($1)):(1.0/$6):(1.0/$6) w p palette ps 2 pt 3 notitle, \
      "slow_secondary.dat" u (-(1.0/$6)*sin($1)*cos($2)):(-(1.0/$6)*sin($1)*sin($2)):(-(1.0/$6)*cos($1)):(1.0/$6):(1.0/$6) w p palette ps 2 pt 3 notitle

# Plot 2: Fast Secondary
stats "fast_secondary.dat" u 6 nooutput
F_MIN = STATS_min
F_MAX = STATS_max
set cbtics calculate_cbtics_interval(F_MIN, F_MAX)


set title "Fast Secondary"
splot "fast_secondary.dat" u ((1.0/$6)*sin($1)*cos($2)):((1.0/$6)*sin($1)*sin($2)):((1.0/$6)*cos($1)):(1.0/$6):(1.0/$6) w p palette ps 2 pt 3 notitle, \
      "fast_secondary.dat" u (-(1.0/$6)*sin($1)*cos($2)):(-(1.0/$6)*sin($1)*sin($2)):(-(1.0/$6)*cos($1)):(1.0/$6):(1.0/$6) w p palette ps 2 pt 3 notitle

# Plot 3: Primary
stats "primary.dat" u 6 nooutput
P_MIN = STATS_min
P_MAX = STATS_max
set cbtics calculate_cbtics_interval(P_MIN, P_MAX)

set title "Primary"
splot "primary.dat" u ((1.0/$6)*sin($1)*cos($2)):((1.0/$6)*sin($1)*sin($2)):((1.0/$6)*cos($1)):(1.0/$6):(1.0/$6) w p palette ps 2 pt 3 notitle, \
      "primary.dat" u (-(1.0/$6)*sin($1)*cos($2)):(-(1.0/$6)*sin($1)*sin($2)):(-(1.0/$6)*cos($1)):(1.0/$6):(1.0/$6) w p palette ps 2 pt 3 notitle

unset multiplot
unset output
reset

"""


slowness_groupdir = """
#########################################################################
#                                                                       #
#   Elastool -- Elastic toolkit for zero and finite-temperature         #
#   elastic constants and mechanical properties calculations            #
#                                                                       #
#   Modified from the Christofel software                               #
#   Copyright (C) 2019-2024 by Chinedu Ekuma                            #
#                                                                       #
#   This program is free software; you can redistribute it and/or       #
#   modify it under the terms of the GNU General Public License         #
#   as published by the Free Software Foundation, version 3 of          #
#   the License.                                                        #
#                                                                       #
#   This program is distributed in the hope that it will be useful,     #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of      #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       #
#   GNU General Public License for more details.                        #
#                                                                       #
#   E-mail: cekuma1@gmail.com                                           #
#                                                                       #
#########################################################################
reset
set terminal pngcairo truecolor enhanced font 'Arial, 22' size 2250,690
set output "slowness_groupdirection.png"

# Get min and max values from the data files
stats "slow_secondary.dat" u 6 nooutput
S_MIN = STATS_min
S_MAX = STATS_max

stats "fast_secondary.dat" u 6 nooutput
if (STATS_min < S_MIN) S_MIN = STATS_min
if (STATS_max > S_MAX) S_MAX = STATS_max

stats "primary.dat" u 6 nooutput
P_MIN = STATS_min
P_MAX = STATS_max

NUM_VALUES = STATS_records
NUM_THETA = 0.5*sqrt(NUM_VALUES)
SKIP = floor(NUM_THETA/12.0)

set palette defined (-1 "red", 0 "white", 1 "blue")
set pm3d depthorder explicit
set view equal xyz
set xyplane 0
set hidden3d front
unset border
unset xtics
unset ytics
unset ztics
set cblabel "Slowness (s/km)"

set view 60,120,1.3

set multiplot layout 1,3

set cbrange[1.0/S_MAX:1.0/S_MIN]

set title "Slow Secondary"
splot "slow_secondary.dat" u ((1.0/$6)*sin($1)*cos($2)):((1.0/$6)*sin($1)*sin($2)):((1.0/$6)*cos($1)):(1.0/$6):(1.0/$6) w p palette ps 2 pt 3 notitle, "slow_secondary.dat" u (-(1.0/$6)*sin($1)*cos($2)):(-(1.0/$6)*sin($1)*sin($2)):(-(1.0/$6)*cos($1)):(1.0/$6):(1.0/$6) w p palette ps 2 pt 3 notitle, "slow_secondary.dat" u ((1.0/$6)*sin($1)*cos($2)):((1.0/$6)*sin($1)*sin($2)):((1.0/$6)*cos($1)):(0.03*($13/$11)):(0.03*($14/$11)):(0.03*($15/$11)) every SKIP:SKIP w vectors lw 2 lc 2 notitle, "slow_secondary.dat" u ((1.0/$6)*sin($1)*cos($2)):((1.0/$6)*sin($1)*sin($2)):((1.0/$6)*cos($1)):(-0.03*($13/$11)):(-0.03*($14/$11)):(-0.03*($15/$11)) every SKIP:SKIP w vectors lw 2 lc 2 notitle, "slow_secondary.dat" u (-(1.0/$6)*sin($1)*cos($2)):(-(1.0/$6)*sin($1)*sin($2)):(-(1.0/$6)*cos($1)):(0.03*($13/$11)):(0.03*($14/$11)):(0.03*($15/$11)) every SKIP:SKIP w vectors lw 2 lc 2 notitle, "slow_secondary.dat" u (-(1.0/$6)*sin($1)*cos($2)):(-(1.0/$6)*sin($1)*sin($2)):(-(1.0/$6)*cos($1)):(-0.03*($13/$11)):(-0.03*($14/$11)):(-0.03*($15/$11)) every SKIP:SKIP w vectors lw 2 lc 2 notitle;

set title "Fast Secondary"
splot "fast_secondary.dat" u ((1.0/$6)*sin($1)*cos($2)):((1.0/$6)*sin($1)*sin($2)):((1.0/$6)*cos($1)):(1.0/$6):(1.0/$6) w p palette ps 2 pt 3 notitle, "fast_secondary.dat" u (-(1.0/$6)*sin($1)*cos($2)):(-(1.0/$6)*sin($1)*sin($2)):(-(1.0/$6)*cos($1)):(1.0/$6):(1.0/$6) w p palette ps 2 pt 3 notitle, "fast_secondary.dat" u ((1.0/$6)*sin($1)*cos($2)):((1.0/$6)*sin($1)*sin($2)):((1.0/$6)*cos($1)):(0.03*($13/$11)):(0.03*($14/$11)):(0.03*($15/$11)) every SKIP:SKIP w vectors lw 2 lc 2 notitle, "fast_secondary.dat" u ((1.0/$6)*sin($1)*cos($2)):((1.0/$6)*sin($1)*sin($2)):((1.0/$6)*cos($1)):(-0.03*($13/$11)):(-0.03*($14/$11)):(-0.03*($15/$11)) every SKIP:SKIP w vectors lw 2 lc 2 notitle, "fast_secondary.dat" u (-(1.0/$6)*sin($1)*cos($2)):(-(1.0/$6)*sin($1)*sin($2)):(-(1.0/$6)*cos($1)):(0.03*($13/$11)):(0.03*($14/$11)):(0.03*($15/$11)) every SKIP:SKIP w vectors lw 2 lc 2 notitle, "fast_secondary.dat" u (-(1.0/$6)*sin($1)*cos($2)):(-(1.0/$6)*sin($1)*sin($2)):(-(1.0/$6)*cos($1)):(-0.03*($13/$11)):(-0.03*($14/$11)):(-0.03*($15/$11)) every SKIP:SKIP w vectors lw 2 lc 2 notitle;

set cbrange[1.0/P_MAX:1.0/P_MIN]

set title "Primary"
splot "primary.dat" u ((1.0/$6)*sin($1)*cos($2)):((1.0/$6)*sin($1)*sin($2)):((1.0/$6)*cos($1)):(1.0/$6):(1.0/$6) w p palette ps 2 pt 3 notitle, "primary.dat" u (-(1.0/$6)*sin($1)*cos($2)):(-(1.0/$6)*sin($1)*sin($2)):(-(1.0/$6)*cos($1)):(1.0/$6):(1.0/$6) w p palette ps 2 pt 3 notitle, "primary.dat" u ((1.0/$6)*sin($1)*cos($2)):((1.0/$6)*sin($1)*sin($2)):((1.0/$6)*cos($1)):(0.03*($13/$11)):(0.03*($14/$11)):(0.03*($15/$11)) every SKIP:SKIP w vectors lw 2 lc 2 notitle, "primary.dat" u ((1.0/$6)*sin($1)*cos($2)):((1.0/$6)*sin($1)*sin($2)):((1.0/$6)*cos($1)):(-0.03*($13/$11)):(-0.03*($14/$11)):(-0.03*($15/$11)) every SKIP:SKIP w vectors lw 2 lc 2 notitle, "primary.dat" u (-(1.0/$6)*sin($1)*cos($2)):(-(1.0/$6)*sin($1)*sin($2)):(-(1.0/$6)*cos($1)):(0.03*($13/$11)):(0.03*($14/$11)):(0.03*($15/$11)) every SKIP:SKIP w vectors lw 2 lc 2 notitle, "primary.dat" u (-(1.0/$6)*sin($1)*cos($2)):(-(1.0/$6)*sin($1)*sin($2)):(-(1.0/$6)*cos($1)):(-0.03*($13/$11)):(-0.03*($14/$11)):(-0.03*($15/$11)) every SKIP:SKIP w vectors lw 2 lc 2 notitle;

unset multiplot
unset output

reset

"""


ray_equar = """
#########################################################################
#                                                                       #
#   Elastool -- Elastic toolkit for zero and finite-temperature         #
#   elastic constants and mechanical properties calculations            #
#                                                                       #
#   Modified from the Christofel software                               #
#   Copyright (C) 2019-2024 by Chinedu Ekuma                            #
#                                                                       #
#   This program is free software; you can redistribute it and/or       #
#   modify it under the terms of the GNU General Public License         #
#   as published by the Free Software Foundation, version 3 of          #
#   the License.                                                        #
#                                                                       #
#   This program is distributed in the hope that it will be useful,     #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of      #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       #
#   GNU General Public License for more details.                        #
#                                                                       #
#   E-mail: cekuma1@gmail.com                                           #
#                                                                       #
#########################################################################
reset
set terminal pngcairo truecolor enhanced font 'Arial, 22' size 2250,750
set output "ray_eqar.png"

set palette defined (-1 "blue", 0 "white", 1 "red")
set view equal xyz
set xyplane 0
unset border
unset xtics
unset ytics
unset ztics
set cblabel "v_g (km/s)"
#set cbtics 0.5
set pm3d depthorder explicit

#set view map
set view 0,0,1.3

set multiplot layout 1,3
set title "Slow Secondary"
splot "slow_secondary.dat" u (($13/$11)*sqrt(2/(1+($15/$11)))):(($14/$11)*sqrt(2/(1+($15/$11)))):11:11 w pm3d notitle;
set title "Fast Secondary"
splot "fast_secondary.dat" u (($13/$11)*sqrt(2/(1+($15/$11)))):(($14/$11)*sqrt(2/(1+($15/$11)))):11:11 w pm3d notitle;
set title "Primary"
splot "primary.dat" u (($13/$11)*sqrt(2/(1+($15/$11)))):(($14/$11)*sqrt(2/(1+($15/$11)))):11:11 w pm3d notitle;
unset multiplot
unset output

reset
"""


ray_stereo = """
#########################################################################
#                                                                       #
#   Elastool -- Elastic toolkit for zero and finite-temperature         #
#   elastic constants and mechanical properties calculations            #
#                                                                       #
#   Modified from the Christofel software                               #
#   Copyright (C) 2019-2024 by Chinedu Ekuma                            #
#                                                                       #
#   This program is free software; you can redistribute it and/or       #
#   modify it under the terms of the GNU General Public License         #
#   as published by the Free Software Foundation, version 3 of          #
#   the License.                                                        #
#                                                                       #
#   This program is distributed in the hope that it will be useful,     #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of      #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       #
#   GNU General Public License for more details.                        #
#                                                                       #
#   E-mail: cekuma1@gmail.com                                           #
#                                                                       #
#########################################################################
reset
set terminal pngcairo truecolor enhanced font 'Arial, 22' size 2250,750
set output "ray_stereo.png"

set palette defined (-1 "blue", 0 "white", 1 "red")
set view equal xyz
set xyplane 0
unset border
unset xtics
unset ytics
unset ztics
set cblabel "v_g (km/s)"
set cbtics 0.5
set pm3d depthorder explicit

#set view map
set view 0,0,1.3

set multiplot layout 1,3
set title "Slow Secondary"
splot "slow_secondary.dat" u ($13/($11+$15)):($14/($11+$15)):11:11 w pm3d notitle;
set title "Fast Secondary"
splot "fast_secondary.dat" u ($13/($11+$15)):($14/($11+$15)):11:11 w pm3d notitle;
set title "Primary"
splot "primary.dat" u ($13/($11+$15)):($14/($11+$15)):11:11 w pm3d notitle;
unset multiplot
unset output

reset
"""



ray_surface_enh = """
#########################################################################
#                                                                       #
#   Elastool -- Elastic toolkit for zero and finite-temperature         #
#   elastic constants and mechanical properties calculations            #
#                                                                       #
#   Modified from the Christofel software                               #
#   Copyright (C) 2019-2024 by Chinedu Ekuma                            #
#                                                                       #
#   This program is free software; you can redistribute it and/or       #
#   modify it under the terms of the GNU General Public License         #
#   as published by the Free Software Foundation, version 3 of          #
#   the License.                                                        #
#                                                                       #
#   This program is distributed in the hope that it will be useful,     #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of      #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       #
#   GNU General Public License for more details.                        #
#                                                                       #
#   E-mail: cekuma1@gmail.com                                           #
#                                                                       #
#########################################################################
reset
set terminal pngcairo truecolor enhanced font 'Arial, 22' size 2250,750
set output "ray_surface_enh.png"

set palette defined (-1 "blue", 0 "white", 1 "red", 2 "black")
set pm3d depthorder explicit
set view equal xyz
set xyplane 0
unset border
unset xtics
unset ytics
unset ztics
set cblabel "log_{10}(A)"
#set cbtics 0.5
set view 60,120,1.2

set multiplot layout 1,3

set cbrange [-1.0:2.0]
set title "Slow Secondary"
splot "slow_secondary.dat" u 13:14:15:(log10($17)) w pm3d notitle, "slow_secondary.dat" u (-$13):(-$14):(-$15):(log10($17)) w pm3d notitle

set title "Fast Secondary"
splot "fast_secondary.dat" u 13:14:15:(log10($17)) w pm3d notitle, "fast_secondary.dat" u (-$13):(-$14):(-$15):(log10($17)) w pm3d notitle

set cbrange [-0.6:1.2]
set title "Primary"
splot "primary.dat" u 13:14:15:(log10($17)) w pm3d notitle, "primary.dat" u (-$13):(-$14):(-$15):(log10($17)) w pm3d notitle

unset multiplot
unset output

reset
"""


ray_surface = """
#########################################################################
#                                                                       #
#   Elastool -- Elastic toolkit for zero and finite-temperature         #
#   elastic constants and mechanical properties calculations            #
#                                                                       #
#   Modified from the Christofel software                               #
#   Copyright (C) 2019-2024 by Chinedu Ekuma                            #
#                                                                       #
#   This program is free software; you can redistribute it and/or       #
#   modify it under the terms of the GNU General Public License         #
#   as published by the Free Software Foundation, version 3 of          #
#   the License.                                                        #
#                                                                       #
#   This program is distributed in the hope that it will be useful,     #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of      #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       #
#   GNU General Public License for more details.                        #
#                                                                       #
#   E-mail: cekuma1@gmail.com                                           #
#                                                                       #
#########################################################################
reset
set terminal pngcairo truecolor enhanced font 'Arial, 22' size 2250,750
set output "ray_surface.png"

set palette defined (-1 "blue", 0 "white", 1 "red")
set pm3d depthorder explicit
set view equal xyz
set xyplane 0
unset border
unset xtics
unset ytics
unset ztics
set cblabel "v_g (km/s)"
#set cbtics 0.5
set view 60,120,1.2

set multiplot layout 1,3

set title "Slow Secondary"
splot "slow_secondary.dat" u 13:14:15:11 w pm3d notitle, "slow_secondary.dat" u (-$13):(-$14):(-$15):11 w pm3d notitle

set title "Fast Secondary"
splot "fast_secondary.dat" u 13:14:15:11 w pm3d notitle, "fast_secondary.dat" u (-$13):(-$14):(-$15):11 w pm3d notitle

set title "Primary"
splot "primary.dat" u 13:14:15:11 w pm3d notitle, "primary.dat" u (-$13):(-$14):(-$15):11 w pm3d notitle

unset multiplot
unset output

reset
"""

group_eqar = """
#########################################################################
#                                                                       #
#   Elastool -- Elastic toolkit for zero and finite-temperature         #
#   elastic constants and mechanical properties calculations            #
#                                                                       #
#   Modified from the Christofel software                               #
#   Copyright (C) 2019-2024 by Chinedu Ekuma                            #
#                                                                       #
#   This program is free software; you can redistribute it and/or       #
#   modify it under the terms of the GNU General Public License         #
#   as published by the Free Software Foundation, version 3 of          #
#   the License.                                                        #
#                                                                       #
#   This program is distributed in the hope that it will be useful,     #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of      #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       #
#   GNU General Public License for more details.                        #
#                                                                       #
#   E-mail: cekuma1@gmail.com                                           #
#                                                                       #
#########################################################################
reset
set terminal pngcairo truecolor enhanced font 'Arial, 22' size 2250,750
set output "group_velocity_eqar.png"

set palette defined (-1 "blue", 0 "white", 1 "red")
set view equal xyz
set xyplane 0
unset border
unset xtics
unset ytics
unset ztics
set cblabel "v_g (km/s)"
#set cbtics 0.5  # Adjust the interval as needed

set arrow from first 1,0 to first 1.2,0 back
set arrow from first 0,1 to first 0,1.25 back

set label 'X' at first 1.1,0.1 front 
set label 'Y' at first 0.1,1.15 front

set view 0,0,1.1

set multiplot layout 1,3

set title "Slow Secondary"
splot "slow_secondary.dat" u (sqrt(2)*sin(0.5*$1)*cos($2)):(sqrt(2)*sin(0.5*$1)*sin($2)):11:11 w pm3d notitle

set title "Fast Secondary"
splot "fast_secondary.dat" u (sqrt(2)*sin(0.5*$1)*cos($2)):(sqrt(2)*sin(0.5*$1)*sin($2)):11:11 w pm3d notitle

set title "Primary"
splot "primary.dat" u (sqrt(2)*sin(0.5*$1)*cos($2)):(sqrt(2)*sin(0.5*$1)*sin($2)):11:11 w pm3d notitle

unset multiplot
unset output

reset
"""

group_rel_eqar = """
#########################################################################
#                                                                       #
#   Elastool -- Elastic toolkit for zero and finite-temperature         #
#   elastic constants and mechanical properties calculations            #
#                                                                       #
#   Modified from the Christofel software                               #
#   Copyright (C) 2019-2024 by Chinedu Ekuma                            #
#                                                                       #
#   This program is free software; you can redistribute it and/or       #
#   modify it under the terms of the GNU General Public License         #
#   as published by the Free Software Foundation, version 3 of          #
#   the License.                                                        #
#                                                                       #
#   This program is distributed in the hope that it will be useful,     #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of      #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       #
#   GNU General Public License for more details.                        #
#                                                                       #
#   E-mail: cekuma1@gmail.com                                           #
#                                                                       #
#########################################################################
reset
set terminal pngcairo truecolor enhanced font 'Arial, 22' size 2250,750
set output "group_velocity_relative_eqar.png"

# Get min and max values from the data files
stats "slow_secondary.dat" u 12 nooutput
S_min = floor(STATS_min)
S_max = ceil(STATS_max)

stats "fast_secondary.dat" u 12 nooutput
if (floor(STATS_min) < S_min) S_min = floor(STATS_min)
if (ceil(STATS_max) > S_max) S_max = ceil(STATS_max)

stats "primary.dat" u 12 nooutput
P_min = floor(STATS_min)
P_max = ceil(STATS_max)

if (-S_min > S_max) S_max = -S_min
if (-P_min > P_max) P_max = -P_min

set palette defined (-1 "blue", 0 "white", 1 "red")
set pm3d depthorder explicit
set view equal xyz
set xyplane 0
unset border
unset xtics
unset ytics
unset ztics
set cblabel "v_g - v_{iso} (%)"

set arrow from first 1,0 to first 1.2,0 back
set arrow from first 0,1 to first 0,1.25 back

set label 'X' at first 1.1,0.1 front 
set label 'Y' at first 0.1,1.15 front

set view 0,0,1.1

set multiplot layout 1,3

set title "Slow Secondary"
set cbrange [-S_max:S_max]
splot "slow_secondary.dat" u (sqrt(2)*sin(0.5*$1)*cos($2)):(sqrt(2)*sin(0.5*$1)*sin($2)):12:12 w pm3d notitle

set title "Fast Secondary"
splot "fast_secondary.dat" u (sqrt(2)*sin(0.5*$1)*cos($2)):(sqrt(2)*sin(0.5*$1)*sin($2)):12:12 w pm3d notitle

set title "Primary"
set cbrange [-P_max:P_max]
splot "primary.dat" u (sqrt(2)*sin(0.5*$1)*cos($2)):(sqrt(2)*sin(0.5*$1)*sin($2)):12:12 w pm3d notitle

unset multiplot
unset output

reset
"""

group_rel_sphere = """
#########################################################################
#                                                                       #
#   Elastool -- Elastic toolkit for zero and finite-temperature         #
#   elastic constants and mechanical properties calculations            #
#                                                                       #
#   Modified from the Christofel software                               #
#   Copyright (C) 2019-2024 by Chinedu Ekuma                            #
#                                                                       #
#   This program is free software; you can redistribute it and/or       #
#   modify it under the terms of the GNU General Public License         #
#   as published by the Free Software Foundation, version 3 of          #
#   the License.                                                        #
#                                                                       #
#   This program is distributed in the hope that it will be useful,     #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of      #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       #
#   GNU General Public License for more details.                        #
#                                                                       #
#   E-mail: cekuma1@gmail.com                                           #
#                                                                       #
#########################################################################
reset
set terminal pngcairo truecolor enhanced font 'Arial, 22' size 2250,750
set output "group_velocity_relative_sphere.png"

# Get min and max values from the data files
stats "slow_secondary.dat" u 12 nooutput
S_min = floor(STATS_min)
S_max = ceil(STATS_max)

stats "fast_secondary.dat" u 12 nooutput
if (floor(STATS_min) < S_min) S_min = floor(STATS_min)
if (ceil(STATS_max) > S_max) S_max = ceil(STATS_max)

stats "primary.dat" u 12 nooutput
P_min = floor(STATS_min)
P_max = ceil(STATS_max)

if (-S_min > S_max) S_max = -S_min
if (-P_min > P_max) P_max = -P_min

set palette defined (-1 "blue", 0 "white", 1 "red")
set pm3d depthorder explicit
set view equal xyz
set xyplane 0
unset border
unset xtics
unset ytics
unset ztics
set cblabel "v_g - v_{iso} (%)"

set view 60,120,1.3
set arrow from first 0,0,1 to first 0,0,1.35 front
set arrow from first 0,1,0 to first 0,1.4,0 front
set arrow from first 1,0,0 to first 1.4,0,0 front

set label 'X' at first 1.7,0,0 front 
set label 'Y' at first 0.2,1.35,0 front
set label 'Z' at first 0,0.07,1.3 front

set multiplot layout 1,3

set title "Slow Secondary"
set cbrange [-S_max:S_max]
splot "slow_secondary.dat" u (sin($1)*cos($2)):(sin($1)*sin($2)):(cos($1)):12 w pm3d notitle, \
      "slow_secondary.dat" u (-sin($1)*cos($2)):(-sin($1)*sin($2)):(-cos($1)):12 w pm3d notitle

set title "Fast Secondary"
splot "fast_secondary.dat" u (sin($1)*cos($2)):(sin($1)*sin($2)):(cos($1)):12 w pm3d notitle, \
      "fast_secondary.dat" u (-sin($1)*cos($2)):(-sin($1)*sin($2)):(-cos($1)):12 w pm3d notitle

set title "Primary"
set cbrange [-P_max:P_max]
splot "primary.dat" u (sin($1)*cos($2)):(sin($1)*sin($2)):(cos($1)):12 w pm3d notitle, \
      "primary.dat" u (-sin($1)*cos($2)):(-sin($1)*sin($2)):(-cos($1)):12 w pm3d notitle

unset multiplot
unset output

reset
"""

group_sphere = """
#########################################################################
#                                                                       #
#   Elastool -- Elastic toolkit for zero and finite-temperature         #
#   elastic constants and mechanical properties calculations            #
#                                                                       #
#   Modified from the Christofel software                               #
#   Copyright (C) 2019-2024 by Chinedu Ekuma                            #
#                                                                       #
#   This program is free software; you can redistribute it and/or       #
#   modify it under the terms of the GNU General Public License         #
#   as published by the Free Software Foundation, version 3 of          #
#   the License.                                                        #
#                                                                       #
#   This program is distributed in the hope that it will be useful,     #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of      #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       #
#   GNU General Public License for more details.                        #
#                                                                       #
#   E-mail: cekuma1@gmail.com                                           #
#                                                                       #
#########################################################################
reset
set terminal pngcairo truecolor enhanced font 'Arial, 22' size 2250,750
set output "group_velocity_sphere.png"

set palette defined (-1 "blue", 0 "white", 1 "red")
set pm3d depthorder explicit
set view equal xyz
set xyplane 0
unset border
unset xtics
unset ytics
unset ztics
set cblabel "v_g (km/s)"
#set cbtics 0.5  # Adjust the interval as needed

set view 60,120,1.3
set arrow from first 0,0,1 to first 0,0,1.35 front
set arrow from first 0,1,0 to first 0,1.4,0 front
set arrow from first 1,0,0 to first 1.4,0,0 front

set label 'X' at first 1.7,0,0 front 
set label 'Y' at first 0.2,1.35,0 front
set label 'Z' at first 0,0.07,1.3 front

set multiplot layout 1,3

set title "Slow Secondary"
splot "slow_secondary.dat" u (sin($1)*cos($2)):(sin($1)*sin($2)):(cos($1)):11 w pm3d notitle, \
      "slow_secondary.dat" u (-sin($1)*cos($2)):(-sin($1)*sin($2)):(-cos($1)):11 w pm3d notitle

set title "Fast Secondary"
splot "fast_secondary.dat" u (sin($1)*cos($2)):(sin($1)*sin($2)):(cos($1)):11 w pm3d notitle, \
      "fast_secondary.dat" u (-sin($1)*cos($2)):(-sin($1)*sin($2)):(-cos($1)):11 w pm3d notitle

set title "Primary"
splot "primary.dat" u (sin($1)*cos($2)):(sin($1)*sin($2)):(cos($1)):11 w pm3d notitle, \
      "primary.dat" u (-sin($1)*cos($2)):(-sin($1)*sin($2)):(-cos($1)):11 w pm3d notitle

unset multiplot
unset output

reset

"""


pol_sphere = """
#########################################################################
#                                                                       #
#   Elastool -- Elastic toolkit for zero and finite-temperature         #
#   elastic constants and mechanical properties calculations            #
#                                                                       #
#   Modified from the Christofel software                               #
#   Copyright (C) 2019-2024 by Chinedu Ekuma                            #
#                                                                       #
#   This program is free software; you can redistribute it and/or       #
#   modify it under the terms of the GNU General Public License         #
#   as published by the Free Software Foundation, version 3 of          #
#   the License.                                                        #
#                                                                       #
#   This program is distributed in the hope that it will be useful,     #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of      #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       #
#   GNU General Public License for more details.                        #
#                                                                       #
#   E-mail: cekuma1@gmail.com                                           #
#                                                                       #
#########################################################################
reset
set terminal pngcairo truecolor enhanced font 'Arial, 22' size 2250,750
set output "phase_polarization.png"

set palette defined (-1 "blue", 0 "white", 1 "red")
set pm3d depthorder explicit
set view equal xyz
set xyplane 0
set hidden3d front
unset border
unset xtics
unset ytics
unset ztics
set cblabel "v_p (km/s)"

# Set the number of ticks on the color bar
set cbtics 1.0  # Adjust the interval as neede

set view 60,120,1.3
set arrow from first 0,0,1 to first 0,0,1.35 front lw 2
set arrow from first 0,1,0 to first 0,1.4,0 front lw 2
set arrow from first 1,0,0 to first 1.4,0,0 front lw 2

set label 'X' at first 1.7,0,0 front 
set label 'Y' at first 0.2,1.35,0 front
set label 'Z' at first 0,0.07,1.3 front

set multiplot layout 1,3

# Determine the min and max values for the color bar range
min_max_range(file) = system(sprintf("awk '{if(NR==1){min=$6; max=$6} else {if($6<min){min=$6} if($6>max){max=$6}}}; END{print min, max}' %s", file))
set cbrange [*:*]

set title "Slow Secondary"
splot "slow_secondary.dat" u (0.98*sin($1)*cos($2)):(0.98*sin($1)*sin($2)):(0.98*cos($1)):6:6 w p palette ps 2 pt 3 notitle, \
      "slow_secondary.dat" u (-0.98*sin($1)*cos($2)):(-0.98*sin($1)*sin($2)):(-0.98*cos($1)):6:6 w p palette ps 2 pt 3 notitle

set title "Fast Secondary"
splot "fast_secondary.dat" u (0.98*sin($1)*cos($2)):(0.98*sin($1)*sin($2)):(0.98*cos($1)):6:6 w p palette ps 2 pt 3 notitle, \
      "fast_secondary.dat" u (-0.98*sin($1)*cos($2)):(-0.98*sin($1)*sin($2)):(-0.98*cos($1)):6:6 w p palette ps 2 pt 3 notitle

set title "Primary"
splot "primary.dat" u (0.98*sin($1)*cos($2)):(0.98*sin($1)*sin($2)):(0.98*cos($1)):6:6 w p palette ps 2 pt 3 notitle, \
      "primary.dat" u (-0.98*sin($1)*cos($2)):(-0.98*sin($1)*sin($2)):(-0.98*cos($1)):6:6 w p palette ps 2 pt 3 notitle

unset multiplot
unset output

reset

"""


pfangle_cube = """
#########################################################################
#                                                                       #
#   Elastool -- Elastic toolkit for zero and finite-temperature         #
#   elastic constants and mechanical properties calculations            #
#                                                                       #
#   Modified from the Christofel software                               #
#   Copyright (C) 2019-2024 by Chinedu Ekuma                            #
#                                                                       #
#   This program is free software; you can redistribute it and/or       #
#   modify it under the terms of the GNU General Public License         #
#   as published by the Free Software Foundation, version 3 of          #
#   the License.                                                        #
#                                                                       #
#   This program is distributed in the hope that it will be useful,     #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of      #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       #
#   GNU General Public License for more details.                        #
#                                                                       #
#   E-mail: cekuma1@gmail.com                                           #
#                                                                       #
#########################################################################
reset
set terminal pngcairo truecolor enhanced font 'Arial, 22' size 2250,750
set output "powerflow_cube.png"

set palette defined (0 "brown", 1 "green")
set pm3d depthorder explicit
set view equal xyz
set xyplane 0
unset border
unset xtics
unset ytics
unset ztics

set cblabel "PF angle (Deg)"

set view 60,120,0.95
set arrow from first -1,-1,+1 to first -1,-1,1.3 front
set arrow from first -1,+1,-1 to first -1,1.4,-1 front
set arrow from first +1,-1,-1 to first 1.5,-1,-1 front

set label 'X' at first 1.45,-0.8,-1 front 
set label 'Y' at first -0.65,1.3,-1 front
set label 'Z' at first -1.1,-0.9,1.2 front

set arrow from +1,+1,+1 to -1,+1,+1 nohead front
set arrow from +1,+1,+1 to +1,-1,+1 nohead front
set arrow from +1,+1,+1 to +1,+1,-1 nohead front

set arrow from -1,+1,+1 to -1,+1,-1 nohead front
set arrow from -1,+1,+1 to -1,-1,+1 nohead front

set arrow from +1,-1,+1 to -1,-1,+1 nohead front
set arrow from +1,-1,+1 to +1,-1,-1 nohead front

set arrow from +1,+1,-1 to -1,+1,-1 nohead front
set arrow from +1,+1,-1 to +1,-1,-1 nohead front

set multiplot layout 1,3

set title "Slow Secondary"
splot "slow_secondary.dat" u 3:4:5:16 w pm3d notitle, \
      "slow_secondary.dat" u (-$3):(-$4):(-$5):16 w pm3d notitle

set title "Fast Secondary"
splot "fast_secondary.dat" u 3:4:5:16 w pm3d notitle, \
      "fast_secondary.dat" u (-$3):(-$4):(-$5):16 w pm3d notitle

set title "Primary"
splot "primary.dat" u 3:4:5:16 w pm3d notitle, \
      "primary.dat" u (-$3):(-$4):(-$5):16 w pm3d notitle

unset multiplot
unset output

reset


"""

pfangle_eqar = """
#########################################################################
#                                                                       #
#   Elastool -- Elastic toolkit for zero and finite-temperature         #
#   elastic constants and mechanical properties calculations            #
#                                                                       #
#   Modified from the Christofel software                               #
#   Copyright (C) 2019-2024 by Chinedu Ekuma                            #
#                                                                       #
#   This program is free software; you can redistribute it and/or       #
#   modify it under the terms of the GNU General Public License         #
#   as published by the Free Software Foundation, version 3 of          #
#   the License.                                                        #
#                                                                       #
#   This program is distributed in the hope that it will be useful,     #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of      #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       #
#   GNU General Public License for more details.                        #
#                                                                       #
#   E-mail: cekuma1@gmail.com                                           #
#                                                                       #
#########################################################################
reset
set terminal pngcairo truecolor enhanced font 'Arial, 22' size 2250,750
set output "powerflow_eqar.png"

set palette defined (0 "gold", 1 "green")
set view equal xyz
set xyplane 0
unset border
unset xtics
unset ytics
unset ztics
set cblabel "PF angle (Deg)"

set arrow from first 1,0 to first 1.2,0 back
set arrow from first 0,1 to first 0,1.25 back

set label 'X' at first 1.1,0.1 front 
set label 'Y' at first 0.1,1.15 front

set view 0,0,1.1

set multiplot layout 1,3

set title "Slow Secondary"
splot "slow_secondary.dat" u (sqrt(2)*sin(0.5*$1)*cos($2)):(sqrt(2)*sin(0.5*$1)*sin($2)):16:16 w pm3d notitle

set title "Fast Secondary"
splot "fast_secondary.dat" u (sqrt(2)*sin(0.5*$1)*cos($2)):(sqrt(2)*sin(0.5*$1)*sin($2)):16:16 w pm3d notitle

set title "Primary"
splot "primary.dat" u (sqrt(2)*sin(0.5*$1)*cos($2)):(sqrt(2)*sin(0.5*$1)*sin($2)):16:16 w pm3d notitle

unset multiplot
unset output

reset
"""

pfangle_sphere = """
#########################################################################
#                                                                       #
#   Elastool -- Elastic toolkit for zero and finite-temperature         #
#   elastic constants and mechanical properties calculations            #
#                                                                       #
#   Modified from the Christofel software                               #
#   Copyright (C) 2019-2024 by Chinedu Ekuma                            #
#                                                                       #
#   This program is free software; you can redistribute it and/or       #
#   modify it under the terms of the GNU General Public License         #
#   as published by the Free Software Foundation, version 3 of          #
#   the License.                                                        #
#                                                                       #
#   This program is distributed in the hope that it will be useful,     #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of      #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       #
#   GNU General Public License for more details.                        #
#                                                                       #
#   E-mail: cekuma1@gmail.com                                           #
#                                                                       #
#########################################################################
reset
set terminal pngcairo truecolor enhanced font 'Arial, 22' size 2250,750
set output "powerflow_sphere.png"

set palette defined (0 "black", 1 "green")
set pm3d depthorder explicit
set view equal xyz
set xyplane 0
unset border
unset xtics
unset ytics
unset ztics
set cblabel "PF angle (Deg)"

set view 60,120,1.3
set arrow from first 0,0,1 to first 0,0,1.35 front
set arrow from first 0,1,0 to first 0,1.4,0 front
set arrow from first 1,0,0 to first 1.4,0,0 front

set label 'X' at first 1.7,0,0 front 
set label 'Y' at first 0.2,1.35,0 front
set label 'Z' at first 0,0.07,1.3 front

set multiplot layout 1,3

set title "Slow Secondary"
splot "slow_secondary.dat" u (sin($1)*cos($2)):(sin($1)*sin($2)):(cos($1)):16 w pm3d notitle, "slow_secondary.dat" u (-sin($1)*cos($2)):(-sin($1)*sin($2)):(-cos($1)):16 w pm3d notitle

set title "Fast Secondary"
splot "fast_secondary.dat" u (sin($1)*cos($2)):(sin($1)*sin($2)):(cos($1)):16 w pm3d notitle, "fast_secondary.dat" u (-sin($1)*cos($2)):(-sin($1)*sin($2)):(-cos($1)):16 w pm3d notitle

set title "Primary"
splot "primary.dat" u (sin($1)*cos($2)):(sin($1)*sin($2)):(cos($1)):16 w pm3d notitle, "primary.dat" u (-sin($1)*cos($2)):(-sin($1)*sin($2)):(-cos($1)):16 w pm3d notitle

unset multiplot
unset output

reset
"""

pfangle_stereo = """
#########################################################################
#                                                                       #
#   Elastool -- Elastic toolkit for zero and finite-temperature         #
#   elastic constants and mechanical properties calculations            #
#                                                                       #
#   Modified from the Christofel software                               #
#   Copyright (C) 2019-2024 by Chinedu Ekuma                            #
#                                                                       #
#   This program is free software; you can redistribute it and/or       #
#   modify it under the terms of the GNU General Public License         #
#   as published by the Free Software Foundation, version 3 of          #
#   the License.                                                        #
#                                                                       #
#   This program is distributed in the hope that it will be useful,     #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of      #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       #
#   GNU General Public License for more details.                        #
#                                                                       #
#   E-mail: cekuma1@gmail.com                                           #
#                                                                       #
#########################################################################
reset
set terminal pngcairo truecolor enhanced font 'Arial, 22' size 2250,750
set output "powerflow_stereo.png"

set palette defined (0 "brown", 1 "green")

set view equal xyz
set xyplane 0
unset border
unset xtics
unset ytics
unset ztics
set cblabel "PF angle (Deg)"

set arrow from first 1,0 to first 1.2,0 back
set arrow from first 0,1 to first 0,1.25 back

set label 'X' at first 1.1,0.1 front 
set label 'Y' at first 0.1,1.15 front

set view 0,0,1.1

set multiplot layout 1,3
set title "Slow Secondary"
splot "slow_secondary.dat" u (tan(0.5*$1)*cos($2)):(tan(0.5*$1)*sin($2)):16:16 w pm3d notitle;
set title "Fast Secondary"
splot "fast_secondary.dat" u (tan(0.5*$1)*cos($2)):(tan(0.5*$1)*sin($2)):16:16 w pm3d notitle;
set title "Primary"
splot "primary.dat" u (tan(0.5*$1)*cos($2)):(tan(0.5*$1)*sin($2)):16:16 w pm3d notitle;
unset multiplot
unset output

reset
"""


readmetxt = """
#########################################################################
#                                                                       #
#   Elastool -- Elastic toolkit for zero and finite-temperature         #
#   elastic constants and mechanical properties calculations            #
#                                                                       #
#   Copyright (C) 2019-2024 by Chinedu Ekuma                            #
#                                                                       #
#   This program is free software; you can redistribute it and/or       #
#   modify it under the terms of the GNU General Public License         #
#   as published by the Free Software Foundation, version 3 of          #
#   the License.                                                        #
#                                                                       #
#   This program is distributed in the hope that it will be useful,     #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of      #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       #
#   GNU General Public License for more details.                        #
#                                                                       #
#   E-mail: cekuma1@gmail.com                                           #
#                                                                       #
#########################################################################


This readme is provided to describe some of the key plots generated by ElasTool toolkit


For the Christoffel plots, you can modify the range of plots in the gnuplot and rerun the elastool code with the postprocessing option, i.e., run_mode=3

The solution of the Christoffel equations generates the following output files:

- `slow_secondary.dat`, `fast_secondary.dat`, `primary.dat`
These files provide acoustic data for various directions defined in the code. They offer insights into the slowest, intermediate, and fastest acoustic modes. 
Each file comprises 17 columns with the following data:
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


The following are the naming conventions for the generated Christoffel png plots:

- `phase` ----- Phase velocity
- `group` ----- Group velocity (absolute value only)
- `ray` ------- Ray surface (directional group velocity)
- `pfangle` --- Power flow angle
- `enh` ------- Enhancement factor - the acoustic energy carried away from a point source in different directions
- `pol` ------- Phase velocity polarization
- `relative` -- Relative to isotropic sound velocities (for phase and group)
- `sphere` ---- Projection onto the unit sphere
- `cube` ------ Projection onto the unit cube
- `eqar` ------ Equal area plane projection - Lambert azimuthal equal-area projection, which is Useful for comparing sound velocity behavior in different directions. Does not preserve angle
- `stereo` ---- Stereographic plane projection, preserving local shape. Useful for studying local shapes or lines of interest
- `radius` ---- For phase velocity, showing data on a sphere with radius scaled by absolute velocity

Other plots generated by the ElasTool toolkit are the followings:

- poisson_ration_contour_directional
- strain_energy_density_
- strain_energy_density_3Dprojection_ (Only for 3D)
- linear_compressibility
- EVGK_polar_
- EVGK_theta_
- EVGK_theta_2Dplotly_ (Same as EVGK_theta_ but outputs only for 2D materials)
- groupvelocity_heatmap_
- EVK_heatmap_
- poisson_ratio_contour_directional_
- EV_ploar_directional_ (Only for 2D materials)
- groupvelocity_contour_ (Only for 2D materials
- poisson_ratio_3Dprojection_ (Only for 3D materials)
- shearmodulus_3Dprojection_ (Only for 3D materials)
- youngmodulus_3Dprojection_ (Only for 3D materials)
"""


phasedir = """
#########################################################################
#                                                                       #
#   Elastool -- Elastic toolkit for zero and finite-temperature         #
#   elastic constants and mechanical properties calculations            #
#                                                                       #
#   Copyright (C) 2019-2024 by Chinedu Ekuma                            #
#                                                                       #
#   This program is free software; you can redistribute it and/or       #
#   modify it under the terms of the GNU General Public License         #
#   as published by the Free Software Foundation, version 3 of          #
#   the License.                                                        #
#                                                                       #
#   This program is distributed in the hope that it will be useful,     #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of      #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       #
#   GNU General Public License for more details.                        #
#                                                                       #
#   E-mail: cekuma1@gmail.com                                           #
#                                                                       #
#########################################################################
#Modify and copy this file into phase_directions.txt and rerun Elastool with run_mode = 3 or using the postprocessing "elastool -pp" 
1 0 0
0 1 0
0 0 1
-1 -1 0
"""



def write_gnuplot_scripts(script_dir):
    # Function to write a script to a file
    def write_script(filename, content):
        with open(os.path.join(script_dir, filename), 'w') as file:
            file.write(content)

#    if os.path.exists(script_dir):
#        shutil.rmtree(script_dir)
#    os.makedirs(script_dir)

    if not os.path.exists(script_dir):
        os.makedirs(script_dir)
        

    # Dictionary of script names and their corresponding content
    scripts = {
        "phase_cube.gnu": phase_cube_script,
        "phase_eqar.gnu": phase_eqar_script,
        "phase_radius.gnu": phase_radius,
        "phase_rel_cube.gnu": phase_rel_cube,
        "phase_rel_eqar.gnu": phase_rel_eqar,
        "phase_rel_sphere.gnu": phase_rel_sphere,
        "phase_rel_stereo.gnu": phase_rel_stereo,
        "phase_sphere.gnu": phase_sphere,
        "phase_stereo.gnu": phase_stereo,
        "enh_cube.gnu": enh_cube,
        "enh_eqar.gnu": enh_eqar,
        "enh_sphere.gnu": enh_sphere,
        "enh_stereo.gnu": enh_stereo,
        "slowness.gnu": slowness,
        "slowness_groupdir.gnu": slowness_groupdir,
        "ray_equar.gnu": ray_equar,
        "ray_stereo.gnu": ray_stereo,
        "ray_surface_enh.gnu": ray_surface_enh,
        "ray_surface.gnu": ray_surface,
        "group_eqar.gnu": group_eqar,
        "group_rel_eqar.gnu": group_rel_eqar,
        "group_rel_sphere.gnu": group_rel_sphere,
        "group_sphere.gnu": group_sphere,
        "pol_sphere.gnu": pol_sphere,
        "pfangle_cube.gnu": pfangle_cube,
        "pfangle_eqar.gnu": pfangle_eqar,
        "pfangle_sphere.gnu": pfangle_sphere,
        "pfangle_stereo.gnu": pfangle_stereo,
        "README_plots.txt":readmetxt,
        "phase_directions.dir":phasedir
    }

    # Write the scripts
    for filename, content in scripts.items():
        write_script(filename, content)

    print("All plots described in README_plots.txt inside property_plots\nScripts for Christoffel plots generated and in '{}' folder".format(script_dir))
    print("**********************************************************************************")


# Usage example
#write_gnuplot_scripts("plot_christoffel")





# Function to write a script to a file
#def write_script(filename, content):
#    with open(os.path.join(script_dir, filename), 'w') as file:
#        file.write(content)

## Write the scripts
#write_script("phase_cube.gnu", phase_cube_script)
#write_script("phase_eqar.gnu", phase_eqar_script)
#write_script("phase_radius.gnu", phase_radius)
#write_script("phase_rel_cube.gnu", phase_rel_cube)
#write_script("phase_rel_eqar.gnu", phase_rel_eqar)
#write_script("phase_rel_sphere.gnu", phase_rel_sphere)
#write_script("phase_rel_stereo.gnu", phase_rel_stereo)
#write_script("phase_sphere.gnu", phase_sphere)
#write_script("phase_stereo.gnu", phase_stereo)
#write_script("enh_cube.gnu", enh_cube)
#write_script("enh_eqar.gnu", enh_eqar)
#write_script("enh_sphere.gnu", enh_sphere)
#write_script("enh_stereo.gnu", enh_stereo)
#write_script("slowness.gnu", slowness)
#write_script("slowness_groupdir.gnu", slowness_groupdir)
#write_script("ray_equar.gnu", ray_equar)
#write_script("ray_stereo.gnu", ray_stereo)
#write_script("ray_surface_enh.gnu", ray_surface_enh)
#write_script("ray_surface.gnu", ray_surface)
#write_script("group_eqar.gnu", group_eqar)
#write_script("group_rel_eqar.gnu", group_rel_eqar)
#write_script("group_rel_sphere.gnu", group_rel_sphere)
#write_script("group_sphere.gnu", group_sphere)
#write_script("pol_sphere.gnu", pol_sphere)
#write_script("pfangle_cube.gnu", pfangle_cube)
#write_script("pfangle_eqar.gnu", pfangle_eqar)
#write_script("pfangle_sphere.gnu", pfangle_sphere)
#write_script("pfangle_stereo.gnu", pfangle_stereo)
##write_script("plot_sphere.gnu", plot_sphere)

#print("Gnuplot scripts for Christoffel plots generated in the 'plot_christoffel' directory.")



import os
import subprocess

def run_gnuplot_scripts(script_dir, result_dir, main_dir):
    """
    Run all Gnuplot scripts in the specified directory and move the generated PNG files to the result directory.

    :param script_dir: Directory containing the Gnuplot scripts.
    :param result_dir: Directory to save the results.
    :param main_dir: Directory where the data files are located.
    """
    # Create the result directory if it doesn't exist
    os.makedirs(result_dir, exist_ok=True)

    # Save the current working directory
    original_dir = os.getcwd()

    # Change to the main directory where data files are located
    os.chdir(main_dir)

    # Loop through all .gnu files in the script directory and run them with gnuplot
    for script in os.listdir(script_dir):
        if script.endswith(".gnu"):
            script_path = os.path.join(script_dir, script)
            print(f"Processing {script_path}...")
            subprocess.run(["gnuplot", script_path])

    # Move the generated PNG files to the result directory
    for file in os.listdir(main_dir):
        if file.endswith(".png"):
            os.rename(file, os.path.join(result_dir, file))

    # Change back to the original working directory
    os.chdir(original_dir)

    print(f"All scripts processed. Results are in {result_dir}.")

# Usage example
#script_dir = 'plot_christoffel'
#result_dir = 'christoffel_results'
#main_dir = '.'  # Assuming the current directory contains the data files
#run_gnuplot_scripts(script_dir, result_dir, main_dir)



