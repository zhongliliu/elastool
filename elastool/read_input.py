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

from os import getcwd
import os
import numpy as np
import pkg_resources
from datetime import datetime
import sys
from write_default_input import write_default_input,write_default_elastool_in,print_default_input_message_0,print_default_input_message_1,display_help
from elastool_elate_browser import ElateAutomation
from material_analysis import MaterialAnalysis

def read_elastic_tensor(file_name):
    # Reading the elastic tensor from the file
    elastic_tensor = np.loadtxt(file_name)
    return elastic_tensor
    
def read_rho_dim(file_name):
    with open(file_name, 'r') as file:
        # Skip the header line
        next(file)
        # Read the data line
        data_line = next(file).strip()

    # Split the line into components and check if it contains exactly two values
    parts = data_line.split()
    if len(parts) != 2:
        raise ValueError(f"Expected 2 values in the line, but got {len(parts)} values.")

    rho_str, dimensional = parts

    # Convert the numeric data back to a float
    rho = float(rho_str)

    return rho, dimensional

        
def process_parameters(input_string):
    parameters = input_string.strip().split(',')
    parameters = [param.strip() for param in parameters]
    return parameters

def read_input():
    cwd = getcwd()
    main_infile = open('%s/elastool.in' % cwd, 'r')
    line = main_infile.readline()
    global indict
    indict = {}
    while line:
        line = main_infile.readline()
        llist = line.split('=')
        
        if llist != ['\n'] and llist != ['']:
            if llist[0][0] != '#':
                inputlist = [i.strip().split() for i in llist]
                
                # Handle elateparameters differently
                if inputlist[0][0] == 'elateparameters':
                    indict['elateparameters'] = process_parameters(llist[1])
                    #parameters = llist[1].strip().split(',')
                    #parameters = [param.strip() for param in parameters]
                    #indict['elateparameters'] = parameters
                elif inputlist[0][0] == 'plotparameters':
                    indict['plotparameters'] = process_parameters(llist[1])
                else:
                    if inputlist[1] == []:
                        with open('../log.elastool', 'a') as logfile:
                            print >>logfile, "Please give the value(s) for: %s" % inputlist[0][0]
                    else:
                        indict[inputlist[0][0]] = inputlist[1]
                        
    run_mode_flag = (len(sys.argv) > 1 and sys.argv[1] == "-0") or ('run_mode' in indict and int(indict['run_mode'][0]) == 0)

    if 'method_stress_statistics' in indict and run_mode_flag:
        write_default_input(indict['method_stress_statistics'][0], cwd)
        print_default_input_message_1()
        sys.exit(0)
    return indict
    
       
def print_pp_message(message):
    # Split the message into lines
    lines = message.split('\n')
    
    # Find the maximum width of the lines
    max_width = max(len(line) for line in lines)
    
    # Print the top of the box
    print('+' + '-' * (max_width + 2) + '+')
    
    # Print each line in the box
    for line in lines:
        print('| ' + line.ljust(max_width) + ' |')
    
    # Print the bottom of the box
    print('+' + '-' * (max_width + 2) + '+')

message = f"Post-processing simulation with \n ElasTool Version {pkg_resources.get_distribution('ElasTool').version} computational toolkit \n for visualizing the elastic parameters of your material. \n Calculation ended at {datetime.now().strftime('%H:%M:%S')} on {datetime.now().strftime('%Y-%m-%d')}"


pp = False
plotly_flag = False  # This will store the state of -T/-F

if len(sys.argv) > 1:
    first_arg = sys.argv[1].lower()
    pp = first_arg in ["-pp", "-postprocess"]

    if pp and len(sys.argv) > 2:
        second_arg = sys.argv[2].lower()
        if second_arg == "-plotly":
            plotly_flag = True
        elif second_arg == "-noplotly":
            plotly_flag = False

def post_process():
    rho, dim = read_rho_dim("massdensity_dim.dat")
    elastic_tensor = read_elastic_tensor("elastic_tensor.dat")
    print_pp_message(message)
    
    if dim == "3D":
        analysis_instance = MaterialAnalysis(elastic_tensor, rho, plot=True, plotly=False)
        analysis_instance.plot_linear_compressibility_3D()
        analysis_instance.plot_orientation_dependent_3D()
        analysis_instance.plot_contour_polar_3D()
        analysis_instance.heatmap_3D()
        analysis_instance.plot_moduli_heatmaps()

                
    elif dim == "2D":
        analysis_instance = MaterialAnalysis(elastic_tensor, rho, plot=True, plotly=False)
        analysis_instance.plot_orientation_dependent_2D()
        analysis_instance.plot_contour_2D()
        # analysis_instance.plot_contour_polar_2D()
        analysis_instance.heatmap_2D()
    elif dim == "1D":
        print("1D systems generally do not have robust spatial dependence.")
    sys.exit(0)
        



if pp:
    post_process()




activate_flag = False               

if len(sys.argv) > 1:
    first_arg = sys.argv[1].lower()
    activate_flag = first_arg in ["-elate"] 

    

if activate_flag:         
    browser = input("Choose a browser; when done press Ctrl+C: (chrome, firefox, edge, safari): ").lower()
    elate_instance = ElateAutomation(browser_name=browser)
    elate_instance.run()
    sys.exit(0)


# Check if any argument is a help command
help_commands = (len(sys.argv) > 1 and (sys.argv[1] == "-help" or sys.argv[1] == "--help" or sys.argv[1] == "--h" or sys.argv[1] == "-h"))
if help_commands:
    display_help()
    sys.exit(0)



cwd = os.getcwd()
elastool_in_exists = os.path.exists(os.path.join(cwd, "elastool.in"))
run_mode_flag_elastoolin = (len(sys.argv) > 1 and sys.argv[1] == "-0")
if run_mode_flag_elastoolin and not elastool_in_exists:
  write_default_elastool_in(cwd)
  print_default_input_message_0()
  sys.exit(0)


indict = read_input()

