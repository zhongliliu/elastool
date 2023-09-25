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
from write_default_input import write_default_input,write_default_elastool_in,print_default_input_message_0,print_default_input_message_1,display_help
import sys


def read_input():
    cwd = getcwd()
    main_infile = open('%s/elastool.in'%cwd, 'r')
    line = main_infile.readline()
    global indict
    indict = {}
    while line:
        line = main_infile.readline()
        llist = line.split('=')
        
        if llist != ['\n'] and llist != ['']:
            if llist[0][0] != '#':
                inputlist = [i.strip().split() for i in llist]
                if inputlist[1] == []:
                    with open('../log.elastool', 'a') as logfile:
                        print >>logfile, "Please give the value(s) for: %s" % inputlist[0][0]
                        
                else:
                    indict[inputlist[0][0]] = inputlist[1]
    run_mode_flag = (len(sys.argv) > 1 and sys.argv[1] == "-0") or ('run_mode' in indict and int(indict['run_mode'][0]) == 0)
      

  
    if 'method_stress_statistics' in indict and run_mode_flag: #if 'method_stress_statistics' in indict and 'run_mode' in indict and int(indict['run_mode'][0]) == 0:
        write_default_input(indict['method_stress_statistics'][0], cwd)
        print_default_input_message_1()
        sys.exit(0)
    return indict

# Check for help flags before reading input
#help_flags = {'-help', '--help', '-h', '--h'}
#if set(sys.argv) & help_flags:
#    print(read_help('help.txt'))
#    sys.exit(0)


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

