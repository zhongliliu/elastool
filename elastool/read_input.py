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
                    
    return indict


indict = read_input()
