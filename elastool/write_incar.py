"""
  Elastool -- Elastic toolkit for zero and finite-temperature elastic constants and mechanical properties calculations

  Copyright (C) 2019-2024 by Zhong-Li Liu and Chinedu Ekuma

  This program is free software; you can redistribute it and/or modify it under the
  terms of the GNU General Public License as published by the Free Software Foundation
  version 3 of the License.

  This program is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
  PARTICULAR PURPOSE.  See the GNU General Public License for more details.

  E-mail: zl.liu@163.com or cekuma1@gmail.com
"""


def write_incar(step, cwd):
    infile = cwd+'/INCARs'
    outfile = 'INCAR'
    tag = 'Step:'
    is_write = False
    is_tag_line_wroten = False

    for line in open(infile, 'r'):
        if tag in line and step in line:
            is_write = True
            is_tag_line_wroten = True
            with open(outfile, 'w') as ofile:
                ofile.write(line)

        if tag in line and step not in line:
            is_write = False
        elif is_write and not is_tag_line_wroten:
            with open(outfile, 'a') as ofile:
                ofile.write(line)

        is_tag_line_wroten = False


if __name__ == '__main__':
    import os
    # make_inputfile('sys.cell', 's.cell', '#castepinput', 4)
    cwd = os.getcwd()
    write_incar('opt', cwd)
