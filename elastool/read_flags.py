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
from datetime import datetime

def get_gnuplot_scripts_path():
    # Get the directory where the elastool package is installed
    try:
        package_dir = os.path.join(os.path.dirname(__file__), 'gnuplot_scripts') #pkg_resources.resource_filename('elastool', '')
        #package_dir = pkg_resources.resource_listdir('elastool', 'gnuplot_scripts')
        #return package_dir
    except ModuleNotFoundError:
        print("elastool package not found.")
        return None

    #gnuplot_scripts_path = os.path.join(package_dir, 'gnuplot_scripts')
    #gnuplot_scripts_path = os.path.join(os.path.dirname(__file__), 'gnuplot_scripts')
    return package_dir

    
#---------------------------------------------------------
#Print out citation
def print_boxed_message(ec_file=None):
    header_footer = "+" + "-" * 78 + "+"
    spacer = "| " + " " * 76 + " |"

    # List of lines to be printed
    lines = [
        (" * CITATIONS *", True),
        ("If you have used Elastool in your research, PLEASE cite:", False),
        ("", False),  # Space after the above line
        ("ElasTool: An automated toolkit for elastic constants calculation, ", False),
        ("Z.-L. Liu, C.E. Ekuma, W.-Q. Li, J.-Q. Yang, and X.-J. Li, ", False),
        ("Computer Physics Communications 270, 108180, (2022)", False),
        ("", False),

        ("", False),  # Blank line for separation
        ("Efficient prediction of temperature-dependent elastic and", False),
        ("mechanical properties of 2D materials, S.M. Kastuar, C.E. Ekuma, Z-L. Liu,", False),
        ("Nature Scientific Report 12, 3776 (2022)", False),
        
        ("", False),  # Blank line for separation
        ("ElasTool v3.0: Efficient computational and visualization", False),
        ("toolkit for elastic and mechanical properties of materials,", False),
        ("C.E. Ekuma, Z-L. Liu,", False),
        ("Computer Physics Communications 300, 109161, (2024)", False)
    ]

    def output_line(line):
        if ec_file:
            ec_file.write(line + "\n")
        else:
            print(line)

    output_line(header_footer)
    
    for line, underline in lines:
        centered_line = line.center(76)
        output_line("| " + centered_line + " |")
        
        if underline:
            underline_str = "-" * len(centered_line)
            output_line("| " + underline_str.center(76) + " |")

    # Print footer of the box
    output_line(header_footer)

max_width = len("|WARNING: This is an empirical approx; validity needs to be checked !! |")

def write_line(ec_file, content, padding=1, border_char="|", filler_char=" "):
    content_width = int(max_width) - (2 * int(padding)) - 2  # Subtract 2 for the border characters
    content = content[:content_width]  # Ensure content doesn't exceed the width
    line = border_char + filler_char*padding + content.ljust(content_width) + filler_char*padding + border_char
    if ec_file:
        ec_file.write(line + "\n")
    else:
        print(line)



def print_banner(version, startend_status, ec_file=None):
    # Get current date and time
    current_time = datetime.now().strftime('%H:%M:%S')
    current_date = datetime.now().strftime('%Y-%m-%d')
    conclusion_msg = f"Calculations {startend_status} at {current_time} on {current_date}"

    # Concatenate the message with the version info
    message = f"ElastoMechanical Simulations\nusing\nElasTool Version: {version}\n{conclusion_msg}"

    write_line(ec_file, '❤' * (max_width - 2), padding=0, border_char='❤', filler_char='❤')  # This will print a line of hearts
    for line in message.split('\n'):
        centered_line = line.center(max_width - 4)  # Subtract 4 for the two border characters and spaces at each end
        write_line(ec_file, centered_line, padding=1, border_char='❤')
    write_line(ec_file, '❤' * (max_width - 2), padding=0, border_char='❤', filler_char='❤')  # This will print another line of hearts



# Function to write the provided script to a file named 'run_christoffel.py'
def write_run_christoffel_script():
    script_content = """
import os
import subprocess

# Directory containing the gnuplot scripts
script_dir = "gnuplot_scripts"

# Directory to save the results
result_dir = "christoffel_results"

# Create the result directory if it doesn't exist
os.makedirs(result_dir, exist_ok=True)

# Change to the script directory
os.chdir(script_dir)

# Loop through all .gnu files and run them with gnuplot
for script in os.listdir():
    if script.endswith(".gnu"):
        print(f"Processing {script}...")
        subprocess.run(["gnuplot", script])

        # Move the generated PNG files to the result directory
        for file in os.listdir():
            if file.endswith(".png"):
                os.rename(file, os.path.join("..", result_dir, file))

print(f"All scripts processed. Results are in {result_dir}.")
"""

    with open("run_christoffel.py", "w") as file:
        file.write(script_content)

    print("run_christoffel.py script has been written.")



