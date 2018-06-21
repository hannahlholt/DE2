'''
replace_whitespace_NACS.py
written by: Hannah Holt

replace_whitespace_NACS.py is a program that takes raw NACS text files downloaded from the FTP site and removes the whitespaces
in the numerical data and replaces iit with tabs for reading in later. This program makes the NACS file compatible with the readData_NACS.py program

'''


from __future__ import print_function

import re
import sys, os
from datetime import datetime

def print_header():
    print('\n******************************************')
    print(os.path.basename(__file__))
    print('Written by: Hannah Holt')
    print('Last Updated: 6/21/2018')
    print('******************************************\n')

print_header()

infile_directory = './NACSfiles/raw_files/'         # directory with raw NACS files
outfile_directory = './NACSfiles/outputs/'
endfile = '.ASC'                                    # endfile tag
input_files = []
output_files = []
header = []
delim_in = ' +'                        # input file delimiter to replace
delim_out = '\t'                       # output file delimiter

if not os.path.exists(infile_directory):
    print("ERROR. Input file directory does not exist.")
    sys.exit()

# if the output file doesnt exist then make it
if not os.path.exists(outfile_directory):
    os.mkdir(outfile_directory)

# get all NACS files in folder
for file in os.listdir(infile_directory):
    if file.endswith(endfile):
        input_files.append(file)
        year = file[:4]
        day = file[4:7]
        output_files.append(year + '_' + day + '_orbit')
    else: continue

numFiles = len(input_files)

print('Reading ' + str(numFiles) + ' files from:\t' + infile_directory)
print('_____________________________________________')

for j in range(0, numFiles):
    fread = open(infile_directory + input_files[j], 'r')
    print('\nOpening file:', input_files[j])

    # get two header lines first
    header.append(fread.readline())
    header.append(fread.readline())

    # get the orbit number from the file for outputting
    line = fread.readline()
    line = line.lstrip(' ')
    line = re.sub(delim_in, delim_out, line)
    orbit = line.split(delim_out)[11]
    output_files[j] += (orbit + '.txt')

    fwrite = open(outfile_directory + output_files[j], 'w+')
    print(header[0], file=fwrite, end='')
    print(header[1], file=fwrite, end='')
    print(line, file=fwrite, end='')

    for line in fread:
        line = line.lstrip(' ')
        line = re.sub(delim_in, delim_out, line)
        print(line, file=fwrite, end='')

    print('Wrote:\t\t', output_files[j])

    fread.close()
    fwrite.close()

print('_____________________________________________')
print('DONE.')
print(str(len(output_files)) + ' files wrote to:\t' + outfile_directory)

