'''
readData_WATS.py
written by: Hannah Holt
last updated: 5/2/18


readData_WATS.py is a program that reads in ASCII data taken from the DE-2 WATS instrument to look at the ETA structure
 using neutral wind and temperature measurements. Each WATS .txt file contains data for one time stamp split between two
consecutive lines. No header lines.

MSIS temperature comparision can be added if wanted.


'''

from __future__ import print_function, division
import sys, os
import numpy as np
from pyglow.pyglow import Point             # for MSIS temperature
from datetime import datetime
import matplotlib.pyplot as plt

def print_header():
    print('\n**********************************************')
    print(os.path.basename(__file__))
    print('Written by: Hannah Holt')
    print('Last updated: 6/21/2018')
    print('**********************************************\n')

def plotData(x_values, y_values, filename="Default.txt"):

    print("Done reading files. Plotting...")

    name = filename.replace("orbit", "")
    name = name.replace(".txt", "")

    if (x_name == 'Latitude') or (x_name == 'Longitude'):
        x_label = 'DE-2 ' + x_name + ' [deg]'
    elif (x_name == 'LST') or (x_name == 'UT'):
        x_label = 'DE-2 ' + x_name + ' [hr]'
    elif x_name == 'Altitude':
        x_label = 'DE-2 ' + x_name + ' [km]'
    else:
        print('ERROR: Bad x axis label.')
        sys.exit()

    if y_name == 'Density':
        y_label = whatplot + ' ' + y_name + ' [1/cm^3]'
    elif y_name == 'Altitude':
        y_label = y_name + ' [km]'
    else:
        print('ERROR: Bad y axis label.')
        sys.exit()

    plt.plot(x_values, y_values, marker='.')


    plt.title("DE-2 Orbit " + name)
    #plt.xlim([-50, 50])
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.show()

    return 0

def make_date(orbitfile):
    '''
        This function makes a datetime object for a specific DE-2 orbit in order to run the MSIS model
        NOTE: NOT INCORPORATING UT YET

        :param orbitfile: the orbit filename
        :param UT: single UT value given in milliseconds

        :return: datetime object (year, month, day, hour, minute, second)
    '''

    if orbitfile == 'orbit1614.txt':
        year = 1981
        month = 11
        day = 20
    else:
        year = 1982
        if orbitfile != 'orbit3838.txt':
            month = 11
            if orbitfile != 'orbit7144.txt' and orbitfile != 'orbit7153.txt':
                day = 21
            elif orbitfile == 'orbit7144.txt':
                day = 19
            else:
                day = 20
        else:
            month = 4
            day = 17

    return datetime(year, month, day)

# ----------------------- VARIABLES TO CHANGE -------------------------------------------------------
# ----------------------------------------------------------------------------------------------------

filenames = ['1982_323_orbit7144.txt']       # 'orbit7153.txt'
whatplot = 'Temperature'            # can be Temperature, Vwind, or Hwind
x_name = "Latitude"                      # can be: Latitude, Longitude, LST, UT, or Altitude
y_name = "Temperature"              # can be: Density or Altitude

add_MSIS = False                   # add MSIS temperature to plot or not
# ----------------------------------------------------------------------------------------------------

print_header()

delimiter = '\t'  # file delimiter
hour = 3600000  # ms in an hr
numFiles = len(filenames)  # number of files to read

# check for good path
for i in range(0, numFiles):
    if not os.path.isfile(filenames[i]):
        print("ERROR. Read file", filenames[i], "does not exist.")
        sys.exit()

# read in file
for j in range(0, numFiles):
    fread = open(filenames[j], 'r')
    print('Opening file', filenames[j])

    # get number of lines and initialize data. Time is in UT [ms], densities are 1/cm^3, errors are %
    line_num = 0
    for line in fread:
        line_num += 1
    fread.close()

    data_lines = int(line_num/2)                                # each data time stamp is split between two liness)

    time = np.zeros(data_lines)
    dens = np.zeros(data_lines)                          # total neutral density [cm^-3]
    temp = np.zeros(data_lines)                             # corrected neutral temperature [K]
    temp_MSIS = np.zeros(data_lines)                        # to compare with MSIS temperature
    alt, lat = np.zeros(data_lines), np.zeros(data_lines)   # altitude and geo lat [km] and [deg]
    long, LST = np.zeros(data_lines), np.zeros(data_lines)  # geo longitude and local solar time [deg] and [hr]

    # now read in the data
    fread = open(filenames[j], 'r')

    counter = 0
    k = 0
    for line in fread:
        line = line.lstrip(' ')         # remove leading whitespace if we have it
        print('k =', k)

        if not (counter % 2):           # we are on first line of data split k % 2 = 0
            for i in range(0, 10):
                if i == 1:
                    time[k] = float((line.split(delimiter)[i])) / hour        # UT time in hours
                elif i == 4 and float(line.split(delimiter)[i]) > 0:
                    dens[k] = float(line.split(delimiter)[i])
                elif i == 6 and float(line.split(delimiter)[i]) > 0:
                    temp[k] = float(line.split(delimiter)[i])
                else: continue

        else:                   # second line of split
            for i in range(0, 10):
                if i == 3:
                    alt[k] = float(line.split(delimiter)[i])
                elif i == 4:
                    lat[k] = float(line.split(delimiter)[i])
                elif i == 5:
                    long[k] = float(line.split(delimiter)[i])
                elif i == 6:
                    LST[k] = float(line.split(delimiter)[i])
                else: continue
            k += 1          # only increment k once per two lines
        counter += 1
    fread.close()



    # Get MSIS temperature data
    if add_MSIS:
        date = make_date(filenames[j])
        for i in range(0, len(temp)):
            pt = Point(date, lat[i], long[i], alt[i])
            result = pt.run_msis()  # call MSIS for these Point objects
            temp_MSIS[i] = result.Tn_msis

    # plt the data....

    print('Done reading files. Plotting...')

    plt.title('DE2 Orbit 1614 - LT 19.0')
    plt.scatter(lat, temp, marker='.', label='DE-2')
    if add_MSIS:
        plt.scatter(lat, temp_MSIS, marker='.', label='MSIS')
    plt.ylabel('Neutral Temperature [K]')
    plt.xlabel('Latitude [deg]')
    plt.legend()
    #plt.ylim([850, 1425])
    #plt.xlim([40, -40])
    plt.show()





