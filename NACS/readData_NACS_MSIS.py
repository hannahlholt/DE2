'''
NEEDS PYTHON 2.7 TO COMPILE AND RUN WITH MSIS

readData.py reads in a specified text file containing UT, gas number densities, and density errors for a certain orbit of the DE-2 satellite.
Uses MSIS to normalize the density values to a common altitude
File contains values for O, He, N2, N and Ar as well as other orbit info
'''

from __future__ import print_function, division
import sys, os
from pyglow.pyglow import Point             # for MSIS total mass density
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt

def print_header():
    print('\n**********************************************')
    print(os.path.basename(__file__))
    print('Written by: Hannah Holt')
    print('Last updated: 6/21/2018')
    print('**********************************************\n')

def plotData(x_values, y_values, filename="Default.txt", split = 0):
    '''
    This function plots the given x and y values for a given gas and orbit number file.

    :param x_values: array of values to plot on x axis. This could be lat, long, LST, etc.
    :param y_values: array of y values for the specific gas. This will usually be density
    :param combo: if this is true then the graph title will change. Defaults is no combo
    :param filename: name of input file for graph title
    :param split: if a value is given then we want to split up the graph between before and after going over the poles
    :return: plots the data and returns 0
    '''

    print("Done reading files. Plotting...")

    # GRAPH LABELING....
    name = filename.replace("orbit", "")
    name = name.replace(".txt", "")

    if normalize:
        name = name + ' NORMALIZED'

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

    if split:
        whichSplit = which_split(LST)

        if whichSplit == 1:
            # the FINAL values to plot!!
            x_plot = x_values[:split]  # takes first half of values up to the split (not including split)
            y_plot = y_values[:split]
            LST_plot = LST[split - 1]
            alt_plot = alt[:split]
            if add_MSIS:
                MSIS_plot = MSIS_rho[:split]

        elif whichSplit == 2:
            x_plot = x_values[split:]  # takes second half of values (includes split value)
            y_plot = y_values[split:]
            LST_plot = LST[split]
            alt_plot = alt[split:]
            if add_MSIS:
                MSIS_plot = MSIS_rho[split:]

    else:
        x_plot = x_values
        y_plot = y_values
        LST_plot = LST[0]
        alt_plot = alt
        if add_MSIS:
            MSIS_plot = MSIS_rho

    # Now plot.
    fig, ax1 = plt.subplots()
    if add_MSIS:
        ax1.plot(x_plot, y_plot, color='r', label='DE-2')
        ax1.plot(x_plot, MSIS_plot, color='g', label='MSIS')
    else:
        ax1.plot(x_plot, y_plot, color='r', label='DE-2')

    ax1.set_xlabel(x_label)
    ax1.set_ylabel(y_label)
    plt.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    plt.grid()

    if y_lower and y_upper:
        plt.ylim([y_lower, y_upper])

    plt.legend()

    ax2 = ax1.twinx()
    ax2.plot(x_plot, alt_plot, color='b', linestyle='--', label='DE-2 Altitude')
    ax2.set_ylabel('Altitude [km]')

    if x_lower and x_upper:
        plt.xlim([x_lower, x_upper])

    plt.title('DE2 Orbit ' + name + ' - LT ' + str(LST_plot))
    plt.show()


    return 0

def normalize_rho(orbitfile, lat, lon, alt, rho, base_alt):
    """
    Uses the MSIS empirical model to return the helium and total mass density in an array

    Arguments:
        orbitfile:  the filename of the DE-2 orbit
        lat     :  array of latitudes  [deg]
        long    :  array of longitudes [deg]
        alt     :  array of altitudes above body surface [km]
        rho     :  array of unnormalized number density values for helium (1/cm^3)
        base_alt:  the fixed altitude to normalize density value

    Returns:
        rho_matrix  : a 2Xn matrix of normalized number densities at the fixed "base altitude" and the expected MSIS value
                      densites are first row - i.e. rho_matrix[0:all]
                      MSIS expected values are second row - i.e. rho_matrix[1:all]
    """
    species = whatplot
    rho_norm = np.zeros(len(rho))
    rho_MSIS_fixed, rho_MSIS_vary = np.zeros(len(rho)), np.zeros(len(rho))
    ratio = np.zeros(len(rho))
    date = make_date(orbitfile)

    for i in range(0, len(rho)):
        # MSIS Point objects - pt1 = fixed alt, pt2 = variable alt
        pt1 = Point(date, lat[i], lon[i], base_alt)
        pt2 = Point(date, lat[i], lon[i], alt[i])

        result1, result2 = pt1.run_msis(), pt2.run_msis()                       # call MSIS for these Point objects
        rho_MSIS_fixed[i] = result1.nn[species]     # particles/cm^3
        rho_MSIS_vary[i] = result2.nn[species]
        ratio[i] = rho_MSIS_fixed[i]/rho_MSIS_vary[i]
        rho_norm[i] = rho[i] * ratio[i]



    plt.subplot(211)
    plt.title(filenames[j])
    plt.plot(alt, rho_MSIS_fixed, label='MSIS rho fixed')
    plt.plot(alt, rho_MSIS_vary, label='MSIS rho vary')
    plt.plot(alt, rho, label='rho')
    plt.plot(alt, rho_norm, label='rho norm')

    plt.legend(loc=1)
    plt.ylabel(whatplot + ' Density 1/cm^3')
    plt.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    plt.grid()

    plt.subplot(212)
    plt.plot(alt, ratio)
    plt.xlabel('Altitude')
    plt.ylabel(r'Ratio MSIS $\rho(z_0))/ \rho(z)$')
    plt.minorticks_on()
    plt.grid()
    plt.show()

    rho_matrix = np.array([rho_norm, rho_MSIS_fixed])
    return rho_matrix         # return both so you can plot

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

def avg(array, num_points):
    '''
    This function averages the given array over the specified number of points
    :param array: a numpy array of values
    :param num_points: an interger that tells how many points you would like to average over.
    :return: an array with len(array) - num_points + 1 values
    '''
    average = np.zeros(len(array) - num_points + 1)

    for i in range(0, len(average)):
        sum = 0.0
        for j in range(0, num_points):
            sum += array[i+j]
        average[i] = sum/num_points

    return average

def find_split(LST_arr):
    '''
    Since DE2 is a polar orbiter, it will drastically change local times (or longitude) when it passes over the poles. We want to find when this happens
    :param LST_arr: the array of LT values for a specific orbit. Find when LT_n+1 - LT_n > x.
    :return: an integer that marks the first location AFTER the split in the array from 0 to len(arr)
    '''

    dLT = 1                 # the LT difference cutoff. I arbitrarily choose this
    split = 0
    for i in range(1, len(LST_arr)):
            if (abs(LST_arr[i] - LST_arr[i-1]) > dLT):
                split = i
                break
    print("LT split = ", LST_arr[split])
    return split

def which_split(LST_arr):
    split = find_split(LST_arr)

    if LST_arr[split] > 10.0:           # we want late local times
        return 2
    else:
        return 1


# ----------------------- VARIABLES TO CHANGE -------------------------------------------------------
# ----------------------------------------------------------------------------------------------------

infile_directory = './NACSfiles/outputs/'
filenames = []
whatplot = "HE"                     # to specify what species to plot, can be HE, N2, N, O, or AR
x_name = "Latitude"                 # can be: Latitude, Longitude, LST, UT, or Altitude
y_name = "Density"                  # can be: Density or Altitude

normalize = True                    # normalize densities to common value (USES MSIS...)
average = False
base_alt = 325                      # altitude to normalize the densities too [km]
x_lower = 50                        # upper and lower x limits for graphing (0 for default) and same for y
x_upper = -50
y_lower = 0
y_upper = 0
add_MSIS = True                     # if you want to add MSIS to final plot to compare graphically
# ----------------------------------------------------------------------------------------------------

print_header()

# get all NACS files in folder
for file in os.listdir(infile_directory):
    filenames.append(file)

delimiter = '\t'  # file delimiter
hour = 3600000  # ms in an hr
numFiles = len(filenames)  # number of files to read

for i in range(0, numFiles):
    if not os.path.isfile(infile_directory + filenames[i]):
        print("ERROR. Read file", filenames[i], "does not exist.")
        sys.exit()

for j in range(0, numFiles):
    fread = open(infile_directory + filenames[j], 'r')
    print('\nOpening file', filenames[j])

    # get number of lines and initialize data. Time is in UT [ms], densities are 1/cm^3, errors are %
    data_lines = 0
    for line in fread:
        data_lines += 1
    fread.close()
    data_lines -= 2                  # don't care about first two lines (header)

    time = np.zeros(data_lines)
    O_dens, O_err = np.zeros(data_lines), np.zeros(data_lines)
    N2_dens, N2_err = np.zeros(data_lines), np.zeros(data_lines)
    He_dens, He_err = np.zeros(data_lines), np.zeros(data_lines)
    N_dens, N_err = np.zeros(data_lines), np.zeros(data_lines)
    Ar_dens, Ar_err = np.zeros(data_lines), np.zeros(data_lines)
    alt, lat, long, LST = np.zeros(data_lines), np.zeros(data_lines), np.zeros(data_lines), np.zeros(data_lines)

    # now read in the data
    fread = open(infile_directory + filenames[j], 'r')

    # skip header
    fread.next()
    fread.next()

    k = 0
    # read in all the values from the file
    for line in fread:
        for i in range(0, 16):
            if i == 0:
                time[k] = float((line.split(delimiter)[i])) / hour
            elif i == 1:
                O_dens[k] = float(line.split(delimiter)[i])
            elif i == 2:
                O_err[k] = float(line.split(delimiter)[i])
            elif i == 3:
                N2_dens[k] = float(line.split(delimiter)[i])
            elif i == 4:
                N2_err[k] = float(line.split(delimiter)[i])
            elif i == 5:
                He_dens[k] = float(line.split(delimiter)[i])
            elif i == 6:
                He_err[k] = float(line.split(delimiter)[i])
            elif i == 7:
                N_dens[k] = float(line.split(delimiter)[i])
            elif i == 8:
                N_err[k] = float(line.split(delimiter)[i])
            elif i == 9:
                Ar_dens[k] = float(line.split(delimiter)[i])
            elif i == 10:
                Ar_err[k] = float(line.split(delimiter)[i])
            elif i == 11:
                continue  # don't care about orbit number
            elif i == 12:
                alt[k] = float(line.split(delimiter)[i])
            elif i == 13:
                lat[k] = float(line.split(delimiter)[i])
            elif i == 14:
                long[k] = float(line.split(delimiter)[i])
            elif i == 15:
                LST[k] = float(line.split(delimiter)[i])
            else:
                continue
        k += 1
    fread.close()

    # find the LT split
    split = find_split(LST)
    print("Split occured at:", split)

    # change x values to what you want to plot
    if x_name == 'LST':
        x_values = LST
    elif x_name == 'Longitude':
        x_values = long
    elif x_name == 'Latitude':
        x_values = lat
    elif x_name == 'Altitude':
        x_values = alt
    elif x_name == 'UT':
        x_values = time
    else:
        print('ERROR: Invalid x axis.')
        sys.exit()

    # same for y...
    if y_name == 'Density':
        if whatplot == 'HE':
            y_values = He_dens
        elif whatplot == 'O':
            y_values = O_dens
        elif whatplot == 'N2':
            y_values = N2_dens
        elif whatplot == 'N':
            y_values = N_dens
        elif whatplot == 'AR':
            y_values = Ar_dens
        else:
            print('ERROR: Invalid gas constituent.')
            sys.exit()
    elif y_name == 'Altitude':
        y_values = alt
    else:
        print('ERROR: Invalid y axis.')
        sys.exit()

    if normalize:
        rho_values = normalize_rho(filenames[j], lat, long, alt, y_values, base_alt)
        MSIS_rho = rho_values[1,:]              # MSIS expected values
        y_values = rho_values[0,:]              # DE normalized values

    plotData(x_values, y_values, filenames[j], split)






