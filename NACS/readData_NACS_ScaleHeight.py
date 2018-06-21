'''
NEEDS PYTHON 2.7 TO COMPILE AND RUN WITH MSIS

readData.py reads in a specified text file containing UT, gas number densities, and density errors for a certain orbit of the DE-2 satellite.
File contains values for O, He, N2, N and Ar as well as other orbit info

VERSION 2.0 UPDATES
    - removed anything that had to do with combining plots. We don't need this tool
    - changed the normalization feature to use a pure helium scale height instead of MSIS.
'''



from __future__ import print_function, division
import sys, os
from pyglow.pyglow import Point  # for MSIS total mass density
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt

def print_header():
    print('\n**********************************************')
    print(os.path.basename(__file__))
    print('Written by: Hannah Holt')
    print('Last updated: 6/21/2018')
    print('**********************************************\n')

def plotData(x_values, y_values, filename="Default.txt", split=0):
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
        if split:
            print('Splitting Local Times...')

            if whichSplit == 1:
                # the FINAL values to plot!!
                x_plot = x_values[:split]  # takes first half of values up to the split (not including split)
                y_plot = y_values[:split]
                LST_plot = LST[split - 1]
                alt_plot = alt[:split]

            elif whichSplit == 2:
                x_plot = x_values[split:]  # takes second half of values (includes split value)
                y_plot = y_values[split:]
                LST_plot = LST[split]
                alt_plot = alt[split:]

        else:
            x_plot = x_values
            y_plot = y_values
            LST_plot = LST[0]
            alt_plot = alt

        # Now plot.
        fig, ax1 = plt.subplots()
        ax1.plot(x_plot, y_plot, color='r', label='DE-2')

        ax1.set_xlabel(x_label)
        ax1.set_ylabel(y_label)
        plt.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
        plt.grid()

        if y_lower and y_upper:
            plt.ylim([y_lower, y_upper])

        ax2 = ax1.twinx()
        ax2.plot(x_plot, alt_plot, color='b', linestyle='--', label='DE-2 Altitude')
        ax2.set_ylabel('Altitude [km]')

        plt.title('DE2 Orbit ' + name + ' - LT ' + str(LST_plot))
        if x_lower and x_upper:
            plt.xlim([x_lower, x_upper])

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
        return_arr  : a 2Xn array of normalized helium number densities at the fixed "base altitude" and the MSIS ratio values
                      densites are first row - i.e. return_arr[0:all]
                      ratios are second row - i.e. return_arr[1:all]
    """
    R = 8.314               # gas const [J/mol K]
    g = 9.81                # acceleration due to grav. [m/s^2]
    T = DE_temp             # Temperature [K]

    if whatplot == 'HE':
        mass = 0.004            # molecular mass of Helium [kg/mol]
    elif whatplot == 'N':
        mass = 0.014
    elif whatplot == 'N2':
        mass = 0.028
    elif whatplot == 'O':
        mass = 0.016
    elif whatplot == 'O2':
        mass = 0.032
    else:
        mass = 40.0

    H = R*T/(mass*g) / 1000        # scale height [km]

    print('scale height', H)

    rho_norm = np.zeros(len(rho))
    ratio = np.zeros(len(rho))

    for i in range(0, len(rho)):

        z_prime = alt[i] - base_alt
        ratio[i] = np.exp(z_prime/H)

        rho_norm[i] = rho[i] * ratio[i]

    plt.subplot(211)
    plt.title(filenames[j])
    plt.plot(alt, rho, label='rho')
    plt.plot(alt, rho_norm, label='rho norm')
    plt.legend(loc=1)
    plt.ylabel(whatplot + ' Density 1/cm^3')
    plt.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    plt.grid()

    plt.subplot(212)
    plt.plot(alt, ratio)
    plt.xlabel('Altitude')
    plt.ylabel(r'$e^{(z - z_0)/H}$')
    plt.minorticks_on()
    plt.grid()
    plt.show()

    return rho_norm

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
            sum += array[i + j]
        average[i] = sum / num_points

    return average

def find_split(LST_arr):
    '''
    Since DE2 is a polar orbiter, it will drastically change local times (or longitude) when it passes over the poles. We want to find when this happens
    :param LST_arr: the array of LT values for a specific orbit. Find when LT_n+1 - LT_n > x.
    :return: an integer that marks the first location AFTER the split in the array from 0 to len(arr)
    '''

    dLT = 1  # the LT difference cutoff. I arbitrarily choose this
    split = 0
    for i in range(1, len(LST_arr)):
        if (abs(LST_arr[i] - LST_arr[i - 1]) > dLT):
            split = i
            break
    return split


# ----------------------- VARIABLES TO CHANGE -------------------------------------------------------
# ----------------------------------------------------------------------------------------------------


filenames = ['./NACSfiles/outputs/1983_028_orbit8232.txt']  # , 'orbit7175.txt', 'orbit7164.txt']#, 'orbit7153.txt', 'orbit7161.txt', 'orbit7164.txt', 'orbit7175.txt'] #['orbit1614.txt', 'orbit3838.txt', ,aph
whatplot = "HE"                 # to specify what species to plot
x_name = "Latitude"             # can be: Latitude, Longitude, LST, UT, or Altitude
y_name = "Density"              # can be: Density or Altitude

normalize = True                 # normalize densities to common value
average = False
base_alt = 325                  # altitude to normalize the densities too [km]
x_lower = 40                   # upper and lower x limits for graphing
x_upper = -40
DE_temp = 1200                  # temperature taken from WATS instrument onboard DE-2
y_lower = 0                     # use zero for autoscale
y_upper = 0
whichSplit = 2
# ----------------------------------------------------------------------------------------------------

print_header()

delimiter = '\t'  # file delimiter
hour = 3600000  # ms in an hr
numFiles = len(filenames)  # number of files to read

for i in range(0, numFiles):
    if not os.path.isfile(filenames[i]):
        print("ERROR. Read file", filenames[i], "does not exist.")
        sys.exit()


for j in range(0, numFiles):
    fread = open(filenames[j], 'r')
    print('Opening file', filenames[j])

    # get number of lines and initialize data. Time is in UT [ms], densities are 1/cm^3, errors are %
    data_lines = 0
    for line in fread:
        data_lines += 1
    fread.close()
    data_lines -= 2  # don't care about first two lines (header)

    time = np.zeros(data_lines)
    O_dens, O_err = np.zeros(data_lines), np.zeros(data_lines)
    N2_dens, N2_err = np.zeros(data_lines), np.zeros(data_lines)
    He_dens, He_err = np.zeros(data_lines), np.zeros(data_lines)
    N_dens, N_err = np.zeros(data_lines), np.zeros(data_lines)
    Ar_dens, Ar_err = np.zeros(data_lines), np.zeros(data_lines)
    alt, lat, long, LST = np.zeros(data_lines), np.zeros(data_lines), np.zeros(data_lines), np.zeros(data_lines)

    # now read in the data
    fread = open(filenames[j], 'r')

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
        elif whatplot == 'Ar':
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
        y_values = normalize_rho(filenames[j], lat, long, alt, y_values, base_alt)

    plotData(x_values, y_values, filenames[j], split)


