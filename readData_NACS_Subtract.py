'''
readData_V3.py
written by: Hannah Holt

NEEDS PYTHON 2.7 TO COMPILE AND RUN

readData.py reads in a specified text file containing UT, gas number densities, and density errors for a certain orbit of the DE-2 satellite.
File contains values for O, He, N2, N and Ar as well as other orbit info

VERSION 2.0 UPDATES
    - removed anything that had to do with combining plots. We don't need this tool
    - changed the normalization feature to use a pure helium scale height instead of MSIS.

VERSION 3.0 UPDATES
    - after talking to Alan Burns and Wenbin Wang, suggestions were made to instead of normalizing He to a specific altitude, compare 2 days before
    and two days after a major ETA event and subtract those from the main event to see
    - added in the ability to plot helium densities vs invariant latitude. This gets rid of magnetic equator errors based on
    the orbits being on different longitudes (even though they have same LT)
'''

from __future__ import print_function, division
import sys, os
from pyglow.pyglow import Point  # for MSIS total mass density
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
import csv


def print_header():
    print('\n**********************************************')
    print(os.path.basename(__file__))
    print('Written by: Hannah Holt')
    print('Last updated: 6/21/2018')
    print('**********************************************\n')

def plotData(x_values, y_values, filenames, LT_name):
    '''
    This function plots the given x and y values for a given gas and orbit number file.

    :param x_values:    array of values to plot on x axis. This could be lat, lon, LST, etc.
    :param y_values:    array of y values for the specific gas. This will usually be density
    :param filenames:    array of input files for graph title

    :return: plots the data and returns 0
    '''

    print("Done reading files. Plotting...")

    name = []

    # GRAPH LABELING....
    for i in range(0, len(filenames)):
        temp = filenames[i]
        temp = temp[4:]
        temp = temp.replace("_orbit", " (#")
        temp = temp.replace(".txt", ")")
        temp = temp.replace("_", "Day ")
        name.append(temp)

    y_label = whatplot + ' ' + y_name + ' Difference ' + ' [1/cm^3]'
    x_label = x_name + ' [deg]'
    graphTitle = 'Density Change of ' + name[1] + ' - ' + name[0] + '\nLT ~ ' + str(LT_name)

    if average:
        plt.plot(x_values, y_values, label='Raw DE2')
        plt.plot(x_values_avg, y_values_avg, label='Avg')
    else:
        plt.plot(x_values, y_values, label='Raw DE2')

    if x_lower != 0:
        plt.xlim(x_lower, x_upper)
    if y_lower != 0:
        plt.ylim(y_lower, y_upper)

    plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(graphTitle)
    plt.legend()
    plt.grid()
    plt.show()

    return 0

def dens_difference(dens_1, dens_2, lat_1, lat_2, lon_1, lon_2, LST_1, LST_2, split_1, split_2):
    """
    Takes in two density arrays and finds the realtive difference between the values for the same latitude

    Arguments:
        dens_1  :   density arr of first file [cm^-3]
        dens_2  :   density arr of second file [cm^-3]
        lat_1   :   latitude arr of first file [deg]
        lat_2   :   latitude arr of second file [deg]
        LST_1   :   local time arr of first file
        LST_2   :   local time arr of second file
        split_1 :   LT split integer for first file
        split_2 :   LT split integer for second file

    Returns:
        lat_diff    :   even though the latitudes will be roughly the same for consecutive files, this returns the average
        dens_diff   :   relative difference between the density values
        LT_name     :   the average LT value for the two orbits. Used later for plotting
    """

    # check to see if we want MLAT instead
    if x_name == "Magnetic Latitude":
        print('Changing to Geomagnetic Coordinates...')
        lat_1 = get_MLAT(lat_1, lon_1)
        lat_2 = get_MLAT(lat_2, lon_2)

    # only use the LT around 18
    if LST_1[split_1] > 12:
        length_1 = len(LST_1[split_1:])
        start1 = split_1
    else:
        length_1 = len(LST_1[:split_1])
        start1 = 0
        print("*** UPDATE: Changing 1st file split...")
        print("New LT1 = ", LST_1[start1])

    if LST_2[split_2] > 12:
        length_2 = len(LST_2[split_2:])
        start2 = split_2
    else:
        length_2 = len(LST_2[:split_2])
        start2 = 0
        print("*** UPDATE: Changing 2nd file split...")
        print("New LT2 = ", LST_2[start2])


    # now compare lengths - difference density array will be size of the smaller one
    if length_1 < length_2:
        size = length_1
    else:
        size = length_2

    # print('1at1[start1], lat2[start2] = ', lat_1[start1], lat_2[start2])
    # print('size = ', size)
    # print('start1 start2 =', start1, start2)
    # print('length to the end of 1 =', len(dens_1[start1:]))
    # print('length to the end of 2 =', len(dens_2[start2:]))

    # need to find latitude bin number where they start to overlap
    if lat_1[start1] > lat_2[start2]:
        for i in range(0, size):
            if lat_1[start1 + i] <= lat_2[start2]:
                start1 += i
                if len(dens_1[start1:start1+size]) < size:
                    size = len(dens_1[start1:start1+size])                                 # size of overlap shrinks
                break
    else:
        for i in range(0, size):
            if lat_2[start2 + i] <= lat_1[start1]:
                start2 += i
                if len(dens_2[start2:start2+size]) < size:
                    size = len(dens_2[start2:start2+size])                                 # size of overlap shrinks
                break

    # print("\nNEW VALUES.....")
    # print('1at1[start1], lat2[start2] = ', lat_1[start1], lat_2[start2])
    # print('size = ', size)
    # print('start1 start2 =', start1, start2)
    # print('length to the end of 1 =', len(dens_1[start1:]))
    # print('length to the end of 2 =', len(dens_2[start2:]))
    #
    # print('lat2[start stop] = ', lat_2[start2: start2+size])
    # print('lat1[start stop] = ', lat_1[start1: start1 + size])

    dens_diff = np.zeros(size)                      # difference between densities at same latitude
    lat_overlap = np.zeros(size)
    LT_value = (LST_1[start1] + LST_2[start2])/2


    # now go through and iterate to find the differences
    for i in range(0, size):
        dens_diff[i] = dens_2[start2 + i] - dens_1[start1 + i]

        # DE satellite will have basically the same latitude for the overlapping parts
        lat_overlap[i] = 0.5 * (lat_2[start2 + i] + lat_1[start1 + i])

    if compare_dens:
        plt.plot(lat_1[start1: start1+size], dens_1[start1:start1+size], label='Orbit 1608')
        plt.plot(lat_2[start2: start2+size], dens_2[start2:start2+size], label='Orbit 1614')
        plt.xlabel('Latitude [deg]')
        plt.ylabel('He Density [1/cm^3]')
        plt.grid()
        plt.legend()
        plt.show()

    return lat_overlap, dens_diff, LT_value

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
    print("LT split = ", LST_arr[split])
    return split

def avg(array, num_points):
    '''
    This function averages the given array over the specified number of points
    :param array: a numpy array of values
    :param num_points: an integer that tells how many points you would like to average over.
    :return: an array with len(array) - num_points + 1 values
    '''
    average = np.zeros(len(array) - num_points + 1)

    for i in range(0, len(average)):
        sum = 0.0
        for j in range(0, num_points):
            sum += array[i+j]
        average[i] = sum/num_points

    return average

def check_lon(lon_1, lat_1, lon_2, lat_2):
    '''

    :param lon_1: longitude values of first orbit file
    :param lat_1: geographical latitude values of first orbit file
    :param lon_2: "       "          second orbit file
    :param lat_2: "       "          second orbit file
    :return: plots the  latitude of the magnetic equator and two orbits (where they overlap) as a function of longitude
    '''

    # read in the magnetic coordinates
    mag_filename = 'mag_eq_lon_lat.txt'
    mag_lon = []
    mag_lat = []
    with open(mag_filename) as csvfile:
        readCSV = csv.reader(csvfile, delimiter=',')
        for line in readCSV:
            mag_lon.append(float(line[0]))
            mag_lat.append(float(line[1]))
    csvfile.close()

    # plot
    fig, ax = plt.subplots()
    label = []                              # legend labels
    for i in range(0, len(filenames)):
        temp = filenames[i]
        temp = temp[4:]
        temp = temp.replace("_orbit", " (#")
        temp = temp.replace(".txt", ")")
        temp = temp.replace("_", "Day ")
        label.append(temp)

    plt.plot(mag_lon, mag_lat, label='Magnetic Eq.')
    plt.plot(lon_1, lat_1, label=label[0])
    plt.plot(lon_2, lat_2, label=label[1])
    plt.xlabel("Geographic Longitude [deg]")
    plt.ylabel("Geographic Latitude [deg]")
    plt.title('DE 2 Orbit Comparison')
    plt.xlim(-180, 180)
    plt.ylim(-90, 90)
    ax.set_xticks(np.arange(-180, 181, step=20))
    plt.grid()
    plt.legend()
    plt.show()

    return 0

def get_MLAT(DE_lat, DE_lon):
    '''

    :param DE_lat: array of DE geographical latitudes for the orbit
    :param DE_lon: corresponding DE geographical longitudes for each of the latitudes
    :return: MLAT:  an array that gives the corresponding magnetic latitude for each of the longitudes
    '''

    #coordinates for the dip equator
    mag_lon_eq = []
    mag_lat_eq = []

    MLAT = np.zeros(len(DE_lat))        # magnetic latitude values to return
    # lat_values = np.zeros(len(DE_lat))

    mag_filename = 'mag_eq_lon_lat.txt'
    with open(mag_filename) as csvfile:
        readCSV = csv.reader(csvfile, delimiter=',')
        for line in readCSV:
            mag_lon_eq.append(float(line[0]))
            mag_lat_eq.append(float(line[1]))
    csvfile.close()

    for i in range(0, len(DE_lon)):
        for j in range(0, len(mag_lon_eq)):
            if int(DE_lon[i]) == mag_lon_eq[j]:
                MLAT[i] = DE_lat[i] - mag_lat_eq[j]
                # lat_values[i] = mag_lat_eq[j]

    # plt.plot(DE_lon, DE_lat, label='DE Lat')
    # plt.plot(DE_lon, lat_values, label='Mag. Eq.')
    # plt.plot(DE_lon, MLAT, label='DE MLat')
    # plt.xlabel("Longitude [deg]")
    # plt.legend()
    # plt.grid()
    # plt.show()

    return MLAT


# ----------------------- VARIABLES TO CHANGE -------------------------------------------------------
# ----------------------------------------------------------------------------------------------------

directory = './NACSfiles/outputs/'
filenames = ['1983_030_orbit8273.txt', '1983_029_orbit8257.txt']                # subtracts first file data from second
whatplot = "HE"                 # to specify what species to plot, can be: O, N2, HE, or AR
x_name = "Latitude"             # can be: 'Latitude' or 'Magnetic Latitude'
y_name = "Density"              # can be: Density
average = True                  # if you want to average the density values
npts = 75                       # number of points to average over
compare_lon = True              # if you want to plot the longitudinal comparison
compare_dens = True             # if you want to compare the individual densities
x_lower = 40                    # upper and lower x limits for graphing, zero = autoscale
x_upper = -40
y_upper = 0
y_lower = 0
# ----------------------------------------------------------------------------------------------------

print_header()

delimiter = '\t'                # file delimiter
hour = 3600000                  # ms in an hr
numFiles = len(filenames)       # number of files to read

if numFiles != 2:
    print('ERROR: Program needs two data files to compare.')
    sys.exit()

for i in range(0, numFiles):
    if not os.path.isfile(directory + filenames[i]):
        print("ERROR. Read file", filenames[i], "does not exist.")
        sys.exit()


for j in range(0, numFiles):
    fread = open(directory + filenames[j], 'r')
    print('.....')
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
    alt, lat, lon, LST = np.zeros(data_lines), np.zeros(data_lines), np.zeros(data_lines), np.zeros(data_lines)


    # ----------- READ IN DATA-----------------------------------------------------------------
    # -----------------------------------------------------------------------------------------

    fread = open(directory + filenames[j], 'r')

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
                lon[k] = float(line.split(delimiter)[i])
            elif i == 15:
                LST[k] = float(line.split(delimiter)[i])
            # elif i == 18:
            #     inv_lat[k] = float(line.split(delimiter)[i])
            else:
                continue
        k += 1

    fread.close()

    # --------------------------------------------------------------------------------

    # we are on the first file
    if j == 0:
        alt_1 = alt
        lat_1 = lat
        lon_1 = lon
        LST_1 = LST
        split_1 = find_split(LST)
        print("First file split occurred at:", split_1, '\n')

        if y_name == 'Density':
            if whatplot == 'HE':
                dens_1 = He_dens
            elif whatplot == 'O':
                dens_1 = O_dens
            elif whatplot == 'N2':
                dens_1 = N2_dens
            elif whatplot == 'N':
                dens_1 = N_dens
            elif whatplot == 'Ar':
                dens_1 = Ar_dens
            else:
                print('ERROR: Invalid gas constituent.')
                sys.exit()
        else:
            print('ERROR: Invalid y axis.')
            sys.exit()

    # we are on the second file
    else:
        alt_2 = alt
        lat_2 = lat
        lon_2 = lon
        LST_2 = LST
        split_2 = find_split(LST)
        print("Second file split occurred at:", split_2, '\n')

        if y_name == 'Density':
            if whatplot == 'HE':
                dens_2 = He_dens
            elif whatplot == 'O':
                dens_2 = O_dens
            elif whatplot == 'N2':
                dens_2 = N2_dens
            elif whatplot == 'N':
                dens_2 = N_dens
            elif whatplot == 'Ar':
                dens_2 = Ar_dens
            else:
                print('ERROR: Invalid gas constituent.')
                sys.exit()
        else:
            print('ERROR: Invalid y axis.')
            sys.exit()



# end j in range (0, numfiles)


# ----------- POST PROCESSING + PLOTTING ---------------------------------------------------
#  -----------------------------------------------------------------------------------------

# now find the difference between the two density values
x_values, y_values, LT_name = dens_difference(dens_1, dens_2, lat_1, lat_2, lon_1, lon_2, LST_1, LST_2, split_1, split_2)

if compare_lon:
    check_lon(lon_1, lat_1, lon_2, lat_2)

if average:
    x_values_avg = avg(x_values, npts)
    y_values_avg = avg(y_values, npts)

plotData(x_values, y_values, filenames, LT_name)
