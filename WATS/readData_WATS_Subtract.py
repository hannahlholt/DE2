'''
readData_WATS.py
written by: Hannah Holt
last updated: 5/2/18


readData_WATS.py is a program that reads in ASCII data taken from the DE-2 WATS instrument to look at the ETA structure
 using neutral wind and temperature measurements. Each WATS .txt file contains data for one time stamp split between two
consecutive lines. No header lines.

VERSION 2.0 UPDATES
    The main update is the program now reads in two WATS data files and looks at the temperature differences between them. T
    The first file is subtracted from the second file. MSIS features have been removed.

'''

from __future__ import print_function, division
import sys, os
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
        tmp = filenames[i]
        tmp = tmp[4:]
        tmp = tmp.replace("_orbit", " (#")
        tmp = tmp.replace(".txt", ")")
        tmp = tmp.replace("_", "Day ")
        name.append(tmp)

    y_label = whatplot + ' ' + y_name + ' Difference ' + ' [1/cm^3]'
    x_label = x_name + ' [deg]'
    graphTitle = 'Temperature Change of ' + name[1] + ' - ' + name[0] + '\nLT ~ ' + str(LT_name)

    if average:
        plt.plot(x_values, y_values, label='Raw DE2')
        plt.plot(x_values_avg, y_values_avg, label='Avg')
    else:
        plt.plot(x_values, y_values, label='Raw DE2')

    if x_lower != 0:
        plt.xlim(x_lower, x_upper)
    if y_lower != 0:
        plt.ylim(y_lower, y_upper)

    #plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(graphTitle)
    plt.legend()
    plt.grid()
    plt.show()

    return 0

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

def temp_difference(temp_1, temp_2, lat_1, lat_2, lon_1, lon_2, LST_1, LST_2, split_1, split_2):
    """
    Takes in two temperature arrays and finds the realtive difference between the values for the same latitude

    Arguments:
        temp_1  :   temperature arr of first file [K]
        temp_2  :   temperature arr of second file [K]
        lat_1   :   latitude arr of first file [deg]
        lat_2   :   latitude arr of second file [deg]
        LST_1   :   local time arr of first file
        LST_2   :   local time arr of second file
        split_1 :   LT split integer for first file
        split_2 :   LT split integer for second file

    Returns:
        lat_diff    :   even though the latitudes will be roughly the same for consecutive files, this returns the average
        temp_diff   :   relative difference between the density values
        LT_name     :   the average LT value for the two orbits. Used later for plotting
    """
    # check to see if we want MLAT instead
    if x_name == "Magnetic Latitude":
        print("Coverting to Geomagentic Coordinates...")
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
    # print('length to the end of 1 =', len(temp_1[start1:]))
    # print('length to the end of 2 =', len(temp_2[start2:]))

    # need to find latitude bin number where they start to overlap
    if lat_1[start1] > lat_2[start2]:
        for i in range(0, length_1):
            if lat_1[start1 + i] <= lat_2[start2]:
                start1 += i
                if len(temp_1[start1:start1+size]) < size:
                    size = len(temp_1[start1:start1+size])                                 # size of overlap shrinks

                break
    else:
        for i in range(0, length_2):
            if lat_2[start2 + i] <= lat_1[start1]:
                start2 += i
                if len(temp_2[start2:start2+size]) < size:
                    size = len(temp_2[start2:start2+size])

                break

    temp_diff = np.zeros(size)                      # difference between densities at same latitude
    lat_overlap = np.zeros(size)
    LT_value = (LST_1[start1] + LST_2[start2])/2


    # print("\nNEW VALUES.....")
    # print('1at1[start1], lat2[start2] = ', lat_1[start1], lat_2[start2])
    # print('size = ', size)
    # print('start1 start2 =', start1, start2)
    # print('length to the end of 1 =', len(temp_1[start1:]))
    # print('length to the end of 2 =', len(temp_2[start2:]))
    #
    # print('lat2[start stop] = ', lat_2[start2: start2+size])
    # print('lat1[start stop] = ', lat_1[start1: start1 + size])

    # now go through and iterate to find the differences
    for i in range(0, (size)):
        #print(i)
        temp_diff[i] = temp_2[start2 + i] - temp_1[start1 + i]

        # DE satellite will have basically the same latitude for the overlapping parts
        lat_overlap[i] = 0.5 * (lat_2[start2 + i] + lat_1[start1 + i])

    if compare:
        plt.plot(lat_1[start1: start1+size], temp_1[start1:start1+size], label='Orbit 7191')
        plt.plot(lat_2[start2: start2+size], temp_2[start2:start2+size], label='Orbit 7153')
        plt.xlabel('Latitude [deg]')
        plt.ylabel('Temperature [K]')
        plt.grid()
        plt.legend()
        plt.show()

    return lat_overlap, temp_diff, LT_value

def delete_bad_values(temp, lat, lon, alt, LST):
    """
    The DE data has some bad temperature values - appearing as zeros in the temperture array after we read everything in
    These must be deleted in order to subtract the two temperature values
    :param temp:    array of temperature values  (with bad values)
    :param lat:
    :param lon:
    :param alt:
    :param LST:
    :return:        arrays for temp, lat, alt, and LST with the bad values deleted
    """

    temp_good, lat_good, lon_good, alt_good, LST_good = [], [], [], [], []

    for i in range(0, len(temp)):
        if temp[i] != 0:
            temp_good.append(temp[i])
            lat_good.append(lat[i])
            lon_good.append(lon[i])
            alt_good.append(alt[i])
            LST_good.append(LST[i])
        else: continue

    return temp_good, lat_good, lon_good, alt_good, LST_good

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
    if not os.path.isfile(mag_filename):
        print("ERROR. Magnetic equator file", mag_filename, "does not exist.")
        sys.exit()

    with open(mag_filename) as csvfile:
        readCSV = csv.reader(csvfile, delimiter=',')
        for line in readCSV:
            mag_lon_eq.append(float(line[0]))
            mag_lat_eq.append(float(line[1]))
    csvfile.close()

    for i in range(0, len(DE_lat)):
        for j in range(0, len(mag_lon_eq)):
            if int(DE_lon[i]) == mag_lon_eq[j]:
                MLAT[i] = DE_lat[i] - mag_lat_eq[j]

    if show_mag:
        plt.scatter(mag_lon_eq, mag_lat_eq, label='Mag. Eq.', marker='.')
        plt.scatter(DE_lon, DE_lat, label='GLAT', marker='.')
        plt.scatter(DE_lon, MLAT, label='MLAT', marker = '.')
        plt.xlabel('Longitude [deg]')
        plt.ylabel('Latitude [deg]')
        plt.legend()
        plt.grid()
        plt.ylim(-180, 180)
        plt.xlim(-100, 100)
        plt.show()


    return MLAT

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


# ----------------------- VARIABLES TO CHANGE -------------------------------------------------------
# ----------------------------------------------------------------------------------------------------
directory = './WATSfiles/'
filenames = ['1982_326_orbit7191.txt', '1982_324_orbit7153.txt']
whatplot = 'Temperature'                # can be Temperature
x_name = "Latitude"                     # can be: 'Latitude' or 'Magnetic Latitude'
y_name = "Temperature"                  # can be: Temperature
average = True                          # do you want to plot the average over it
npts = 25                               # number of points to average array over
show_mag = False                        # show plot of geographic to geomagentic coordinates
show_lon_check = True                  # show where the orbits are longitudinally with respect to the mag eq.
compare = True                          # compare the two tempertaure plots before subtracting from each other
x_lower = 40
x_upper = -40
y_upper = 0
y_lower = 0

# ----------------------------------------------------------------------------------------------------

print_header()

delimiter = '\t'  # file delimiter
hour = 3600000  # ms in an hr
numFiles = len(filenames)  # number of files to read


# check for good paths/ variables
for i in range(0, numFiles):
    if not os.path.isfile(directory + filenames[i]):
        print("ERROR. Read file", filenames[i], "does not exist.")
        sys.exit()
if y_name != 'Temperature':
    print('ERROR: Invalid y axis.')
    sys.exit()
if x_name != 'Latitude' and x_name != 'Magnetic Latitude':
    print('ERROR: Invalid x axis.')
    sys.exit()

# read in file
for j in range(0, numFiles):
    fread = open(directory + filenames[j], 'r')
    print('Opening file', filenames[j])

    # get number of lines and initialize data. Time is in UT [ms], densities are 1/cm^3, errors are %
    line_num = 0
    for line in fread:
        line_num += 1
    fread.close()

    data_lines = int(line_num/2)                                # each data time stamp is split between two liness)

    time = np.zeros(data_lines)
    dens = np.zeros(data_lines)                                 # total neutral density [cm^-3]
    temp = np.zeros(data_lines)                                 # corrected neutral temperature [K]
    alt, lat = np.zeros(data_lines), np.zeros(data_lines)       # altitude and geo lat [km] and [deg]
    lon, LST = np.zeros(data_lines), np.zeros(data_lines)       # geo longitude and local solar time [deg] and [hr]


    # ----------- READ IN DATA-----------------------------------------------------------------
    # -----------------------------------------------------------------------------------------
    fread = open(directory + filenames[j], 'r')

    counter = 0
    k = 0

    for line in fread:
        line = line.lstrip(' ')         # remove leading whitespace if we have it

        if not (counter % 2):           # we are on first line of data split k % 2 = 0
            for i in range(0, 10):
                if i == 1:
                    time[k] = float((line.split(delimiter)[i])) / hour        # UT time in hours
                elif i == 4:
                    dens[k] = float(line.split(delimiter)[i])

                # although temperatures have some bad values (inlcuding zeros, we will take care of these
                # later in post processing)
                elif i == 6 and float(line.split(delimiter)[i]) < 9000 and float(line.split(delimiter)[i]) != -9:
                    temp[k] = float(line.split(delimiter)[i])
                else:
                    continue

        else:                               # second line of split
            for i in range(0, 10):
                if i == 3:
                    alt[k] = float(line.split(delimiter)[i])
                elif i == 4:
                    lat[k] = float(line.split(delimiter)[i])
                elif i == 5:
                    lon[k] = float(line.split(delimiter)[i])
                elif i == 6:
                    LST[k] = float(line.split(delimiter)[i])
                else: continue
            k += 1          # only increment k once per two lines
        counter += 1
    fread.close()

    # -----------------------------------------------------------------------------------------

    # ---------------- GET GOOD VALUES --------------------------------------------------------
    # we are on the first file
    if j == 0:
        temp1, lat1, lon1, alt1, LST1 = delete_bad_values(temp, lat, lon, alt, LST)
        split1 = find_split(LST1)
        print("First file split occurred at:", split1, '\n')


    # we are on the second file
    else:
        temp2, lat2, lon2, alt2, LST2 = delete_bad_values(temp, lat, lon, alt, LST)
        split2 = find_split(LST2)
        print("Second file split occurred at:", split2, '\n')


# end j in range (0, numfiles)


# ----------- POST PROCESSING + PLOTTING ---------------------------------------------------
#  -----------------------------------------------------------------------------------------

x_values, y_values, LT_name = temp_difference(temp1, temp2, lat1, lat2, lon1, lon2, LST1, LST2, split1, split2)


if show_lon_check:
    check_lon(lon1, lat1, lon2, lat2)

if average:
    x_values_avg = avg(x_values, npts)
    y_values_avg = avg(y_values, npts)

plotData(x_values, y_values, filenames, LT_name)






