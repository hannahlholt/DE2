'''
NEEDS PYTHON 2.7 TO COMPILE AND RUN WITH MSIS

NACS_massDens.py reads in a specified text file containing UT, gas number densities, and density errors for a certain orbit of the DE-2 satellite.
This program plots the average mean molecular weight as a function of latitude for the specified orbit file.
'''


from __future__ import print_function, division
import sys, os
import numpy as np
import matplotlib.pyplot as plt


filenames = ['orbit7164.txt']
m_O, m_N2, m_He, m_N, m_Ar = 16, 28, 4, 14, 40


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
    data_lines -= 2                  # don't care about first two lines (header)

    O_dens = np.zeros(data_lines)
    N2_dens = np.zeros(data_lines)
    He_dens = np.zeros(data_lines)
    N_dens = np.zeros(data_lines)
    Ar_dens = np.zeros(data_lines)
    alt, lat = np.zeros(data_lines), np.zeros(data_lines)

    # now read in the data
    fread = open(filenames[j], 'r')

    # skip header
    fread.next()
    fread.next()

    k = 0
    # read in all the values from the file, all densities are number densities [1/cm^3]
    for line in fread:
        for i in range(0, 16):
            if i == 0:
                continue
            elif i == 1:
                O_dens[k] = float(line.split(delimiter)[i])
            elif i == 2:
                continue
            elif i == 3:
                N2_dens[k] = float(line.split(delimiter)[i])
            elif i == 4:
                continue
            elif i == 5:
                He_dens[k] = float(line.split(delimiter)[i])
            elif i == 6:
                continue
            elif i == 7:
                N_dens[k] = float(line.split(delimiter)[i])
            elif i == 8:
                continue
            elif i == 9:
                Ar_dens[k] = float(line.split(delimiter)[i])
            elif i == 10:
                continue
            elif i == 11:
                continue
            elif i == 12:
                alt[k] = float(line.split(delimiter)[i])
            elif i == 13:
                lat[k] = float(line.split(delimiter)[i])
            else:
                continue
        k += 1
    fread.close()

    # find weighted average molecular mass

    avg_mass = np.zeros(data_lines)
    for i in range(0, data_lines):
        sum_tot = O_dens[i] + N2_dens[i] + He_dens[i] + N_dens[i] + Ar_dens[i]
        avg_mass[i] = (m_O*O_dens[i] + m_N2*N2_dens[i] + m_He*He_dens[i] + m_N*N_dens[i] + m_Ar*Ar_dens[i]) / sum_tot

    # now plot average mean molecular mass (amu) vs. lat

    name = filenames[j].replace("orbit", "")
    name = name.replace(".txt", "")

    fig, ax1 = plt.subplots()
    ax1.plot(lat, avg_mass, color='r', label='Mean Molecular Mass')
    ax1.set_ylabel('Avg. Molecular Mass [amu]')
    ax1.set_xlabel('Latitude [deg]')
    plt.legend(loc=1)
    plt.grid()

    ax2 = ax1.twinx()
    ax2.plot(lat, alt, color='b', linestyle='--', label='DE-2 Altitude')
    ax2.set_ylabel('Altitude [km]')
    plt.legend(loc=4)
    plt.title('DE2 Orbit ' + name + ' Average Molecular Mass')
    plt.show()



