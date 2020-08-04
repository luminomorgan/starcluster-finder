#!/usr/bin/python
"""

  This code loads a list of spectra and S-PLUS filters and convolves the spectra
  with the respective S-PLUS filters, using the Simpson's numerical intregration. 
  The output are the synthetic magnitudes. 


                                MVCD - 18/08/2018

"""
import sys
import time
import numpy as np
import os
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import simps
import pandas as pd

def calc_lpivot(l, T):

    """

    It calculates the lambda pivot, based on the following website

    http://www.astro.ljmu.ac.uk/~ikb/research/mags-fluxes/

    Arguments:

    Transmission curves array -> T
    lambdas array -> l
    
    Output: 

    l_pivot
    
                                         MVCD - 12/03/2018
    """

    idx = np.where(l > 0.)[0]
    return np.sqrt(simps(T[idx], l[idx]) / simps(T[idx] / l[idx]**2, l[idx]))

def load_jplus_filter(infile):
    """
    Load filter transmission

    """

    data = [line.strip() for line in open(infile)]
    l_filter = []
    s_filter = []
    for line in data:
        p = line.split()
        if p[0] != '#':
            l_filter.append(float(p[0]))
            s_filter.append(float(p[1]))

    l_filter = np.array(l_filter)
    s_filter = np.array(s_filter)
    return l_filter,s_filter

def load_SPLUS_set_filters(infile_list, path_filters, dl):

    """
      This routine loads the SPLUS filter set and outputs into 
      a matrix. 

    """
    # Load list of filters

    file_filters = [line.strip() for line in open(path_filters + infile_list)]

    
    nline = np.zeros(len(file_filters)) # number of lambdas of each filter
    for i in range(len(file_filters)):
        l, s = load_jplus_filter(path_filters + file_filters[i])
        l2 = np.arange(float(int(min(l))), float(int(max(l))), dl) 
        nline[i] = len(l2)

    s_filter_matrix = np.zeros(len(file_filters) * int(max(nline))).reshape(len(file_filters), int(max(nline)))
    l_filter_matrix = np.zeros(len(file_filters) * int(max(nline))).reshape(len(file_filters), int(max(nline)))

    # Filling matrices with interpolated QEs in 1\AA

    for i in range(len(file_filters)):
        l, s = load_jplus_filter(path_filters + file_filters[i])
        l_new = np.arange(float(int(min(l))), float(int(max(l))), dl) 
        f_new = np.interp(l_new, l, s, left = 0.0, right = 0.0)
        s_filter_matrix[i][:len(l_new)] = f_new
        l_filter_matrix[i][:len(l_new)] = l_new

    return l_filter_matrix, s_filter_matrix, file_filters

def read_ascii_table(infile, header_string, columns, type_variable):

    columns = np.array(columns)
    type_variable = np.array(type_variable)

    # Arrays columns and type_variables have different lengths

    if len(columns) != len(type_variable): 
        print("[read_ascii_table] len(columns) != len(type_variables)")
        exit()

    # Load table
    f = open(infile,'r')
    data = f.readlines()
    f.close()

    matrix = []
    for i in range(len(columns)):
        matrix.append([])

    n0 = -1
    for line in data:
        p = line.split()
        nline = -1
        n0 += 1
        #print n0, p[0]
        if p[0][0] != header_string: 
            for i in range(len(columns)):
                nline += 1
                if type_variable[i] == 0: # string
                    matrix[nline].append(str(p[columns[i]]))
                if type_variable[i] == 1: # float
                    matrix[nline].append(float(p[columns[i]]))
                if type_variable[i] == 2: # integer
                    matrix[nline].append(int(p[columns[i]]))
                if type_variable[i] > 2 or type_variable[i] < 0:
                    print("[read_ascii_table] type_variables -> bad values!!")
                    exit()
    return matrix

##############################################

flag_interp = 0 # (0/1) = (no interp / interp)

path_filters = "./filter/"
list_filters = "SPLUS_filters.list"

#list_spec = 'NGSLnewsorted.list'
#path_spec = './NGSLsorted_new/'

list_spec = "spec.list"
path_spec = "./spec/"

## Wavelength range (\AA)

l_min = 500.
l_max = 100000.
dl = 1. # wavelength resolution
l = np.arange(l_min, l_max + 1., 1.)

## Normalization range

l_norm_i = 6100.
l_norm_f = 6200.

## Lowest mag uncertainty

lowest_mag_uncertainty = 0.03

## Lowest reduced chi2 

chi2_red_min = 3e-2

## Minimum spectral SN

SN_min = 1. 

## Output file

outfile = 'convolved_mags_' + time.strftime("%d%m%Y") + '.dat'

## Spectral window for SN ratio estimate

lambda_SN_i = 4010. 
lambda_SN_f = 4060.
  
## lightspeed in AA/s

c_AAs = 2.99792458e18     

# Plotting spec + synthetic mags

flag_plot = 0. # (0/1) = (DO NOT plot/ plot)

########################################################################

if __name__ == '__main__':

    # Make the lambda array

    l_out = np.arange(l_min, l_max + 1., dl)

    print(">> Reading input list<<")
    print()
    print("list=", list_spec)
    print()
    #
    f = open(list_spec,'r')
    data = f.readlines()
    f.close()
    #
    path_sed = []
    file_spec = []
    for line in data:
        p = line.split()
        if p[0] != '#':
            file_spec.append(str(p[0]))

    file_spec = np.array(file_spec)
    print("Number of spectra = ", len(file_spec))
    
    # Loading SPLUS filters

    l_filters, qe_filters, file_filters_SPLUS = load_SPLUS_set_filters(list_filters, path_filters, dl)
    n_filters = len(l_filters[:,1])
    print("'n_filters = ", n_filters)
    print(file_filters_SPLUS)

    # Calculate lambda pivot

    lambda_pivot = np.zeros(len(file_filters_SPLUS))
    for i in range(len(file_filters_SPLUS)):
        lambda_pivot[i] = calc_lpivot(l_filters[i], qe_filters[i])

    # Declare synthetic mags matrix

    mag_AB = np.zeros(len(l_filters) * len(file_spec)).reshape(len(file_spec), len(l_filters))
    mag_AB_scipy = np.zeros(len(l_filters) * len(file_spec)).reshape(len(file_spec), len(l_filters))

    for ii in range(len(file_spec)):

        print(ii, path_spec + file_spec[ii])

        # Loading spectrum

        columns = [0, 1]
        type_variables = [1, 1]
        data = read_ascii_table(path_spec + file_spec[ii], '#', columns, type_variables)

        # Define arrays of lambda and flux

        l0 = np.array(data[0], dtype = float)
        f0 = np.array(data[1], dtype = float)
        f0 *= 1e-17 # 1e-17 erg/s/A/cm2

        # Interpolating the spectrum at the red part

        f = np.interp(l, l0, f0)

        # Declare output array

        for j in range(len(l_filters)):

            idx = np.where(l_filters[j] > 0.)[0] # valid lambdas

            # SCIPY convolution

            filt_int  = np.interp(l, l_filters[j][idx], qe_filters[j][idx])              #Interpolate to common wavelength axis
            I1        = simps(f * filt_int * l, l)                     #Denominator
            I2        = simps(filt_int / l, l)                     #Numerator
            fnu       = (I1 / I2) / c_AAs                          #Average flux density
            mag_AB_scipy[ii, j] = -2.5 * np.log10(fnu) - 48.6              #AB magnitude

        flux_AB_scipy = 10**(-0.4 * (mag_AB_scipy[ii, :] + 48.6))
        flux_AB_scipy *= c_AAs / (lambda_pivot ** 2)
        flux_AB_scipy /= 1e-17 # in 1e-17 units

        # Plotting...

        if flag_plot == 1:

            fig = plt.figure()
            ax = fig.add_subplot(111)
            plt.plot(l, f / 1e-17, color = 'red', alpha = 0.5)
            plt.plot(lambda_pivot, flux_AB_scipy, 'o', color = 'green', alpha = 0.7)
            plt.xlabel('Wavelength(A)')
            plt.ylabel('Flux (1e-17 erg/s/A/cm2)')
            plt.xlim([2500., 11500.])
            plt.savefig(file_spec[ii] + '.png', bbox_inches = 'tight')
            plt.clf()
            plt.cla()
            #plt.show()

    # Output 

    file_filters_SPLUS[0] = 'file_spec ' + file_filters_SPLUS[0]
    df = pd.DataFrame(data = mag_AB_scipy, index = file_spec, columns = file_filters_SPLUS)
    df.to_csv(outfile, sep='\t')

    exit()
