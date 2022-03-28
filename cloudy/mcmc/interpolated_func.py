# list of custom interpolation functions for various projects.
import numpy as np
import astropy.table as tab
from scipy.interpolate import interp2d
from scipy.interpolate import RegularGridInterpolator
import matplotlib.pyplot as plt
import emcee
import corner
import os
import glob

#--- model interpolation for Aiswarya's project (23 March 2022)
def get_interp_func_nT_aiswarya(model_path, ions_to_use, Q_uvb, uvb = 'KS18'):

    logT = np.around(np.arange(3.8, 6, 0.01), decimals = 2)
    #get nH array
    logT_try = 4
    model_try = model_path + '/try_Q{}_T{:.0f}.fits'.format(Q_uvb, logT_try*100)
    model = tab.Table.read(model_try)
    lognH = np.log10(np.array(model['hden']))

    interpolation_function_list = []
    for ion in ions_to_use:
        z = np.zeros((len(lognH), len(logT)))
        for i in range(len(logT)):
            model = model_path + '/try_Q{}_T{:.0f}.fits'.format(Q_uvb, logT[i]*100)
            d = tab.Table.read(model)
            z [:, i] = d[ion]
        f = interp2d(lognH, logT, z.T)
        interpolation_function_list.append(f)

    return interpolation_function_list


#----model interpolation
def get_interp_func_nT(model_path, ions_to_use, Q_uvb,uvb = 'KS18'):
    logT = np.around(np.arange(4.0, 5.50, 0.02), decimals = 2)
    #get nH array
    logT_try = 4
    model_try = model_path + '/try_{}_Q{}_logT{:.0f}_NHI{:.0f}.fits'.format(uvb, Q_uvb,logT_try*100,NHI*100)
    model = tab.Table.read(model_try)
    lognH = np.log10(np.array(model['hden']))

    interpolation_function_list = []
    for ion in ions_to_use:
        z = np.zeros((len(lognH), len(logT)))
        for i in range(len(logT)):
            model = model_path + '/try_{}_Q{}_logT{:.0f}_NHI{:.0f}.fits'.format(uvb, Q_uvb,(logT[i])*100,NHI*100)
            d = tab.Table.read(model)
            z [:, i] = d[ion]
        f = interp2d(lognH, logT, z.T)
        interpolation_function_list.append(f)

    return interpolation_function_list


def get_nZT_array(model_filepath, identifier_redshift, identifier_logNH):
    # hard coded things .....
    # trial files
    """
    # the files are named with following code
        logZ_ref        =   (logZ+4)*100
        logNHI_ref      =   logNHI*100
        logT_ref        =   (logT)*100
        z_ref           =   (z_re*1000000.)
        input_file      =   '/scratch/vikram/cloudy_run/abhisek/try_{}_Q{}_Z{:.0f}_NHI{:.0f}_logT{:.0f}_z_{:.0f}.in'.format(uvb, uvb_Q, logZ_ref,logNHI_ref,logT_ref,z_ref)
    """

    # input_file      =   '/scratch/vikram/cloudy_run/abhisek/try_{}_Q{}_Z{:.0f}_NHI{:.0f}_logT{:.0f}_z_{:.0f}.in'.format(uvb, uvb_Q, logZ_ref,logNHI_ref,logT_ref,z_ref)

    logNHI_ref = identifier_logNH * 100
    z_ref = identifier_redshift * 1000000.

    list_of_files = glob.glob(model_filepath + '/*NHI{:.0f}*z_{:.0f}.fits'.format(logNHI_ref, z_ref))

    logT_array = []
    logZ_array = []
    for i in list_of_files:
        read_logZ = int(i.split('Z')[1].split('_')[0])
        logZ_array.append(read_logZ / 100 - 4)  # the file naming convention

        read_logT = int(i.split('logT')[1].split('_')[0])
        logT_array.append(read_logT / 100)

    logZ_array = sorted(list(set(logZ_array)))
    logT_array = sorted(list(set(logT_array)))

    # get nH array
    sample_data = tab.Table.read(list_of_files[0])
    nH_array = list(sample_data['hden'])

    return nH_array, logZ_array, logT_array


def get_interp_func_nZT(model_path, ions_to_use, identifier_redshift, identifier_logNH, uvb = 'KS18', uvb_Q ='18'):

    nH_array, logZ_array, logT_array = get_nZT_array(model_path, identifier_redshift, identifier_logNH)
    print(logZ_array)
    # hardcoded for filenames
    logNHI_ref = identifier_logNH * 100
    z_ref = (identifier_redshift * 1000000.)

    """
    we need data in this format 
    data [i, j, k]  = fuction (x[i], y[j] , z[k])
    to interpolate using
    my_interpolating_function = RegularGridInterpolator((x, y, z), data)
    """

    #creating data for each ion
    interpolation_function_list = []
    for ion in ions_to_use:
        data = np.zeros((len(nH_array), len(logZ_array), len(logT_array)))
        for j in range(len(logZ_array)):
            for k in range(len(logT_array)):
                logZ_ref = (logZ_array[j]+4)*100
                logT_ref = logT_array[k]*100
                model = model_path + '/try_{}_Q{}_Z{:.0f}_NHI{:.0f}_logT{:.0f}_z_{:.0f}.fits'.format(uvb, uvb_Q, logZ_ref,logNHI_ref,logT_ref,z_ref)
                d = tab.Table.read(model)
                data[:, j, k] = d[ion]

        f = RegularGridInterpolator((np.log10(nH_array), logZ_array, logT_array), data)
        interpolation_function_list.append(f)

    print(model)
    return interpolation_function_list

"""
# test 3D interpolation
model_path  = '/home/vikram/data/cloudy/Cloudy'
NHI = 15.11
redshift = 0.002377
ions_to_use = ['Si+2', 'C+2']
flist = get_interp_func_nZT(model_path=model_path, ions_to_use=ions_to_use, identifier_logNH=NHI,
    identifier_redshift=redshift)
"""

