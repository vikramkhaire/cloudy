
import numpy as np
import astropy.table as tab
from scipy.interpolate import interp1d
from scipy.interpolate import interp2d
import matplotlib.pyplot as plt
import emcee
import corner

import os
import glob

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


def get_interp_func_nZT(model_path, ions_to_use, output_filepath, identifier_redshift, identifier_logNH):
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

    #input_file      =   '/scratch/vikram/cloudy_run/abhisek/try_{}_Q{}_Z{:.0f}_NHI{:.0f}_logT{:.0f}_z_{:.0f}.in'.format(uvb, uvb_Q, logZ_ref,logNHI_ref,logT_ref,z_ref)

    logNHI_ref      =   identifier_logNH*100
    z_ref           =   identifier_redshift*1000000.

    list_of_files = glob.glob(output_filepath + '/*NHI{:.0f}*z_{:.0f}.fits'.format(logNHI_ref, z_ref))

    logT_array = []
    logZ_array = []
    for i in list_of_files:
        read_logZ = int(i.split('Z')[1].split('_')[0])
        logZ_array.append(read_logZ/100 -4) # the file naming convention

        read_logT = int(i.split('logT')[1].split('_')[0])
        logT_array.append(read_logT/100)

    logZ_array = sorted(list(set(logZ_array)))
    logT_array = sorted(list(set(logT_array)))

    # get nH array
    sample_data = tab.Table.read(list_of_files[0])
    nH_array  = list(sample_data['hden'])













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

