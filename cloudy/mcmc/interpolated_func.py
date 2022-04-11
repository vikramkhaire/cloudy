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
import sys


def find_nearest(array, value):
    idx = (np.abs(array - value)).argmin()
    return array[idx], idx

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
def get_interp_func_nT(model_path, ions_to_use, Q_uvb, uvb = 'KS18'):
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


def create_missing_models (missing_logT_array, missing_logZ_array, example_file, variable  = 'T'):
    # under linear interpolation
    """
    :param logT_array:
    :param logZ_array:
    :param example_file:
    :param variable: interpolating on T since there are supposed to be more models with T than Z
    :return:
    """

    # getting file info
    logNHI_ref = ((example_file.split('/')[-1]).split('NHI')[-1]).split('_')[0]
    z_ref = ((example_file.split('/')[-1]).split('z_')[-1]).split('.fits')[0]
    model_path = example_file.split('/')[0]
    first_part = (example_file.split('/')[-1]).split('Z')[0]

    data = tab.Table.read(example_file)
    column_names  = data.colnames


    if variable == 'T':
        for Z, logT in zip(missing_logZ_array, missing_logT_array):
            logZ_ref = (Z+4)*100
            find_T_array = []
            list_of_files = glob.glob(model_path + '/' + first_part + 'Z{:.0f}*_z_{}.fits'.format(logZ_ref, z_ref))
            for i in list_of_files:
                read_logT = int(i.split('logT')[1].split('_')[0])
                find_T_array.append(read_logT / 100)

            print(find_T_array)
            
            logT_array = np.sort(np.array(set(find_T_array)))
            valueT1, indexT1 = find_nearest(logT_array, logT)
            logT_array_dummy = np.delete(logT_array, indexT1)
            valueT2, indexT2 =find_nearest(logT_array_dummy, logT)

            # linear interpolation
            # y(x) = y1 + (x-x1)*(y2-y1)/(x2-x1) ==> N(T)  = N1 + [(T-T1)/(T2-T1)] (N2-N1)
            fact= (10**logT - 10**valueT1)/ (10**valueT2 - 10**valueT1)
            # filling
            print('imputing missing values for logZ {} and logT {} and fact {}'.format(Z, logT, fact))

            new_data = tab.Table()
            table1 = tab.Table.read(model_path +'/' + first_part +
                                    'Z{:.0f}_NHI{}_logT{:.0f}_z_{}.fits'.format(logZ_ref, logNHI_ref, valueT1*100, z_ref))
            table2 = tab.Table.read(model_path +'/' + first_part +
                                    'Z{:.0f}_NHI{}_logT{:.0f}_z_{}.fits'.format(logZ_ref, logNHI_ref, valueT2*100, z_ref))
            for col in column_names:
                new_data[col] = table1[col] + fact * (table2[col] - table1[col])

            filename_to_write =  model_path+ '/' + first_part + \
                                 'Z{:.0f}_NHI{}_logT{:.0f}_z_{}_dummy.fits'.format(logZ_ref, logNHI_ref, logT*100, z_ref)
            new_data.write(filename_to_write, overwrite = True)


    return

def get_nZT_array(model_filepath, identifier_redshift, uvb = 'KS18', uvb_Q ='18'):
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

    z_ref = identifier_redshift * 1000000.

    # assuming there is unique cloud (i.e same NHI) at each z value
    #logNHI_ref = identifier_logNH * 100
    #list_of_files = glob.glob(model_filepath + '/*NHI{:.0f}*z_{:.0f}.fits'.format(logNHI_ref, z_ref))
    list_of_files = glob.glob(model_filepath + '/*z_{:.0f}.fits'.format(z_ref))


    logT_array = []
    logZ_array = []
    logNHI_values = []
    for i in list_of_files:
        read_logZ = int(i.split('Z')[1].split('_')[0])
        logZ_array.append(read_logZ / 100 - 4)  # the file naming convention

        read_logT = int(i.split('logT')[1].split('_')[0])
        logT_array.append(read_logT / 100)

        read_logNHI = int(i.split('NHI')[1].split('_')[0])
        logNHI_values.append(read_logNHI)

    # see if NHI values are real
    unique_NHI_values = list(set(logNHI_values))
    if len(unique_NHI_values) !=1:
        print('There are ', len(unique_NHI_values), ' values of NHI. Stopping!')
        sys.exit()
    else:
        logNHI = unique_NHI_values[0]/100


    logZ_array = sorted(list(set(logZ_array)))
    logT_array = sorted(list(set(logT_array)))

    # find if the files are missing
    missing_files =  len(logT_array) * len(logZ_array) - len(list_of_files)
    if missing_files > 0:
        print('Warning: There are {} files missing'.format(missing_files))
        # find the files that are missing
        logNHI_ref = logNHI * 100
        missing_T_array = []
        missing_Z_array = []
        for Z in logZ_array:
            logZ_ref = (Z+4)*100
            for T in logT_array:
                logT_ref = T*100
                file_check = model_filepath + '/try_{}_Q{}_Z{:.0f}_NHI{:.0f}_logT{:.0f}_z_{:.0f}.fits'.format(uvb, uvb_Q, logZ_ref, logNHI_ref, logT_ref, z_ref)
                if not os.path.exists(file_check):
                    print('found missing model for log Z {} and log T {}'.format(Z, T))
                    missing_Z_array.append(Z)
                    missing_T_array.append(T)

        # fill the missing values
        create_missing_models(missing_T_array, missing_Z_array, example_file=list_of_files[0])


    # get nH array
    sample_data = tab.Table.read(list_of_files[0])
    nH_array = list(sample_data['hden'])


    return nH_array, logZ_array, logT_array, logNHI


def get_interp_func_nZT(model_path, ions_to_use, identifier_redshift, uvb = 'KS18', uvb_Q ='18'):

    nH_array, logZ_array, logT_array, identifier_logNH = get_nZT_array(model_path, identifier_redshift)
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
                try:
                    d = tab.Table.read(model)
                except:
                    print('searching for imputed dummy file (logZ, logT)', logT_ref/100 -4, logT_ref/100)
                    model = model_path + '/try_{}_Q{}_Z{:.0f}_NHI{:.0f}_logT{:.0f}_z_{:.0f}_dummy.fits'.\
                        format(uvb, uvb_Q, logZ_ref, logNHI_ref, logT_ref, z_ref)
                    d = tab.Table.read(model)

                data[:, j, k] = d[ion]

        f = RegularGridInterpolator((np.log10(nH_array), logZ_array, logT_array), data, bounds_error = False)
        interpolation_function_list.append(f)

    #print(model)
    return interpolation_function_list

"""
# test 3D interpolation
model_path  = '/home/vikram/data/cloudy/Cloudy'
NHI = 14.77
redshift = 0.12389
ions_to_use = ['Si+2', 'C+2']
flist = get_interp_func_nZT(model_path=model_path, ions_to_use=ions_to_use, identifier_logNH=NHI,
    identifier_redshift=redshift)

"""


