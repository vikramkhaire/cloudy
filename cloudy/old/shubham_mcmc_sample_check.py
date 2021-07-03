#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 11:03:00 2021

@author: jarvis-astro
"""


import numpy as np
import astropy.table as tab
from scipy.interpolate import interp1d
from scipy.interpolate import interp2d
import matplotlib.pyplot as plt
import matplotlib as mpl
import astropy.table as tab
import matplotlib.pyplot as plt


#----data
def get_true_model(model_Q, Q= 18):
    """
    :param model: The data where Q18 model is stored
    :return: a row of ion column densities at n_H = 1e-4 cm^-2
    """
    model = model_Q.split('_Q')[0] + '_Q{}.fits'.format(Q)
    data = tab.Table.read(model)
    true_ion_col = data [data['hden'] == 1e-4]
   # print(true_ion_col)

    return true_ion_col

#----model interpolation

def get_interp_func(model_path, ions_to_use, Q_uvb = 18, uvb = 'KS18'):
    logZ = np.around(np.arange(-2, 1, 0.2), decimals = 2) # hardcoded
    #get nH array
    logZ_try = -1
    model_try = model_path + '/try_{}_Q{}_Z{:.0f}.fits'.format(uvb, Q_uvb, (logZ_try+4)*100)
    model = tab.Table.read(model_try)
    lognH = np.log10(np.array(model['hden']))

    interpolation_function_list = []
#    number_of_ions = len(ions_to_use)
    for ion in ions_to_use:
        z = np.zeros((len(lognH), len(logZ)))
        for i in range(len(logZ)):
            model = model_path + '/try_{}_Q{}_Z{:.0f}.fits'.format(uvb, Q_uvb, (logZ[i]+4)*100)
            model_q = tab.Table.read(model)
            sorted_model = model_q[ion]
            model_touple = ()
#            for j in range(len(ions_to_use)):
#                model_touple += (sorted_model[ions_to_use[j]],)
            sorted_model[sorted_model == 0 ] = 1e-15 # for avoiding log10 (0) error
            z [:, i] = np.log10(sorted_model) #--- for log - log interpolation
#        logf = interp2d(lognH, logZ, np.log10(model_touple), fill_value='extrapolate')
        logf = interp2d(lognH, logZ, z.T)

        interpolation_function_list.append(logf)
#    return logf
    return interpolation_function_list




def sample_plot(npy_file, interp_logf, true_ions, true_col, true_sigma, figname ='test.png', samples = 400,
                reference_log_metal = -1):

    # setting the figure
    font = {'family': 'serif', 'weight': 'normal', 'size': 14}
    plt.rc('font', **font)
    mpl.rcParams['axes.linewidth'] = 1.5

    figure_size = [7, 6]
    fig, ax = plt.subplots(1, 1, figsize=(figure_size[0], figure_size[1]))

    flat_samples = np.load(npy_file)
    inds = np.random.randint(len(flat_samples), size=samples)
    

#    model_col  = np.array(col)
    
    for ind in inds:
        sample = flat_samples[ind]
        lognH, logZ = sample
#        print(lognH)
#        print(logZ)
        col = []
        for i in range(len(true_col)):
#            print(i)
            col_mod = 10 ** interp_logf[i](lognH,logZ)[0]
            col.append(col_mod)
        col_m = np.array(col)
        # scale the column densities by the metallicity Z
#        metal_scaling_linear = 10 ** logZ / 10 ** reference_log_metal
#        model_col = col_m * metal_scaling_linear
        ax.plot(true_ions, np.log10(col_m),  "C1", alpha=0.2)

    ax.errorbar(true_ions, true_col, yerr = true_sigma, color = 'k', elinewidth =2, marker = '.',  markersize = 15,
                linestyle = '', zorder =1000)

    ax.set_xlabel('Ion')
    ax.set_ylabel(r'N$_{\rm ion}$ (cm$^{-2}$)')
    # ax.set_xscale('log')
#    ax.set_yscale('log')

    # deco
    ax.tick_params(direction='in', length=5, width=1.5)
    ax.tick_params(direction='in', which='minor', length=3.5, width=1.5)
    # ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    # decorating the plot
    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(1.5)

    fig.savefig(figname, bbox_inches='tight')




def get_mcmc_sample(model_path, ions_to_use, data_col, sigma_col, npyfile,  true_Q =18, figname = 'test.png', same_error = False):
    # run_mcmc(model_Q= model, ions_to_use= ions)
    # ------------------ here is a way to run code
    number_of_ions = len(ions_to_use)

#    data_col_all = get_true_model(model_Q, Q=true_Q)
    # converting astropy table row to a list
#    data_col = []
#    for name in ions_to_use:
#        data_col.append(data_col_all[name][0])


    np.random.seed(0)
#    if same_error:
#        sigma_col = 0.2 * np.ones(number_of_ions)
#    else:
#        sigma_col = np.random.uniform(0.1, 0.3, number_of_ions)

#    print(np.log10(data_col), sigma_col)

    interp_logf = get_interp_func(model_path = model_path, ions_to_use = ions_to_use)

#    interp_logf = get_interp_func(model, list(ions_to_use))
    sample_plot(npy_file= npyfile, interp_logf= interp_logf, true_ions= ions_to_use, true_col= data_col,
                true_sigma= sigma_col, figname = figname)

    return

#uvb_q= 18
Q_uvb = 18
uvb = 'KS18'
#model = []
model_path  = '/home/jarvis-astro/cloudy_run/metal_NH18_85'
logZ = np.around(np.arange(-2, 1, 0.2), decimals = 2)
    #get nH array
#logZ_try = -1
#model_try = model_path + '/try_KS18_Q{}_Z{:.0f}.fits'.format(Q_uvb, (logZ_try+4)*100)
#model_try_1 = tab.Table.read(model_try)
#lognH = np.log10(np.array(model_try_1['hden']))
#model = []
#for i in range(len(logZ)):
#    model_Q = model_path + '/try_KS18_Q{}_Z{:.0f}.fits'.format(Q_uvb, (logZ[i]+4)*100)
#    model.append(model_Q)
ions_to_use = ['Si+', 'C+', 'Si+2', 'C+2']
data_col = np.array([[17.47, 17.67, 16.78, 16.51], [17.44, 16.49, 17.47, 16.61], [17.87, 16.87, 16.43, 17.85], [17.08, 16.09, 17.07, 17.91], [16.38, 16.25, 16.92, 16.38], [17.62, 17.21, 16.7,  16.54], [16.11, 16.95, 17.91, 16.48], [16.27, 16.8,  17.27, 17.15], [17.47, 16.86, 17.79, 17.16], [16.44, 16.58, 16.29, 17.8 ], [16.98, 17.32, 16.81, 17.57], [16.95, 17.07, 17.46, 17.83], [17.56, 17.97, 17.08, 17.74], [17.28, 17.27, 16.52, 17.94], [17.18, 17.93, 16.14, 16.81], [17.67, 16.97, 17.95, 17.82], [17.15, 16.44, 16.58, 16.08], [16.25, 16.3,  16.43, 17.22], [16.58,17.,   17.16, 16.95], [17.56, 16.5,  17.8,  17.78]])
sigma_col = np.array([[0.43, 0.3,  0.34, 0.21], [0.22, 0.36, 0.23, 0.21], [0.31, 0.45, 0.48, 0.29], [0.46, 0.2,  0.21, 0.22], [0.27, 0.45, 0.31, 0.48], [0.42, 0.33, 0.5,  0.42], [0.27, 0.46, 0.3,  0.29], [0.44, 0.43, 0.47, 0.3 ], [0.36, 0.4,  0.26, 0.4 ], [0.41, 0.37, 0.23, 0.4 ], [0.42, 0.45, 0.22, 0.4 ], [0.38, 0.48, 0.22, 0.28], [0.37, 0.26, 0.46, 0.33], [0.32, 0.34, 0.31, 0.39], [0.22, 0.25, 0.42, 0.24], [0.5,  0.21, 0.39, 0.43], [0.46, 0.24, 0.49, 0.42], [0.34, 0.36, 0.33, 0.33], [0.25, 0.34, 0.27, 0.22], [0.32, 0.43, 0.34, 0.23]])
#ions_to_use = np.array(ions_to_use)[np.argsort(data_col)]
#print(ions_to_use, ': sorted ions_to_use')
#data_col = data_col[np.argsort(data_col)]
#print(data_col, ': sorted data_col')
#sigma_col = sigma_col[np.argsort(data_col)]
#print(sigma_col, ': sorted sigma_col')
true_Q =18

npyfile = '/home/jarvis-astro/cloudy_run/figures/without_limits/KS18_Q{}.npy'.format(Q_uvb)

#get_mcmc_sample(model_path, ions_to_use, data_col, sigma_col, npyfile, figname= 'test18_85_sample.png')

for j in range(0, len(random_sigma_SiII)):
get_mcmc_sample(model_path, ions_to_use, data_col, sigma_col, npyfile, figname= 'test18_85_sample.png')
