#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 16:03:08 2021

@author: jarvis-astro
"""


import numpy as np
import astropy.table as tab
from scipy.interpolate import interp1d
from scipy.interpolate import interp2d
import matplotlib.pyplot as plt
import emcee
import corner


#----data
def get_true_model(model_path, Q= 18, logZ = -1, uvb = 'KS18'):
    """
    :param model: The data where Q18 model is stored
    :return: a row of ion column densities at n_H = 1e-4 cm^-2
    """
    model = model_path + '/try_{}_Q{}_Z{:.0f}.fits'.format(uvb, Q, (logZ+4)*100)
    data = tab.Table.read(model)
    true_ion_col = data [data['hden'] == 1e-4]
   # print(true_ion_col)

    return true_ion_col

#----model interpolation
def get_interp_func(model_path, ions_to_use, Q_uvb, uvb = 'KS18'):
    logZ = np.around(np.arange(-2, 1, 0.2), decimals = 2) # hardcoded
    #get nH array
    logZ_try = -1
    model_try = model_path + '/try_{}_Q{}_Z{:.0f}.fits'.format(uvb, Q_uvb, (logZ_try+4)*100)
    model = tab.Table.read(model_try)
    lognH = np.log10(np.array(model['hden']))
#    logSII = np.log10(np.array(model['S+']))

    interpolation_function_list = []
    for ion in ions_to_use:
        z = np.zeros((len(lognH), len(logZ)))
        for i in range(len(logZ)):
#            if logSII.any() < 14.4:
                model = model_path + '/try_{}_Q{}_Z{:.0f}.fits'.format(uvb, Q_uvb, (logZ[i]+4)*100)
                d = tab.Table.read(model)
                d[ion][d[ion] == 0 ] = 1e-15 # for avoiding log10 (0) error
                z [:, i] = np.log10(d[ion]) #--- for log - log interpolation
        f = interp2d(lognH, logZ, z.T)
        interpolation_function_list.append(f)

    return interpolation_function_list


#----for mcmc
def log_likelihood(theta, interp_logf, obs_ion_col, col_err):
    """
    For a gaussian distributed errors
    :param theta: parameters [nH, Z]
    :param x: data x
    :param y: data y
    :param yerr: data err
    :return:
    """
    lognH, logZ =  theta
    # get metal ion column density for n_H and Z = 0.1
    col = []
    for i in range(len(obs_ion_col)):
        #print('==>', i, lognH, logT)
        #print(interp_logf[i](lognH, logT), i, lognH, logT)
        col_mod = interp_logf[i](lognH, logZ)[0]
        col.append(col_mod)

    model_col  = np.array(col)

    lnL = -0.5 * np.sum(np.log(2 * np.pi * col_err ** 2) + (obs_ion_col - model_col) ** 2 / col_err ** 2)

    return lnL

#def log_prior(theta, interp_logf, obs_ion_col):
def log_prior(theta, interp_logf_for_lower_limit, data_col_lower_limit):
    lognH, logZ =  theta

    # finding column densities of ions used as limits in prior
    lower_limit_col_mod = []
    for i in range(len(interp_logf_for_lower_limit)):
        logcol_lower = interp_logf_for_lower_limit[i](lognH, logZ)[0]
        lower_limit_col_mod.append(logcol_lower)
    lower_limit_col = np.array(lower_limit_col_mod)
    


    # todo for shubham
    """
    if_condition below should take all elements one by one and compare with the limit values
    if its only on element SII
    then if_condition is "limit_col[0] < limit_col[0]"
    """
#    for i in range(len(interp_logf_for_limits)):
#        condition = limit_col[i] < data_col_limits[i]
#        condition_array = condition.all()
    condition_lower = lower_limit_col > data_col_lower_limit
    condition_array_lower = condition_lower.all()
    

    

    # flat prior
    if -6 < lognH < -2 and -3 < logZ < 1 and condition_array_lower == True:
            return 0.0
    # flat prior
#    if -6 < lognH < -2 and -3 < logZ < 1 and  if_condition:
#            return 0.0
    return -np.inf

def log_posterior(theta, interp_func, interp_func_for_lower_limit, data_col_lower_limit, data_col, sigma_col):
    log_p = log_prior(theta, interp_logf_for_lower_limit = interp_func_for_lower_limit, data_col_lower_limit=data_col_lower_limit) + \
            log_likelihood(theta, interp_logf = interp_func, obs_ion_col = data_col, col_err = sigma_col)

    return log_p


def run_mcmc(model_path, Q_uvb, ions_to_use, ions_to_use_for_lower_limit = None, data_col=None, sigma_col=None, data_col_lower_limit=None, true_Q =18, uvb = 'KS18', figname = 'testT.pdf', same_error = False):
    # run_mcmc(model_Q= model, ions_to_use= ions)
    # ------------------ here is a way to run code


#    data_col = np.array([13.83, 15.38, 14.35, 14.61, 14.47, 14.27])
#    sigma_col = np.array([0.32, 0.51, 0.04, 0.67, 0.76, 0.12])
#    limit_col = np.array([])

#    print(data_col, sigma_col)
    
    
#    ionname = list()
#    data_col = list()
#    sigma_col = list()
#    f = open('/home/jarvis-astro/cloudy_run/ion list.txt','r')
#    for line in f:
#        line = line.strip('"')
#        columns = line.split('"')
#        ionname_val = str(columns[0])
#        ionname.append(ionname_val)
#        data_col_val = float(columns[1])
#        data_col.append(data_col_val)
#        sigma_col_val = float(columns[2])
#        sigma_col.append(sigma_col_val)
#    f.close() 
#    ionname = np.asarray(ionname,dtype=str)
#    data_col = np.asarray(data_col,dtype=float)
#    sigma_col = np.asarray(sigma_col,dtype=float)
#    print('Ion Name: ', ionname)
    print('Observed Column Density: ', data_col)
    print('Error in Observed Column Density: ', sigma_col)


    interp_logf = get_interp_func(model_path = model_path, ions_to_use = ions_to_use, Q_uvb = Q_uvb, uvb = uvb)
    interp_logf_for_lower_limit = get_interp_func(model_path = model_path, ions_to_use = ions_to_use_for_lower_limit,
                                             Q_uvb = Q_uvb, uvb = uvb)


    # Here we'll set up the computation. emcee combines multiple "walkers",
    # each of which is its own MCMC chain. The number of trace results will
    # be nwalkers * nsteps

    ndim = 2  # number of parameters in the model
    nwalkers = 50  # number of MCMC walkers
    nsteps = 5000  # number of MCMC steps to take

    # set theta near the maximum likelihood, with
    n_guess = np.random.uniform(-5, -3, nwalkers)
    z_guess = np.random.uniform(-2, 0, nwalkers)
    np.random.seed(0)
    starting_guesses = np.vstack((n_guess, z_guess)).T  # initialise at a tiny sphere

    # Here's the function call where all the work happens:
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior, args=(interp_logf, interp_logf_for_lower_limit, data_col_lower_limit, data_col, sigma_col))
    sampler.run_mcmc(starting_guesses, nsteps, progress=True)

    # find out number of steps
    tau = sampler.get_autocorr_time()  # number of steps needed to forget the starting position
    #print(tau)
    thin = int(np.mean(tau) / 2)  # use this number for flattning the sample as done below
    #thin = 100
    flat_samples = sampler.get_chain(discard=thin * 20, thin= 5, flat=True)
    # we are discarding some initial steps roughly 5 times the autocorr_time steps
    # then we thin by about half the autocorrelation time steps for plotting => one does not have to do this step

    labels = ['log nH (in cm$\mathregular{^{-3}}$)', 'log Z (in $Z_{\odot}$)']
    #uvb_q= int((model_Q.split('try_Q')[-1]).split('.fits')[0])

    if Q_uvb == true_Q:
        fig = corner.corner(flat_samples, labels=labels, quantiles=[0.16, 0.5, 0.84],
            show_titles=True, title_kwargs={"fontsize": 12})
    else:
        fig = corner.corner(flat_samples, labels=labels, quantiles=[0.16, 0.5, 0.84],
            show_titles=True, title_kwargs={"fontsize": 12})

    fig.savefig(figname)

    plt.close()

    for i in range(ndim):
        mcmc = np.percentile(flat_samples[:, i], [16, 50, 84])
        q = np.diff(mcmc)
        print(labels[i], '=', mcmc[1], q[0], q[1])


    return flat_samples, ndim






ions_to_use= ['Si+', 'N+2', 'C+']
data_col = np.array([14.47, 14.5, 14.79])
sigma_col = np.array([0.14, 0.02, 0.02])
ions_to_use_for_upper_limit = ['N+4', 'S+', 'S+2', 'S+3', 'S+5', 'Fe+']
data_col_upper_limit = np.array([13.4, 14.4, 14.0, 13.7, 12.8, 11.4])
ions_to_use_for_lower_limit = ['C+2', 'Si+2']
data_col_lower_limit = np.array([14.5, 14.0])

true_Q =18

outpath = '/home/jarvis-astro/cloudy_run/figures'
model_path  = '/home/jarvis-astro/cloudy_run/metal_NH18'
outfile = outpath + '/metal_NH18_2D.fits'

uvb_array = ['KS18']
Q_array= [18]

out_tab =  tab.Table()
for uvb, q in zip(uvb_array, Q_array):
    name =uvb + '_Q{}'.format(q)
    figname = outpath + '/' + name + '.pdf'

    flat_samples, ndim = run_mcmc(model_path= model_path, Q_uvb=q, ions_to_use=ions_to_use, sigma_col=sigma_col,
                                  ions_to_use_for_lower_limit=ions_to_use_for_lower_limit, data_col=data_col,
                                  data_col_lower_limit=data_col_lower_limit, true_Q=true_Q, figname=figname, uvb = uvb)
    # to efficiently save numpy array
    save_file_name = outpath + '/' + name
    np.save(save_file_name, flat_samples)

    out =[[q]]
    for i in range(ndim):
        mcmc = np.percentile(flat_samples[:, i], [16, 50, 84])
        q = np.diff(mcmc)
        out.append([mcmc[1]])
        out.append([q[0]])
        out.append([q[1]])

    print(out)
    t = tab.Table(out, names = ('Q', 'nH', 'n16', 'n84', 'Z', 'Z16', 'Z84'))
    out_tab = tab.vstack((out_tab, t))



uvb_column = ['Q18']
out_tab.add_column(uvb_column, name = 'uvb')

out_tab.write(outfile, overwrite = True)


"""
ions_to_use= ['C+3', 'N+3', 'Si+3', 'O+5', 'C+2']
true_Q =18
outpath = '/home/vikram/cloudy_run/figures/2DLLS'
model_path  = '/home/vikram/cloudy_run/metal_NH19'
outfile = outpath + '/NH19_metal_2D.fits'
uvb_array = ['KS18', 'KS18', 'KS18', 'KS18', 'KS18', 'KS18', 'KS18', 'P19', 'FG20', 'HM12']
Q_array= [14, 15, 16, 17, 18, 19, 20, 18, 18, 18]
out_tab =  tab.Table()
for uvb, q in zip(uvb_array, Q_array):
    name =uvb + '_Q{}'.format(q)
    figname = outpath + '/' + name + '.pdf'
    flat_samples, ndim = run_mcmc(model_path= model_path, Q_uvb=q, ions_to_use=ions_to_use, true_Q=true_Q,
        figname=figname, uvb = uvb)
    # to efficiently save numpy array
    save_file_name = outpath + '/' + name
    np.save(save_file_name, flat_samples)
    out =[[q]]
    for i in range(ndim):
        mcmc = np.percentile(flat_samples[:, i], [16, 50, 84])
        q = np.diff(mcmc)
        out.append([mcmc[1]])
        out.append([q[0]])
        out.append([q[1]])
    print(out)
    t = tab.Table(out, names = ('Q', 'nH', 'n16', 'n84', 'Z', 'Z16', 'Z84'))
    out_tab = tab.vstack((out_tab, t))
uvb_column = ['Q14', 'Q15', 'Q16', 'Q17', 'Q18', 'Q19', 'Q20', 'P19', 'FG20', 'HM12']
out_tab.add_column(uvb_column, name = 'uvb')
out_tab.write(outfile, overwrite = True)
"""