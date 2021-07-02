#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  1 19:46:44 2021

@author: jarvis-astro
"""

import numpy as np
import random
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

    interpolation_function_list = []
    for ion in ions_to_use:
        z = np.zeros((len(lognH), len(logZ)))
        for i in range(len(logZ)):
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
    lognH, logZ=  theta
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

def log_prior(theta):
    lognH, logZ =  theta
    # flat prior
    if -6 < lognH < -2 and -3 < logZ < 1 :
        return 0.0
    return -np.inf

def log_posterior(theta, interp_func, data_col, sigma_col):
    log_p = log_prior(theta) + \
            log_likelihood(theta, interp_logf = interp_func, obs_ion_col = data_col, col_err = sigma_col)

    return log_p


def run_mcmc(model_path, Q_uvb, ions_to_use, data_col=None, sigma_col=None, true_Q =18, uvb = 'KS18', figname = 'testT.pdf', same_error = False):
    # run_mcmc(model_Q= model, ions_to_use= ions)
    # ------------------ here is a way to run code
    truths = [-4, -1]  # (lognH, logZ, logT) true values
    number_of_ions = len(ions_to_use)

#    data_col = np.array([13.83, 15.38, 14.35, 14.61, 14.47, 14.27])
#    sigma_col = np.array([0.32, 0.51, 0.04, 0.67, 0.76, 0.12])

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

    # Here we'll set up the computation. emcee combines multiple "walkers",
    # each of which is its own MCMC chain. The number of trace results will
    # be nwalkers * nsteps

    ndim = 2  # number of parameters in the model
    nwalkers = 50  # number of MCMC walkers
    nsteps = 5000  # number of MCMC steps to take

    # set theta near the maximum likelihood, with
    n_guess = np.random.uniform(-5, -3, nwalkers)
    z_guess = np.random.uniform(-2, 0, nwalkers)
    np.random.seed(1)
    starting_guesses = np.vstack((n_guess, z_guess)).T  # initialise at a tiny sphere

    # Here's the function call where all the work happens:
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior, args=(interp_logf, data_col, sigma_col))
    sampler.run_mcmc(starting_guesses, nsteps, progress=True)

    # find out number of steps
    tau = sampler.get_autocorr_time()  # number of steps needed to forget the starting position
    #print(tau)
    thin = int(np.mean(tau) / 2)  # use this number for flattning the sample as done below
    #thin = 100
    flat_samples = sampler.get_chain(discard=thin * 20, thin= 5, flat=True)
    # we are discarding some initial steps roughly 5 times the autocorr_time steps
    # then we thin by about half the autocorrelation time steps for plotting => one does not have to do this step

    labels = ['log nH (cm$\mathregular{^{-3}}$)', 'log Z ($Z_{\odot}$)']
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






# generating random 20 observations for all 4 ions
random_col_SiII = []
for i in range(0,20):
    #n = random.randint(16,18)
    n = round(random.uniform(16,18), 2)
    random_col_SiII.append(n)
#print(random_col_SiII)

random_col_CII = []
for i in range(0,20):
    #n = random.randint(16,18)
    n = round(random.uniform(16,18), 2)
    random_col_CII.append(n)
#print(random_col_CII)

random_col_SiIII = []
for i in range(0,20):
    #n = random.randint(16,18)
    n = round(random.uniform(16,18), 2)
    random_col_SiIII.append(n)
#print(random_col_SiIII)

random_col_CIII = []
for i in range(0,20):
    #n = random.randint(16,18)
    n = round(random.uniform(16,18), 2)
    random_col_CIII.append(n)
#print(random_col_CIII)

random_sigma_SiII = []
for i in range(0,20):
    #n = random.randint(16,18)
    n = round(random.uniform(0.2,0.5), 2)
    random_sigma_SiII.append(n)
#print(random_sigma_SiII)

random_sigma_CII = []
for i in range(0,20):
    #n = random.randint(16,18)
    n = round(random.uniform(0.2,0.5), 2)
    random_sigma_CII.append(n)
#print(random_sigma_CII)

random_sigma_SiIII = []
for i in range(0,20):
    #n = random.randint(16,18)
    n = round(random.uniform(0.2,0.5), 2)
    random_sigma_SiIII.append(n)
#print(random_sigma_SiIII)

random_sigma_CIII = []
for i in range(0,20):
    #n = random.randint(16,18)
    n = round(random.uniform(0.2,0.5), 2)
    random_sigma_CIII.append(n)
#print(random_sigma_CIII)


col = np.array([random_col_SiII, random_col_CII, random_col_SiIII, random_col_CIII])
#print(col)
#print(col[:,0])
col_1 = []
for i in range(len(random_col_SiII)):
#    col = []
    col_2 = col[:,i]
    col_1.append(col_2)
#    print(col)
data_col = np.array(col_1)
print('data_col: ', data_col)

sigma = np.array([random_sigma_SiII, random_sigma_CII, random_sigma_SiIII, random_sigma_CIII])
#print(sigma)
#print(sigma[:,0])
sigma_1 = []
for i in range(len(random_sigma_SiII)):
#    sigma = []
    sigma_2 = sigma[:,i]
    sigma_1.append(sigma_2)
#    print(col)
sigma_col = np.array(sigma_1)
print('sigma_col: ', sigma_col)


ions_to_use= ['Si+', 'C+', 'Si+2', 'C+2']
#data_col = np.array([16.37, 17.82, 16.96, 17.16])
#sigma_col = np.array([0.57, 0.46, 1, 1])
true_Q =18

outpath = '/home/jarvis-astro/cloudy_run/figures'
model_path  = '/home/jarvis-astro/cloudy_run/metal_NH18_85'
#outfile = outpath + '/metal_NH18_85_2D' + '.fits'

uvb_array = ['KS18']
Q_array= [18]

#out_tab =  tab.Table()
for j in range(0, len(random_sigma_SiII)):
    for uvb, q in zip(uvb_array, Q_array):
        name = uvb + '_Q{}'.format(q)
        #    for j in range(0, len(random_sigma_SiII) - 1):
        figname = outpath + '/' + name + '_' + str(j) +  '.pdf'

        flat_samples, ndim = run_mcmc(model_path= model_path, Q_uvb=q, ions_to_use=ions_to_use,
                                      data_col=data_col[j], sigma_col=sigma_col[j], true_Q=true_Q,
                                      figname=figname, uvb = uvb)
        # to efficiently save numpy array
        save_file_name = outpath + '/' + name + '_' + str(j)
        np.save(save_file_name, flat_samples)

        out =[[q]]
        for i in range(ndim):
            mcmc = np.percentile(flat_samples[:, i], [16, 50, 84])
            q = np.diff(mcmc)
            out.append([mcmc[1]])
            out.append([q[0]])
            out.append([q[1]])
        out_tab =  tab.Table()
        print(out)
        t = tab.Table(out, names = ('Q', 'nH', 'n16', 'n84', 'Z', 'Z16', 'Z84'))
        out_tab = tab.vstack((out_tab, t))


    outfile = outpath + '/metal_NH18_85_2D' + '_' + str(j) + '.fits'
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
