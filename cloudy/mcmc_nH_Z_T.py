
import numpy as np
import astropy.table as tab
from scipy.interpolate import interp1d
from scipy.interpolate import interp2d
import matplotlib.pyplot as plt
import emcee
import corner


#----data
def get_true_model(model_path, Q= 18, logT = 4,uvb = 'KS18'):
    """
    :param model: The data where Q18 model is stored
    :return: a row of ion column densities at n_H = 1e-4 cm^-2
    """
    NHI=15.08
    model = model_path + '/try_{}_Q{}_logT{:.0f}_NHI{:.0f}.fits'.format(uvb, Q,logT*100,NHI*100)
    data = tab.Table.read(model)
    true_ion_col = data [data['hden'] == 1e-4]
   # print(true_ion_col)

    return true_ion_col

#----model interpolation
def get_interp_func(model_path, ions_to_use, Q_uvb,uvb = 'KS18'):
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


#----for mcmc
def log_likelihood(theta, interp_logf, obs_ion_col, col_err, reference_log_metal = -1.0):
    """
    For a gaussian distributed errors
    :param theta: parameters [nH, Z]
    :param x: data x
    :param y: data y
    :param yerr: data err
    :return:
    """
    lognH, logZ, logT =  theta
    # get metal ion column density for n_H and Z = 0.1
    col = []
    for i in range(len(obs_ion_col)):
        #print('==>', i, lognH, logT)
        #print(interp_logf[i](lognH, logT), i, lognH, logT)
        col_mod = interp_logf[i](lognH, logT)[0]
        col.append(col_mod)

    # scale the column densities by the metallicity Z
    metal_scaling_linear = 10 ** logZ / 10 ** reference_log_metal
    model_col = np.log10(np.array(col) * metal_scaling_linear)

    lnL = -0.5 * np.sum(np.log(2 * np.pi * col_err ** 2) + (obs_ion_col - model_col) ** 2 / col_err ** 2)

    return lnL

def log_prior(theta):
    lognH, logZ, logT =  theta
    # flat prior
    if -4.7 < lognH < -2.0 and -2.0 < logZ < 1.0 and 4.0 < logT < 5.5:   # Better result when ranges are same as grids.
        return 0.0
    return -np.inf

def log_posterior(theta, interp_func, data_col, sigma_col):
    log_p = log_prior(theta) + \
            log_likelihood(theta, interp_logf = interp_func, obs_ion_col = data_col, col_err = sigma_col)

    return log_p


def run_mcmc(model_path, Q_uvb, ions_to_use, data_col, sigma_col, uvb = 'KS18', figname = 'testT.pdf'):
    # run_mcmc(model_Q= model, ions_to_use= ions)
    # ------------------ here is a way to run code
    truths = [-4, -1, 4]  # (lognH, logZ, logT) true values
    number_of_ions = len(ions_to_use)

    print(data_col, sigma_col)

    interp_logf = get_interp_func(model_path = model_path, ions_to_use= ions_to_use, Q_uvb= Q_uvb,uvb=uvb)

    # Here we'll set up the computation. emcee combines multiple "walkers",
    # each of which is its own MCMC chain. The number of trace results will
    # be nwalkers * nsteps

    ndim = 3  # number of parameters in the model
    nwalkers = 50  # number of MCMC walkers
    nsteps = 10000  # number of MCMC steps to take

    # set theta near the maximum likelihood, with
    n_guess = np.random.uniform(-4.68, -2.0, nwalkers)
    z_guess = np.random.uniform(-2.0, 1.0, nwalkers)
    T_guess = np.random.uniform(4.0, 5.5, nwalkers)
    starting_guesses = np.vstack((n_guess, z_guess, T_guess)).T  # initialise at a tiny sphere

    # Here's the function call where all the work happens:
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior, args=(interp_logf, np.array(data_col), np.array(sigma_col)))
    sampler.run_mcmc(starting_guesses, nsteps, progress=True)

    # find out number of steps
    tau = sampler.get_autocorr_time()  # number of steps needed to forget the starting position
    #print(tau)
    thin = int(np.mean(tau) / 2)  # use this number for flattning the sample as done below
    #thin = 100
    flat_samples = sampler.get_chain(discard=thin * 20, thin= 5, flat=True)
    # we are discarding some initial steps roughly 5 times the autocorr_time steps
    # then we thin by about half the autocorrelation time steps for plotting => one does not have to do this step

    labels = [r'log$n_H$', 'log Z', 'log T']
    #uvb_q= int((model_Q.split('try_Q')[-1]).split('.fits')[0])

    #if Q_uvb == true_Q:
      #  fig = corner.corner(flat_samples, labels=labels, quantiles=[0.16, 0.5, 0.84],
     #       show_titles=True, title_kwargs={"fontsize": 12})
    #else:
    fig = corner.corner(flat_samples, labels=labels, quantiles=[0.16, 0.5, 0.84],
            show_titles=True, title_kwargs={"fontsize": 12})

    fig.savefig(figname)
    plt.close()

    for i in range(ndim):
        mcmc = np.percentile(flat_samples[:, i], [16, 50, 84])
        q = np.diff(mcmc)
        print(labels[i], '=', mcmc[1], q[0], q[1])


    return flat_samples, ndim

"""
# testing on KITP simulations
ions_to_use= ['C+','C+2','N+','N+2','N+4','O+5','Si+','Si+2','Si+3']
data_col =[14.81,15.71,14.17,15.12,13.62,14.92,11.88,13.54,13.23]
sigma_col=[0.38,0.16,0.17,0.19,0.37,2.42,0.15,0.05,0.28]

true_Q =18
outpath = '/home/abhisek/Dropbox/Dropbox/uvb_const/Cloudy_results/15.08/temp/outHM12'
model_path  = '/home/abhisek/Dropbox/Dropbox/uvb_const/Cloudy_results/15.08/temp/HM12'

NHI=15.08
uvb = 'HM12'
Q_uvb = 18
name = uvb + '_T{}_NHI_{:.0f}'.format(Q_uvb,NHI*100)
name = uvb + '_Q{}_NHI_{:.0f}'.format(Q_uvb,NHI*100)
figname = outpath + '/' + name + '.pdf'

flat_samples, ndim = run_mcmc(model_path=model_path, data_col=data_col, sigma_col=sigma_col, Q_uvb=Q_uvb, ions_to_use=ions_to_use,figname=figname, uvb=uvb)

"""
