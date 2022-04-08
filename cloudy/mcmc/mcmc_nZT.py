import numpy as np
import matplotlib.pyplot as plt
import emcee
import corner
import math
import multiprocessing as mp

#----for mcmc
def log_likelihood(theta, interp_logf, data_col, sigma_col, reference_log_metal = -1.0, Z_scaling = True):
    """
    :param theta:  parameters [nH, Z and T]
    :param interp_logf: interpolation function
    :param data_col: data column density array
    :param sigma_col: errors on data_col
    :param reference_log_metal: to scale metallicity if Z_scaling below is True
    :param Z_scaling:
    :return: log of likelihood
    """

    """
    For a gaussian distributed errors
    :param theta: parameters [nH, Z and T]
    :param x: data x
    :param y: data y
    :param yerr: data err
    :return:
    """
    lognH, logZ, logT =  theta

    if Z_scaling:
        # get metal ion column density for n_H and Z = 0.1
        col = []
        for i in range(len(data_col)):
            col_mod = interp_logf[i](lognH, logT)[0]
            col.append(col_mod)

        # scale the column densities by the metallicity Z
        metal_scaling_linear = 10 ** logZ / 10 ** reference_log_metal
        model_col = np.log10(np.clip(col, 1e-10, 1e22) * metal_scaling_linear)

    else:
        col = []
        for i in range(len(data_col)):
            col_mod = interp_logf[i]([lognH, logZ, logT])[0]
            col.append(col_mod)

        # replacing nan values with column of 1e-10
        #col = [1e-10 if math.isnan(x) else x for x in col]
        col= np.nan_to_num(col, nan= 1e-10, copy = False)

        model_col = np.log10(np.clip(col, 1e-10, 1e22)) # clip to prevent log 0 error
        #print(model_col, [lognH, logZ, logT]) #....for debug

    lnL = -0.5 * np.sum(np.log(2 * np.pi * sigma_col ** 2) + (data_col - model_col) ** 2 / sigma_col ** 2)

    return lnL

def log_prior(theta):
    lognH, logZ, logT =  theta
    # flat prior
    if -6 < lognH < 0 and -3.0 < logZ < 1 and 4 < logT < 6:   # Better result when ranges are same as grids.
        return 0.0
    return -np.inf

def log_posterior(theta, interp_func, data_col, sigma_col, Z_scaling = True):
    """
    :param theta:
    :param interp_func:
    :param data_col:
    :param sigma_col:
    :param Z_scaling:
    :return:
    """
    log_p = log_prior(theta) + \
            log_likelihood(theta, interp_func, data_col, sigma_col, Z_scaling = Z_scaling)

    return log_p


def run_mcmc(data_col, sigma_col, interp_logf, nwalkers = 50, nsteps =5000, ndim =3,  Z_scaling = True, figname = 'testT.pdf', parallel = False):
    """
    :param data_col: array with column densities of different metals
    :param sigma_col: corresponding errors
    :param interp_logf: the interpolation function list to give column density of same ions (in same order) from cloudy models
    :param figname: figure name for the corner plot
    :return: flattened chain and number of dimensions
    """

    # ------------------ here is a way to run code
    print(data_col, sigma_col)


    # Here we'll set up the computation. emcee combines multiple "walkers",
    # each of which is its own MCMC chain. The number of trace results will
    # be nwalkers * nsteps

    # set theta near the maximum likelihood, with
    n_guess = np.random.uniform(-4, -3, nwalkers)

    # initializing for two densities
    #n_guess = np.concatenate((np.random.uniform(-3.6, -3.5, int(nwalkers/2)), np.random.uniform(-1, 2, int(nwalkers/2))))

    z_guess = np.random.uniform(-2, 0, nwalkers)
    T_guess = np.random.uniform(4.1, 5.5, nwalkers)
    starting_guesses = np.vstack((n_guess, z_guess, T_guess)).T  # initialise at a tiny sphere


    if parallel:
        # somehow it makes it very slow :)
        with mp.Pool() as pool:
            # Here's the function call where all the work happens:
            sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior,
                args=(interp_logf, np.array(data_col), np.array(sigma_col), Z_scaling), pool=pool)
            sampler.run_mcmc(starting_guesses, nsteps, progress=True)
    else:
        # Here's the function call where all the work happens:
        sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior,
            args=(interp_logf, np.array(data_col), np.array(sigma_col), Z_scaling))
        sampler.run_mcmc(starting_guesses, nsteps, progress=True)


    # find out number of steps
    tau = sampler.get_autocorr_time()  # number of steps needed to forget the starting position
    #print(tau)
    thin = int(np.mean(tau) / 2)  # use this number for flattning the sample as done below
    #thin = 100
    flat_samples = sampler.get_chain(discard=thin * 10, thin= 2, flat=True)
    # we are discarding some initial steps roughly 5 times the autocorr_time steps
    # then we thin by 5 (or if thin =  thin; half the autocorrelation time steps) for plotting => one does not have to do this step

    labels = [r'log $n_H$', 'log Z', 'log T']
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

