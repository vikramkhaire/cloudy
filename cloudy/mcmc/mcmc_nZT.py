
import numpy as np
import astropy.table as tab
from scipy.interpolate import interp1d
from scipy.interpolate import interp2d
import matplotlib.pyplot as plt
import emcee
import corner


#----for mcmc
def log_likelihood(theta, interp_logf, obs_ion_col, col_err, reference_log_metal = -1.0, Z_scaling = True):
    """
    For a gaussian distributed errors
    :param theta: parameters [nH, Z]
    :param x: data x
    :param y: data y
    :param yerr: data err
    :return:
    """
    lognH, logZ, logT =  theta

    if Z_scaling:
        # get metal ion column density for n_H and Z = 0.1
        col = []
        for i in range(len(obs_ion_col)):
            col_mod = interp_logf[i](lognH, logT)[0]
            col.append(col_mod)

        # scale the column densities by the metallicity Z
        metal_scaling_linear = 10 ** logZ / 10 ** reference_log_metal
        model_col = np.log10(np.array(col) * metal_scaling_linear)

    else:
        model_col = []
        for i in range(len(obs_ion_col)):
            col_mod = interp_logf[i](lognH, logT, logZ)[0]
            model_col.append(col_mod)

    lnL = -0.5 * np.sum(np.log(2 * np.pi * col_err ** 2) + (obs_ion_col - model_col) ** 2 / col_err ** 2)

    return lnL

def log_prior(theta):
    lognH, logZ, logT =  theta
    # flat prior
    if -4.7 < lognH < -2.0 and -2.0 < logZ < 1.0 and 4.0 < logT < 5.5:   # Better result when ranges are same as grids.
        return 0.0
    return -np.inf

def log_posterior(theta, interp_func, data_col, sigma_col, Z_scaling = True):
    log_p = log_prior(theta) + \
            log_likelihood(theta, interp_logf = interp_func, obs_ion_col = data_col, col_err = sigma_col, Z_scaling = Z_scaling)

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
