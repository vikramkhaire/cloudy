import numpy as np
import astropy.table as tab
from scipy.interpolate import interp1d
from scipy.interpolate import interp2d
import matplotlib.pyplot as plt
import emcee
import corner


#----model interpolation
def get_interp_func(model_path, ions_to_use, Q_uvb, uvb = 'KS18', logZ_try = -1):
    logZ = np.around(np.arange(-3, 1, 0.05), decimals = 2) # hardcoded
    #get nH array
    model_try = model_path + '/try_{}_Q{}_Z{:.0f}.fits'.format(uvb, Q_uvb, (logZ_try+4)*100)
    model = tab.Table.read(model_try)
    lognH = np.log10(np.array(model['hden']))
    print(len(lognH) )

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

    print('Observed Column Density: ', data_col)
    print('Error in Observed Column Density: ', sigma_col)


    interp_logf = get_interp_func(model_path = model_path, ions_to_use = ions_to_use, Q_uvb = Q_uvb, uvb = uvb)

    # Here we'll set up the computation. emcee combines multiple "walkers",
    # each of which is its own MCMC chain. The number of trace results will
    # be nwalkers * nsteps

    ndim = 2  # number of parameters in the model
    nwalkers = 100  # number of MCMC walkers
    nsteps = 10000  # number of MCMC steps to take

    # set theta near the maximum likelihood, with
    n_guess = np.random.uniform(-4, -2, nwalkers)
    z_guess = np.random.uniform(-2, 0, nwalkers)
    np.random.seed(1)
    starting_guesses = np.vstack((n_guess, z_guess)).T  # initialise at a tiny sphere

    # Here's the function call where all the work happens:-----------------
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior, args=(interp_logf, data_col, sigma_col))
    sampler.run_mcmc(starting_guesses, nsteps, progress=True)

    # find out number of steps
    tau = sampler.get_autocorr_time()  # number of steps needed to forget the starting position
    #print(tau)
    thin = int(np.mean(tau) / 2)  # use this number for flattening the sample as done below
    #thin = 100
    flat_samples = sampler.get_chain(discard=thin * 20, thin= thin, flat=True)
    # we are discarding some initial steps roughly 5 times the autocorr_time steps
    # then we thin by about half the autocorrelation time steps for plotting => one does not have to do this step

    #--------------Plotting-------------
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



