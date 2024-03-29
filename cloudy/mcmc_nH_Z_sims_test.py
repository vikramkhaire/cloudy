import numpy as np
import astropy.table as tab
from scipy.interpolate import interp1d
from scipy.interpolate import interp2d
import matplotlib.pyplot as plt
import emcee
import corner


#----model interpolation
def get_interp_func(model_path, ions_to_use, Q_uvb, uvb = 'KS18', logZ_try = 0):
    logZ = np.around(np.arange(-0.5, 0.5, 0.05), decimals = 2) # hardcoded
    #get nH array
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
    n_guess = np.random.uniform(-4.5, -3, nwalkers)
    z_guess = np.random.uniform(-0.5, 0.5, nwalkers)
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
    flat_samples = sampler.get_chain(discard=thin * 20, thin= 5, flat=True)
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



"""
##Availabel column densities,'H','H+','He','He+','He+2','C','C+','C+2','C+3','C+4','C+5','N','N+','N+2',
'N+3','N+4','O','O+','O+2','O+3','O+4','O+5','O+6','O+7','S','S+','S+2','S+3','S+4','S+5','Si','Si+',
'Si+2','Si+3','Si+4','Mg','Mg+','Mg+2','Ne','Ne+','Ne+2','Ne+3','Ne+4','Ne+5','Ne+6','Ne+7','Ne+8','Fe','Fe+','Fe+2','Na','Na+','Na+2
"""


ions_to_use = ['C+','C+2','N+','N+2','N+4','O','O+5','Si+','Si+2','Si+3']
data_col=np.array([14.81,15.71,14.17,15.12,13.62,11.55,14.92,11.88,13.54,13.23])
sigma_col=np.array([0.38,0.16,0.17,0.19,0.37,0.31,2.42,0.15,0.05,0.28])

#ions_to_use = ['C+','C+2','N+']
#data_col=np.array([14.81,15.71,14.17])
#sigma_col=np.array([0.38,0.16,0.17])
#ions_to_use= ['Si+', 'N+2', 'C+']
#data_col = np.array([14.47,14.5,14.79])
#sigma_col = np.array([0.14,0.02,0.02])
true_Q =18

outpath = '/home/vikram/Dropbox/uvb_const/Cloudy_run_simulated/output_vik'
model_path  = '/home/vikram/Dropbox/uvb_const/Cloudy_run_simulated/HI_15.08'
outfile = outpath + '/metal_NH18_2D.fits'

uvb_array = ['KS18','HM12']
Q_array= [18, 18]

out_tab =  tab.Table()
for uvb, q in zip(uvb_array, Q_array):
    name =uvb + '_Q{}'.format(q)
    figname = outpath + '/' + name + '.pdf'

    flat_samples, ndim = run_mcmc(model_path= model_path, Q_uvb=q, ions_to_use=ions_to_use,
                                  data_col=data_col, sigma_col=sigma_col, true_Q=true_Q,
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



uvb_column = ['Q18']
out_tab.add_column(uvb_column, name = 'uvb')

out_tab.write(outfile, overwrite = True)

