import numpy as np
import matplotlib.pyplot as plt
import emcee
import corner
from cloudy.mcmc.interpolated_func import  get_interp_func_nZT

z_array= np.array([0.004409, 0.005602, 0.042275, 0.043318, 0.059285, 0.060158,
       0.063275, 0.077493, 0.077701, 0.078068, 0.094864, 0.098787,
       0.113918, 0.123596, 0.12389 , 0.124783, 0.135467, 0.140754,
       0.146789, 0.161068, 0.166588, 0.170062, 0.187731, 0.292317,
       0.310529, 0.349368, 0.360841, 0.386094, 0.39346 , 0.42188 ,
       0.423919, 0.424307, 0.44678 ])

model_path  = '/home/vikram/data/cloudy/Cloudy'
output_path = ''

for redshift in z_array:
    figname = output_path + 'z_{}'

    # find ions_to_use for absorber at z= redshift

    ions_to_use =

    data_col_log = []
    sigma_col_log = []

    # get interpolated functions
    func_list = get_interp_func_nZT(model_path=model_path, ions_to_use=ions_to_use, identifier_redshift=redshift)


    flat_samples, ndim = run_mcmc(data_col=data_col_log, sigma_col=sigma_col_log, interp_logf=func_list,
        figname=figname + '.pdf', Z_scaling = False)

    # file to save mcmc chain
    save_file_name = figname
    np.save(save_file_name, flat_samples)


