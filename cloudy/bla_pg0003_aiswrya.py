from cloudy.mcmc.mcmc_nZT import run_mcmc
from cloudy.mcmc.interpolated_func import get_interp_func_nT_aiswarya as interpol
import numpy as np

model_path = '/home/vikram/iist_projects/aiswarya/output/'

"""
# observations
C+	    13.41	0.03
C+2	    13.62	0.09
Si+	    12.46	0.05
Si+2	12.94	0.1
O+5	    13.79	0.04
"""

ions_to_use = ['C+', 'C+2', 'Si+', 'Si+2', 'O+5']
data_col_log = [13.41, 13.62, 12.46, 12.94, 13.79]
sigma_col_log = [0.03, 0.09, 0.05, 0.10, 0.04]

Q_uvb= '18'
uvb = 'KS18'
figname = 'all_ions_full'

interp_logf = interpol(model_path = model_path, ions_to_use= ions_to_use, Q_uvb= Q_uvb, uvb=uvb)

flat_samples, ndim  = run_mcmc(data_col = data_col_log, sigma_col =sigma_col_log, interp_logf = interp_logf, figname = figname + '.pdf')

# file to save mcmc chain
save_file_name = figname
np.save(save_file_name, flat_samples)
