from cloudy.mcmc.mcmc_nH_Z_T import run_mcmc
from cloudy.mcmc.interpolated_func import get_interp_func_nT_aiswarya as interpol

model_path = '/home/vikram/iist_projects/aiswarya'
ions_to_use =
data_col =
sigma_col =

interp_logf = interpol(model_path = model_path, ions_to_use= ions_to_use, Q_uvb= Q_uvb,uvb=uvb)


