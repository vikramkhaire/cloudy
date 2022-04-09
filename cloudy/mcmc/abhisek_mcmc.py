#------------------
# for reducing  numpy threads
import os
#os.environ["NUMEXPR_NUM_THREADS"] = "1"
#os.environ["OMP_NUM_THREADS"] = "1"
#os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
#os.environ["OPENBLAS_NUM_THREADS"] = "1"
#os.environ["MKL_NUM_THREADS"] = "1"
#-----------------
import multiprocessing as mp
import numpy as np
import matplotlib.pyplot as plt
import glob
from cloudy.mcmc.interpolated_func import  get_interp_func_nZT
from cloudy.mcmc.mcmc_nZT import run_mcmc
import astropy.table as tab
import sys
import time

def run_parallel(redshift):
    # file paths:
    model_path = '/home/vikram/data/cloudy/Cloudy'
    observation_path = '/home/vikram/data/cloudy/observations'
    output_path = '/home/vikram/data/cloudy/output'

    # find ions_to_use for absorber at z= redshift
    z_ref = redshift*1e6
    list_of_files = glob.glob(observation_path + '/z_{:.0f}*.dat'.format(z_ref))
    if len(list_of_files) !=1:
        print('more observations at single redshift', list_of_files)
        sys.exit()
    else:
        observed_file = list_of_files[0]

    data = tab.Table.read(observed_file, format= 'ascii')
    ions_to_use =data['ions']
    print(list(ions_to_use))
    data_col_log = data['N']
    sigma_col_log = data['eN']

    # setting figure name
    figname = output_path + '/' + (observed_file.split('/')[-1]).split('.dat')[0]
    print('Output File name and Path', figname)

    final_filename = figname + '.pdf'
    if not os.path.exists(final_filename):
        print('file does not exits', final_filename)

        # get interpolated functions
        try:
            func_list = get_interp_func_nZT(model_path=model_path, ions_to_use=ions_to_use,
                identifier_redshift=redshift)

            flat_samples, ndim = run_mcmc(data_col=data_col_log, sigma_col=sigma_col_log, interp_logf=func_list,
                figname=figname + '.pdf', Z_scaling=False, parallel=False)

            # file to save mcmc chain
            save_file_name = figname
            np.save(save_file_name, flat_samples)
            print('*************** saved output for', redshift)
        except Exception as e:
            print(e)
            print('exception at ', redshift)
    else:
        print('Calculation exists, see', final_filename)



    return


z_array= np.array([0.004409, 0.005602, 0.043318, 0.059285, 0.060158,
       0.063275, 0.077493, 0.077701, 0.078068, 0.094864, 0.098787,
       0.113918, 0.123596, 0.12389 , 0.124783, 0.135467, 0.140754,
       0.146789, 0.161068, 0.166588, 0.170062, 0.187731, 0.292317,
       0.310529, 0.349368, 0.360841, 0.386094, 0.39346 , 0.42188 ,
       0.423919, 0.424307, 0.44678 ])




#pool = mp.Pool(processes=3)
#pool.imap_unordered(run_parallel, z_array)


if __name__ == '__main__':
    with mp.Pool(6) as p:
        p.map(run_parallel, z_array)
"""
if __name__=='__main__':
    starttime = time.time()
    pool = mp.Pool(processes=3)
    pool.map(run_parallel, z_array)
    pool.close()
    print('That took {} seconds'.format(time.time() - starttime))

"""

#with mp.Pool(processes=4) as pool:
#    for i in pool.imap_unordered(run_parallel, z_array):
#        print(i)
# print same numbers in arbitrary order
#for i in pool.imap_unordered(f, range(10)):
#    print(i)
#results = [pool.apply_async(run_parallel, args=(redshift,)) for redshift in z_array]
#output = [p.get() for p in results]


"""
# testing mcmc code
model_path  = '/home/vikram/data/cloudy/Cloudy'
output_path = '/home/vikram/data/cloudy/output_test'

# find ions_to_use for absorber at z= redshift

redshift = 0.12389
figname = output_path + '/'+ 'z_{}'.format(redshift)


ions_to_use = ['O+5', 'N+4', 'C+3', 'Si+2', 'Si+3' ]

data_col_log = [14.16, 12.98, 13.82, 12.35, 12.81]
sigma_col_log = [0.17, 0.2, 0.1, 0.11, 0.21]

# get interpolated functions
func_list = get_interp_func_nZT(model_path=model_path, ions_to_use=ions_to_use, identifier_redshift=redshift)

flat_samples, ndim = run_mcmc(data_col=data_col_log, sigma_col=sigma_col_log, interp_logf=func_list,
    figname=figname + '.pdf', Z_scaling=False)

# file to save mcmc chain
save_file_name = figname
np.save(save_file_name, flat_samples)
"""

