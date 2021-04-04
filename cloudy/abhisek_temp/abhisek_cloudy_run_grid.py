#------------------
# for reducing  numpy threads
import os
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
#-----------------
from cloudy_run import write_input
from cloudy_run import run
from cloudy_run import store_table
from cloudy_run import cloudy_params_defaults
import multiprocessing as mp
import numpy as np
from write_uvb_in_cloudy_format import write_uvb_in_cloudy_format
import astropy.table as tab


def run_parallel(logNHI,logZ, uvb_Q, uvb):
    # for vikram
    cloudy_path = '/home/abhisek/Soft/c17.02'
    fname = (logZ+4)*100
    fname1=logNHI*100
    input_File = '/home/abhisek/Desktop/PHD_final/work/with_vikram/cloudy_test/try_{}_Q{}_Z{:.0f}_NHI{}.in'.format(uvb, uvb_Q, fname,fname1)
    print(uvb, 'Q=', uvb_Q, 'Z=', logZ, 'logNHI=', logNHI)

    # write input file and run cloudy
    ions, params = cloudy_params_defaults(uvb = uvb, uvb_Q=uvb_Q, log_hden=[-4.5, -1.5, 0.02], stop_NHI = logNHI, T = None, metal = logZ,
                                          sequential = True,CMB='CMB',z=0.25)
    write_input(input_File, *ions, **params)
    run(cloudy_path=cloudy_path, input_file=input_File)

    # write output tables
    output_filename = input_File.split('.in')[0] +  '.spC'
    fits_filename = input_File.split('.in')[0] + '.fits'
    store_table(ions=ions, output_file=output_filename, fits_filename=fits_filename)

    return

# runnning in parallel
#uvb_array= [14, 15, 16, 17, 18, 19, 20]
logZ_array = np.around(np.arange(-2.5, .55, 0.05), decimals = 2)
logNHI_array = np.around(np.arange(15, 16.5, 0.05), decimals = 2)
uvb = ['KS18', 'HM12']
uvb_Q = [ 18]

logZ = []
logNHI = []
uvb_models =[]
the_Q_values = []
for background in uvb:
    if background == 'KS18':
        for q in uvb_Q:
            for metal in logZ_array:
                for NHI in logNHI_array:
                    uvb_models.append(background)
                    the_Q_values.append(q)
                    logZ.append(metal)
                    logNHI.append(NHI)
    else:
        q = 18
        for metal in logZ_array:
                for NHI in logNHI_array:
                    uvb_models.append(background)
                    the_Q_values.append(q)
                    logZ.append(metal)
                    logNHI.append(NHI)

#-----write uvb fg and hm in cloudy format first
#path = '/mnt/quasar2/vikram/cloudy_run/metal_NH19_new'

#kwagrs = {'uvb' : 'P19', 'z' : 0.2}
#uvb_files(path, **kwagrs)

#kwagrs = {'uvb' : 'FG20', 'z' : 0.2}
#uvb_files(path, **kwagrs)


pool = mp.Pool(processes=90)
results = [pool.apply_async(run_parallel, args=(NHI, Z, Q, mod,)) for  NHI, Z, Q, mod in zip(logNHI, logZ, the_Q_values, uvb_models)]
output = [p.get() for p in results]




