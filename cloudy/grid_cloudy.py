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


def run_parallel(logNHI,logZ, uvb_Q, uvb,logT,z_re):
    cloudy_path     =   '/home/abhisek/soft/c17.02'
    #cloudy_path     =   '/home/abhisek/Soft/c17.02'
    logZ_ref        =   (logZ+4)*100
    logNHI_ref      =   logNHI*100
    logT_ref        =   (logT)*100
    z_ref           =   (z_re*1000000.)
    input_File      =   '/home/abhisek/soft/cld/run/try_{}_Q{}_Z{:.0f}_NHI{:.0f}_logT{:.0f}_z_{:.0f}.in'.format(uvb, uvb_Q, logZ_ref,logNHI_ref,logT_ref,z_ref)
    #input_File      =   '/home/abhisek/mega/NHI_check/test_NHI/T_4_CMB/cloudy/T_4_CMB/run/try_{}_Q{}_Z{:.0f}_NHI{:.0f}_logT{:.0f}_z{:.0f}.in'.format(uvb, uvb_Q, logZ_ref,logNHI_ref,logT_ref,z_ref)
    print(uvb, 'Q=', uvb_Q, 'Z=', logZ, 'logNHI=', logNHI,'T=', logT,'z=',z_re)

    # write input file and run cloudy
    ions, params    =   cloudy_params_defaults(uvb = uvb, uvb_Q=uvb_Q, log_hden=[-6.0, -2.0, 0.02], stop_NHI = logNHI, metal = logZ,T = 10**logT,sequential = True,z=z_re,CMB='CMB')
    write_input(input_File, *ions, **params)
    run(cloudy_path =   cloudy_path, input_file=input_File)

    # write output tables
    output_filename =   input_File.split('.in')[0] +  '.spC'
    fits_filename   =   input_File.split('.in')[0] + '.fits'
    store_table(ions=ions, output_file=output_filename, fits_filename=fits_filename)

    return

# runnning in parallel
#uvb_array= [14, 15, 16, 17, 18, 19, 20]
#input from absorber tables [to be updated]

logNHI_array        =   np.array([15.11, 15.02, 16.05, 14.48, 14.69, 15.25, 19.23, 15.61, 15.89,
                                    15.28, 14.39, 14.86, 15.41, 15.1 , 14.33, 14.46, 15.84, 14.86,
                                    15.24, 15.42, 13.82, 14.5 , 14.77, 18.  , 16.19, 16.07, 15.16,
                                    16.25, 14.78, 13.99, 18.  , 15.81, 14.63, 14.77, 14.56, 14.54,
                                    15.17, 14.6 , 14.83, 15.23, 14.98, 13.57, 15.75, 13.68, 13.75,
                                    13.83, 15.94, 15.65, 15.53, 14.37, 15.16, 15.16, 15.21, 15.8 ,
                                    16.61, 15.96, 17.21, 16.13, 16.38, 15.08, 15.04, 13.46, 15.35,
                                    16.21, 15.96, 15.49, 15.68, 16.04, 14.46, 15.09, 13.25, 13.7 ,
                                    15.5 , 14.63, 14.35, 15.2 , 15.49, 15.66, 16.48, 14.88, 14.73,
                                    14.69, 15.94, 14.26, 17.54, 16.26, 14.12, 15.49, 14.36, 16.22])

z_re_array        =     np.array([0.002377, 0.003975, 0.004088, 0.004409, 0.005602, 0.005725,
                                    0.00639 , 0.01752 , 0.018064, 0.018146, 0.042275, 0.043318,
                                    0.046107, 0.054491, 0.059285, 0.060158, 0.062241, 0.063275,
                                    0.063515, 0.063775, 0.077493, 0.077701, 0.078068, 0.080837,
                                    0.083628, 0.083645, 0.094864, 0.097767, 0.098787, 0.113918,
                                    0.121342, 0.122881, 0.123596, 0.12389 , 0.124783, 0.135467,
                                    0.138527, 0.140754, 0.146789, 0.147135, 0.161068, 0.166588,
                                    0.167125, 0.170062, 0.187572, 0.187731, 0.192434, 0.194423,
                                    0.199618, 0.202627, 0.207116, 0.219098, 0.219449, 0.224832,
                                    0.225963, 0.238786, 0.249843, 0.250136, 0.282195, 0.292317,
                                    0.310529, 0.322338, 0.323049, 0.323387, 0.324526, 0.325663,
                                    0.328164, 0.347961, 0.349368, 0.360841, 0.363402, 0.365126,
                                    0.3786  , 0.386094, 0.39346 , 0.40133 , 0.417573, 0.420792,
                                    0.421322, 0.42188 , 0.423919, 0.424307, 0.446093, 0.44678 ,
                                    0.463201, 0.468468, 0.494986, 0.527928, 0.533232, 0.619111])

"""
for testing

#logNHI_array        =   np.array([12.88])
#z_re_array          =   np.array([0.364536])
logNHI_array        =   np.array([12.88, 12.93, 13.11, 16.61, 17.21, 17.54, 18.  , 18.  ,
       					19.23, 15.84, 15.89])
z_re_array          =   np.array([0.364536, 0.193186, 0.423067, 0.225963, 0.249843,
      				 0.463201, 0.080837, 0.121342, 0.00639 , 0.062241, 0.018064])
"""
#Array ranges hard coded
logT_array          =   np.around(np.arange(3.5,6.31,0.02), decimals = 2)
#logZ_array          =   np.around(np.arange(-3.0,-2.0,3.0), decimals = 2)
logZ_array	    =   np.array([-3.0,-1.0, -0.5,0.0]) #for AM
#logZ_array	    =   np.array([0.25,0.50,0.75,1.0]) for VK


uvb          = ['KS18']
uvb_Q        = [ 18]

logT         = []
z_red        = []
logZ         = []
logNHI       = []
uvb_models   = []
the_Q_values = []

for i in range(len(z_re_array)):
    for background in uvb:
        if background == 'KS18':
            for q in uvb_Q:
                for metal in logZ_array:
                    for temp in logT_array:
                        uvb_models.append(background)
                        the_Q_values.append(q)
                        logZ.append(metal)
                        logT.append(temp)
                        logNHI.append(logNHI_array[i])
                        z_red.append(z_re_array[i])


        else:   
            q = 18
            for redshift in z_re_array:
                for metal in logZ_array:
                    for temp in logT_array:
                        uvb_models.append(background)
                        the_Q_values.append(q)
                        logZ.append(metal)
                        logT.append(temp)
                        logNHI.append(logNHI_array[i])
                        z_red.append(z_re_array[i])

#-----write uvb fg and hm in cloudy format first
#path = '/mnt/quasar2/vikram/cloudy_run/metal_NH19_new'

#kwagrs = {'uvb' : 'P19', 'z' : 0.2}
#uvb_files(path, **kwagrs)

#kwagrs = {'uvb' : 'FG20', 'z' : 0.2}
#uvb_files(path, **kwagrs)
"""
f = open("try.txt","a")
s=[]
for  NHI, Z, Q, mod,T,z_re in zip(logNHI, logZ, the_Q_values, uvb_models,logT,z_red):
    t=(NHI,Z,Q,mod,T,z_re)
    s.append(t)
with open('try.txt','w') as f:
    for a in s:
        f.writelines('%s %s %s %s %s %s\n' % a)


"""

pool = mp.Pool(processes=46)
results = [pool.apply_async(run_parallel, args=(NHI, Z, Q, mod,T,z_re)) for  NHI, Z, Q, mod,T,z_re in zip(logNHI, logZ, the_Q_values, uvb_models,logT,z_red)]
output = [p.get() for p in results]
