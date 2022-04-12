#------------------
import os
import numpy as np
import glob
import astropy.table as tab

output_path = '/home/vikram/data/cloudy/output'
list_of_files = glob.glob(output_path + '/' + '*.npy')

redshift_array = []
logZ_array = []
logZ_err1= []
logZ_err2 = []
logT_array = []
logT_err1= []
logT_err2 = []
lognH_array = []
lognH_err1= []
lognH_err2 = []

for chain_file in list_of_files:
    redshift = int((chain_file.split('/z_')[-1]).split('_')[0])/1e6
    redshift_array.append(redshift)
    print('for z {}'.format(redshift))

    # read chain
    flat_samples = np.load(chain_file)

    # thin chain by factor 2
    flat_samples =  flat_samples[::2, :]

    i=0
    mcmc = np.percentile(flat_samples[:, i], [16, 50, 84])
    q = np.diff(mcmc)
    lognH_array.append(mcmc[1])
    lognH_err1.append(q[0])
    lognH_err2.append(q[1])

    i=1
    mcmc = np.percentile(flat_samples[:, i], [16, 50, 84])
    q = np.diff(mcmc)
    logZ_array.append(mcmc[1])
    logZ_err1.append(q[0])
    logZ_err2.append(q[1])

    i=2
    mcmc = np.percentile(flat_samples[:, i], [16, 50, 84])
    q = np.diff(mcmc)
    logT_array.append(mcmc[1])
    logT_err1.append(q[0])
    logT_err2.append(q[1])



# write table
table_file = output_path + '/results.fits'
res = tab.Table(
    [redshift_array, lognH_array, lognH_err1, lognH_err2, logZ_array, logZ_err1, logZ_err2, logT_array, logT_err1, logT_err2],
    names=('z', 'lognH', 'lognHerr1', 'lognHerr2', 'logZ', 'logZerr1', 'logZerr2', 'logT', 'logTerr1', 'logTerr2' ))
res.write(table_file, overwrite = True)

print(res)

