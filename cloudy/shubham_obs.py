from cloudy import mcmc_nH_Z_observations as mcmc
import numpy as np
import astropy.table as tab

ions_to_use= ['Si+', 'N+2', 'C+']
data_col = np.array([14.47,14.5,14.79])
sigma_col = np.array([0.14,0.02,0.02])
true_Q =18

outpath = '/home/vikram/cloudy_run/shubham/inference'
model_path  = '/home/vikram/cloudy_run/shubham'
outfile = outpath + '/metal_NH18_2D.fits'

uvb_array = ['KS18']
Q_array= [18]

out_tab =  tab.Table()
for uvb, q in zip(uvb_array, Q_array):
    name =uvb + '_Q{}'.format(q)
    figname = outpath + '/' + name + '.pdf'

    flat_samples, ndim = mcmc.run_mcmc(model_path= model_path, Q_uvb=q, ions_to_use=ions_to_use,
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

