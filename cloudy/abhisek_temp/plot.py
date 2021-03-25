import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt
import numpy as  np
import matplotlib as mpl
import astropy.table as tab
from astropy.table import Table
import matplotlib.pyplot as plt
import os


def make_plot(ifile1,ifile2,fig_name,cmap):

# setting the figure
    font = {'family': 'serif', 'weight': 'normal', 'size': 10}
    plt.rc('font', **font)
    mpl.rcParams['axes.linewidth'] = 1.5
    fig = plt.figure()
    fig.set_size_inches(30,5)
    ax = fig.add_subplot(121)

    tab2=np.genfromtxt(ifile1)
    tab1=np.genfromtxt(ifile2)
    hden=np.arange(-5,1,1)[np.newaxis]
    hden_t=hden.T
    data=np.abs(tab2-tab1)
    dataf=np.log10(data+1)
    ions = ["H", "H+",
                "He", "He+", "He+2",
                "C", "C+", "C+2", "C+3", "C+4", "C+5",
                "N", "N+", "N+2", "N+3", "N+4",
                "O", "O+", "O+2", "O+3", "O+4", "O+5", "O+6", "O+7",
                "S", "S+", "S+2", "S+3", "S+4", "S+5",
                "Si", "Si+", "Si+2", "Si+3", "Si+4",
                "Mg", "Mg+", "Mg+2",
                "Ne", "Ne+", "Ne+2", "Ne+3", "Ne+4", "Ne+5", "Ne+6", "Ne+7", "Ne+8",
                "Fe", "Fe+", "Fe+2",
                "Na", "Na+", "Na+2"]

    dataout=Table(dataf,names=ions)
    output='final.csv'
    dataout.write(output,format='csv',overwrite=True)
    read=pd.read_csv('final.csv')
    read.index=[r'$n_H$=-5',r'$n_H$=-4',r'$n_H$=-3',r'$n_H$=-2',r'$n_H$=-1',r'$n_H$=0']
    ax = sb.heatmap(read, xticklabels=1,cmap=cmap, cbar_kws={'label': 'Log (Difference b/w Uni and Bi-directional incident)'})
    ax.yaxis.set_tick_params(rotation=0)
    #plt.ylabel(r"Difference b/w UNI and BI-Directional")
    plt.xlabel(r"Ions")
    ax.tick_params(direction='out', which='major', length=5, width=1.50)
    ax.tick_params(direction='out', which='major', length=5, width=1.50)
    fig.savefig(fig_name, bbox_inches='tight',dpi=600)
    #os.remove("final.csv")
    return


"""
example:
from plot import make_plot as mk
mk('opt.spC','noopt.spC','plot_name.pdf','cmap')

t1=Table.read('opt_Q18.fits')
t2=Table.read('noopt_Q18.fits')
for s in ions:
	for i in range (0, 5):
		plt.scatter((t1[s]-t2[s])[i], t1['hden'][i],marker='s',facecolors='none')

Cmap options are:
rainbow,
MRmap_r', 'Dark2', 'Dark2_r', 'GnBu', 'GnBu_r', 'Greens', 'Greens_r', 'Greys', 'Greys_r', 'OrRd', 'OrRd_r', 'Oranges', 'Oranges_r', 'PRGn', 'PRGn_r', 'Paired', 'Paired_r', 'Pastel1', 'Pastel1_r', 'Pastel2', 'Pastel2_r', 'PiYG', 'PiYG_r', 'PuBu', 'PuBuGn', 'PuBuGn_r', 'PuBu_r', 'PuOr', 'PuOr_r', 'PuRd', 'PuRd_r', 'Purples', 'Purples_r', 'RdBu', 'RdBu_r', 'RdGy', 'RdGy_r', 'RdPu', 'RdPu_r', 'RdYlBu', 'RdYlBu_r', 'RdYlGn', 'RdYlGn_r', 'Reds', 'Reds_r', 'Set1', 'Set1_r', 'Set2', 'Set2_r', 'Set3', 'Set3_r', 'Spectral', 'Spectral_r', 'Wistia', 'Wistia_r', 'YlGn', 'YlGnBu', 'YlGnBu_r', 'YlGn_r', 'YlOrBr', 'YlOrBr_r', 'YlOrRd', 'YlOrRd_r', 'afmhot', 'afmhot_r', 'autumn', 'autumn_r', 'binary', 'binary_r', 'bone', 'bone_r', 'brg', 'brg_r', 'bwr', 'bwr_r', 'cividis', 'cividis_r', 'cool', 'cool_r', 'coolwarm', 'coolwarm_r', 'copper', 'copper_r', 'crest', 'crest_r', 'cubehelix', 'cubehelix_r', 'flag', 'flag_r', 'flare', 'flare_r', 'gist_earth', 'gist_earth_r', 'gist_gray', 'gist_gray_r', 'gist_heat', 'gist_heat_r', 'gist_ncar', 'gist_ncar_r', 'gist_rainbow', 'gist_rainbow_r', 'gist_stern', 'gist_stern_r', 'gist_yarg', 'gist_yarg_r', 'gnuplot', 'gnuplot2', 'gnuplot2_r', 'gnuplot_r', 'gray', 'gray_r', 'hot', 'hot_r', 'hsv', 'hsv_r', 'icefire', 'icefire_r', 'inferno', 'inferno_r', 'jet', 'jet_r', 'magma', 'magma_r', 'mako', 'mako_r', 'nipy_spectral', 'nipy_spectral_r', 'ocean', 'ocean_r', 'pink', 'pink_r', 'plasma', 'plasma_r', 'prism', 'prism_r', 'rainbow', 'rainbow_r', 'rocket', 'rocket_r', 'seismic', 'seismic_r', 'spring', 'spring_r', 'summer', 'summer_r', 'tab10', 'tab10_r', 'tab20', 'tab20_r', 'tab20b', 'tab20b_r', 'tab20c', 'tab20c_r', 'terrain', 'terrain_r', 'turbo', 'turbo_r', 'twilight', 'twilight_r', 'twilight_shifted', 'twilight_shifted_r', 'viridis', 'viridis_r', 'vlag', 'vlag_r', 'winter', 'winter_r'
"""
