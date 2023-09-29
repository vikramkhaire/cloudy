import matplotlib.pyplot as plt
import matplotlib
import astropy.table as t
from astropy.io import ascii
import glob
import numpy as np
import matplotlib.colors as col
import seaborn as sns
import os
from astropy import coordinates
import astropy.units as u

# run this code in the directory downloaded from Danforth sample's MAST webpage
def find_nearest_id(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

#------------------------------------- Reading Danfoth's whole sample
# grab files from
files=glob.glob('*/*linelist.txt', recursive=True)

# I want to make a big table from all files of 82 QSOs
# read first file to store the meta file
data = t.Table.read(files[0], format='ascii', fast_reader=False) # reading a file 
alltable=data[:0].copy() # for keeping header OR just keeping the same table format
count=0
for file in files:
    data = t.Table.read(file, format='ascii', fast_reader=False)
    data['z'] = float(data.meta['comments'][16].split("/")[0].split("=")[-1])
    data.meta=None
    data['qname']=os.path.dirname(file) # for keeping qso name in the alltable
    count=count+len(data)
    alltable=t.vstack([alltable, data], metadata_conflicts='silent')

print(count, 'length of the whole dataset')

# store the table for future use
alltable.write('alltable_alllines.dat', format='ascii', overwrite=True)


#------------------------------- Finding Si contamination
# read Sapna'a file
stack = t.Table.read('485_OVI_stack_lya_abs.txt', format = 'ascii')
#SiII => 1193.2897, SiIII => 1206.500
wave_siII = 1193.2897
wave_siIII = 1206.5
c = 3.0e5 #km/s
# counters
count_siII =0
count_siIII = 0
count = 0
# arrays to store contamination info
SiII_cont_line = []
SiII_cont_delv = []
SiII_cont_num = []

SiIII_cont_line = []
SiIII_cont_delv = []
SiIII_cont_num = []

# line by line search
for i in range(len(stack)):
    #print(i)
    qname =  stack['Dirname'][i]
    zlya = stack['z_abs_lya'][i]
    delta_v = 3*2.355*stack['b_lya'][i] # 3 times fwhm of the lya
    delta_wave = delta_v* 1200*(1+zlya)/c # 1200 is approximation for rest wavelength of both Si lines
    #print(delta_v, stack['b_lya'][i])
    sort_table = alltable[alltable['qname'] == qname]
    if len(sort_table) == 0:
        print("Warning quasar name is not recognised")
    # wavelength region of Si II    
    w_obs = (1+ zlya)*wave_siII
    w_start =  w_obs - 0.5*delta_wave
    w_end = w_obs + 0.5 *delta_wave
    check_siII = sort_table[sort_table['col1'] < w_end]
    check_siII = check_siII[check_siII['col1'] > w_start]
    
    flag  = 0
    if len(check_siII)>0:
        count_siII = count_siII +1
        #print('line is contaminated')
        # number of lines contaminating
        SiII_cont_num.append(len(check_siII))
        if len(check_siII) > 1:
            j = find_nearest_id(np.array(check_siII['col1']), w_obs)
        else:
            j=0
        # line id
        SiII_cont_line.append(check_siII['col2'][j])
        # wavelength
        del_v_cont = (w_obs - check_siII['col1'][j])*c/w_obs
        SiII_cont_delv.append(del_v_cont)
        flag = 1
    else:
        SiII_cont_num.append(0)
        SiII_cont_delv.append(9999)
        SiII_cont_line.append('None')
        
    
    # wavelength region of Si II    
    w_obs = (1+ zlya)*wave_siIII
    w_start =  w_obs - 0.5*delta_wave
    w_end = w_obs + 0.5 *delta_wave
    check_siIII = sort_table[sort_table['col1'] < w_end]
    check_siIII = check_siIII[check_siIII['col1'] > w_start]
    
    if len(check_siIII)>0:
        count_siIII = count_siIII +1
        #print('line is contaminated', check_siIII)
        # number of lines contaminating
        SiIII_cont_num.append(len(check_siIII))
        if len(check_siIII) > 1:
            j = find_nearest_id(np.array(check_siIII['col1']), w_obs)
        else:
            j=0
        # line id
        SiIII_cont_line.append(check_siIII['col2'][j])
        # wavelength
        del_v_cont = (w_obs - check_siIII['col1'][j])*c/w_obs
        SiIII_cont_delv.append(del_v_cont)
        
        
        if flag ==1:
            count = count +1
            
    else:
        SiIII_cont_num.append(0)
        SiIII_cont_delv.append(9999)
        SiIII_cont_line.append('None')

# adding the columns to the table
stack.add_columns([SiII_cont_num, SiII_cont_line, SiII_cont_delv], names = ('SiII_cont_num', 'SiII_cont_line', 'SiII_cont_delv'))
stack.add_columns([SiIII_cont_num, SiIII_cont_line, SiIII_cont_delv], names = ('SiIII_cont_num', 'SiIII_cont_line', 'SiIII_cont_delv'))

# writing final table
stack.write('with_Si_contamination.txt', format = 'ascii')

print('# siII contmination ', count_siII)
print('# siIII contmination ', count_siIII)
print('# both contmination ', count)
