import netCDF4
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import glob
import scipy.io as sio

original_files = [line.rstrip('\n') for line in open('ahitest.txt')]
clear_files = sorted(glob.glob("data/clear/*.nc"))
pass2_files = sorted(glob.glob("data/pass2/*.nc"))
TL0_files = sorted(glob.glob("data/TL0/*.nc"))
TL1_files = sorted(glob.glob("data/TL1/*.nc"))
TL2_files = sorted(glob.glob("data/TL2/*.nc"))


orig_filenames = [x.split('/')[6] for x in original_files]
clear_filenames = [x.split('/')[-1] for x in clear_files]
pass2_filenames = [x.split('/')[-1] for x in pass2_files]

def read_var(cdf,variable):
    data = np.squeeze(cdf[variable][:])
    return data


points = [(1631,2621), (1695,2316), (785,2316), (3676,1006), (4250,3767), (1800, 1690)]
#y_lims = [(299,302),   (290,295),   (292.5,295.5), (296,299),   (294,297),   (290,295),   (298,303),   (284.5,287.5), (298,302) ]
times = []

hour_ind = []
x_lims = []

for i,pass2_file in enumerate(pass2_filenames):
    if pass2_file[10:12] == '00':
        times.append(pass2_file[8:12])
        hour_ind.append(i)
    if pass2_file[8:12] == '0000':
        x_lims.append(i)

x_lims = tuple(x_lims)

sst = {}
clear = {}
pass2 = {}
TL0 = {}
TL1 = {}
TL2 = {}
h = {}

for i in points:    
    sst[i] = []
    clear[i] = []
    pass2[i] = []
    TL0[i] = []
    TL1[i] = []
    TL2[i] = []
    h[i] = []

total_files = len(pass2_files)
for i in range(total_files):

    base_file = pass2_files[i].split('/')[2]
        
    pass2nc = netCDF4.Dataset(pass2_files[i])
    pass2_vals = read_var(pass2nc,'brightness_temperature_11um2')
    pass2nc.close()

    TL0nc = netCDF4.Dataset(TL0_files[i])
    TL0_vals = read_var(TL0nc,'brightness_temperature_11um2')
    TL0nc.close()

    TL1nc = netCDF4.Dataset(TL1_files[i])
    TL1_vals = read_var(TL1nc,'brightness_temperature_11um2')
    TL1nc.close()

    TL2nc = netCDF4.Dataset(TL2_files[i])
    TL2_vals = read_var(TL2nc,'brightness_temperature_11um2')
    TL2nc.close()

    filename = clear_files[clear_filenames.index(base_file)]

    clearnc = netCDF4.Dataset(filename)
    clear_vals = read_var(clearnc,'brightness_temperature_11um2')
    clearnc.close()

    filename = original_files[orig_filenames.index(base_file)]

    orignc=netCDF4.Dataset(filename)
    sst_val = read_var(orignc,'sea_surface_temperature')
    orignc.close()

    print filename

    for j in points:
        
        sst[j].append(sst_val[j[0],j[1]])
        
        TL0[j].append(TL0_vals[j[0],j[1]])

        TL1[j].append(TL1_vals[j[0],j[1]])

        TL2[j].append(TL2_vals[j[0],j[1]])

        if(clear_vals[j[0],j[1]] == 255):
            clear[j].append(sst[j][i])
        else:
        	clear[j].append(np.nan)

        if(pass2_vals[j[0],j[1]] == 255):
            pass2[j].append(sst[j][i])
        else:
            pass2[j].append(np.nan)
    

for j in points:
    for i in hour_ind:
        h[j].append(pass2[j][i])

x_inds = range(total_files)
for j,i in enumerate(points):
    
    fig = plt.figure(figsize=(12,10))
    fig.canvas.set_window_title(str(i))
    plt.title("2017 - 01 - 09   Location: " + str(i), fontsize=20)
    plt.grid()
    plt.plot(x_inds,sst[i], 'b', label="Original SST")
    plt.plot(x_inds,sst[i], 'b.')
    plt.plot(x_inds,clear[i], 'bo', fillstyle='none', markeredgewidth=2, markersize=12, label="Clear Pass 1")
    plt.plot(x_inds, pass2[i], '.r', markersize=15, label="Clear Pass 2")
    plt.plot(x_inds, TL0[i], c='g', label="TL0")
    plt.plot(x_inds, TL1[i], c='m', label="TL1")
    plt.plot(x_inds, TL2[i], c='c', label="TL2")
    plt.xticks(hour_ind, times, rotation=-45, fontsize=15)
    plt.xlim(x_lims)
    plt.legend(fontsize=15)
    plt.yticks(fontsize=15)
    
    plt.show()
