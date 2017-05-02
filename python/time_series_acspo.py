import netCDF4
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import glob
import scipy.io as sio
import sys

if len(sys.argv) < 2:
    print "usage: python time_series.py <point_file>"
    sys.exit()

point_file = sys.argv[1]

original_files = [line.rstrip('\n') for line in open('../ahitest.txt')]
collated2_files = sorted(glob.glob("../data/collated_mat2/*.nc"))
reinstated_files = sorted(glob.glob("../data/reinstated/*.nc"))
pass1_files = sorted(glob.glob("../data/clear/*.nc"))

orig_filenames = [x.split('/')[-1] for x in original_files]
collated2_filenames = [x.split('/')[-1] for x in collated2_files]
pass1_filenames = [x.split('/')[-1] for x in pass1_files]

def read_var(cdf,variable):
    data = np.squeeze(cdf[variable][:])
    return data

def read_points(point_file):
    points = []
    with open(point_file, 'r') as f:
        for line in f:
            point = map(int,line.split(','))
            points.append((point[0],point[1]))
    return points

#points = [(1738,2880), (1533,1649), (4255,1488),   (3975,1405), (4443,3479), (4516,3858), (1885,4360), (4763,2138),   (1906,507)]
points = read_points(point_file)
#y_lims = [(299,302),   (290,295),   (292.5,295.5), (296,299),   (294,297),   (290,295),   (298,303),   (284.5,287.5), (298,302) ]
times = []

hour_ind = []
x_lims = []

for i,collated2_file in enumerate(collated2_filenames):
    if collated2_file[10:12] == '00':
        times.append(collated2_file[8:12])
        hour_ind.append(i)
    if collated2_file[8:12] == '0000':
        x_lims.append(i)

x_lims = tuple(x_lims)

sst = {}
acspo = {}
collated2 = {}
reinstated = {}
pass1 = {}
h = {}

for i in points:    
    sst[i] = []
    collated2[i] = []
    acspo[i] = []
    reinstated[i] = []
    pass1[i] = []
    h[i] = []

total_files = len(collated2_files)
for i in range(total_files):

    base_file = collated2_files[i].split('/')[-1]
        
    collatednc = netCDF4.Dataset(collated2_files[i])
    collated2_vals = read_var(collatednc,'sea_surface_temperature')
    collatednc.close()

    reinstatednc = netCDF4.Dataset(reinstated_files[i])
    reinstated_vals = read_var(reinstatednc,'sea_surface_temperature')
    reinstatednc.close()

    filename = original_files[orig_filenames.index(base_file)]

    orignc=netCDF4.Dataset(filename)
    sst_val = read_var(orignc,'sea_surface_temperature')
    acspo_val = read_var(orignc, 'l2p_flags')
    orignc.close()

    acspo_val = np.bitwise_and(acspo_val,-16384).astype(bool)

    filename = pass1_files[pass1_filenames.index(base_file)]

    orignc=netCDF4.Dataset(filename)
    clear_val = read_var(orignc,'brightness_temperature_11um2')
    orignc.close()


    for j in points:
        
        sst[j].append(sst_val[j[0],j[1]])
        
        collated2[j].append(collated2_vals[j[0],j[1]])

        reinstated[j].append(reinstated_vals[j[0],j[1]])

        if(not acspo_val[j[0],j[1]]):
            acspo[j].append(sst[j][i])
        else:
        	acspo[j].append(np.nan)

        if clear_val[j[0],j[1]] == 255:
            pass1[j].append(sst[j][i])
        else:
            pass1[j].append(np.nan)
    

for j in points:
    for i in hour_ind:
        h[j].append(collated2[j][i])

x_inds = range(total_files)
for j,i in enumerate(points):
    
    fig = plt.figure(figsize=(12,10))
    fig.canvas.set_window_title(str(i))
    plt.title("2017 - 01 -08   Location: " + str(i), fontsize=20)
    plt.grid()
    plt.plot(x_inds,sst[i], 'b', label="Original SST")
    plt.plot(x_inds,sst[i], 'b.')
    plt.plot(x_inds,acspo[i], 'bo', fillstyle='none', markeredgewidth=2, markersize=12, label="ACSPO Mask")
    plt.plot(x_inds, reinstated[i], '.r', markersize=15, label="Reinstated")
    plt.plot(x_inds,collated2[i], '-k', linewidth=3, markersize=15, label="Collated 2")
    plt.plot(x_inds,pass1[i], '*',  markersize=20, label="pass1")
    plt.plot(hour_ind,h[i], '.', c='c', markersize=20, label="Hour")
    plt.xticks(hour_ind, times, rotation=-45, fontsize=15)
    #plt.xlim(x_lims)
    #plt.ylim(y_lims[j])
    plt.legend(fontsize=15)
    plt.yticks(fontsize=15)
    
    plt.show()
