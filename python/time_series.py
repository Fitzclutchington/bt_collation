import netCDF4
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import glob
import scipy.io as sio

original_files = [line.rstrip('\n') for line in open('ahitest.txt')]
#original_files = sorted(glob.glob("../data/2016-11-15/*.nc"))
clear_files = sorted(glob.glob("data/clear/*.nc"))
smooth_files = sorted(glob.glob("data/smooth/*.nc"))
approx_files = sorted(glob.glob("data/approx/*.nc"))
collated_files = sorted(glob.glob("data/collated_mat/*.nc"))
approx2_files = sorted(glob.glob("data/approx2/*.nc"))
collated2_files = sorted(glob.glob("data/collated_mat2/*.nc"))
reinstated_files = sorted(glob.glob("data/reinstated/*.nc"))
pass2_files = sorted(glob.glob("data/pass2/*.nc"))
smooth_collate_files = sorted(glob.glob("data/smooth_collate/*.nc"))


original_files = [line.rstrip('\n') for line in open('ahitest.txt')]

orig_filenames = [x.split('/')[6] for x in original_files]
pass2_filenames = [x.split('/')[2] for x in pass2_files]
approx_filenames = [x.split('/')[2] for x in approx_files]
clear_filenames = [x.split('/')[2] for x in clear_files]
smooth_filenames = [x.split('/')[2] for x in smooth_files]
smooth_collate_filenames = [x.split('/')[2] for x in smooth_collate_files]

filenames = [x.split('/')[2] for x in smooth_files]
original_filenames = [x.split('/')[6] for x in original_files]
"""
TL0_files = []
TL1_files = []
TL2_files = []
for i in range(len(pass2_files)):
    TL0_files.append("data/TL0/TL0_"+str(i)+".nc")
    TL1_files.append("data/TL1/TL1_"+str(i)+".nc")
    TL2_files.append("data/TL2/TL2_"+str(i)+".nc")
"""
filter_window_lag = 1
second_pass_lag = 25/2
smooth_lag = 19/2
collated_inds = []

def read_var(cdf,variable):
    data = np.squeeze(cdf[variable][:])
    return data


points = [(1738,2880), (1533,1649), (4255,1488),   (3975,1405), (4443,3479), (4516,3858), (1885,4360), (4763,2138),   (1906,507)]
points = [(3922,1271), (3920,1288), (3939,1282), (4270, 1547), (4270,1543), (4269,1542) ]
#y_lims = [(299,302),   (290,295),   (292.5,295.5), (296,299),   (294,297),   (290,295),   (298,303),   (284.5,287.5), (298,302) ]
times = []

hour_ind = []
x_lims = []

for i,approx_file in enumerate(approx_filenames):
    if approx_file[10:12] == '00':
        times.append(approx_file[8:12])
        hour_ind.append(i)
    if approx_file[8:12] == '0000':
        x_lims.append(i)

x_lims = tuple(x_lims)

sst = {}
second_pass = {}
clear = {}
smooth = {}
approx = {}
collated = {}
approx2 = {}
collated2 = {}
smooth_collate = {}
reinstated = {}
h = {}

for i in points:    
    sst[i] = []
    clear[i] = []
    smooth[i] = []
    collated[i] = []
    approx[i] = []
    collated2[i] = []
    approx2[i] = []
    second_pass[i] = []
    smooth_collate[i] = []
    reinstated[i] = []
    h[i] = []

total_files = len(approx_files)
for i in range(total_files):

    base_file = approx_files[i].split('/')[2]
    
    approxnc = netCDF4.Dataset(approx_files[i])
    approx_val = read_var(approxnc,'sea_surface_temperature')
    approxnc.close()

    smoothnc = netCDF4.Dataset(smooth_files[i])
    smooth_val = read_var(smoothnc,'sea_surface_temperature')
    smoothnc.close()
    
    approxnc = netCDF4.Dataset(approx2_files[i])
    approx2_val = read_var(approxnc,'sea_surface_temperature')
    approxnc.close()

    filename = smooth_collate_files[i]
    smoothnc=netCDF4.Dataset(filename)
    smooth_collate_val = read_var(smoothnc,'sea_surface_temperature')
    smoothnc.close()
    
    filename = pass2_files[pass2_filenames.index(base_file)]

    pass2nc=netCDF4.Dataset(filename)
    sst_clear_pass2 = read_var(pass2nc,'brightness_temperature_11um2')
    pass2nc.close()

    filename = clear_files[clear_filenames.index(base_file)]

    collatednc = netCDF4.Dataset(collated_files[i])
    collated_vals = read_var(collatednc,'sea_surface_temperature')
    collatednc.close()
        
    collatednc = netCDF4.Dataset(collated2_files[i])
    collated2_vals = read_var(collatednc,'sea_surface_temperature')
    collatednc.close()

    reinstatednc = netCDF4.Dataset(reinstated_files[i])
    reinstated_vals = read_var(reinstatednc,'sea_surface_temperature')
    reinstatednc.close()

    filename = original_files[orig_filenames.index(base_file)]

    orignc=netCDF4.Dataset(filename)
    sst_val = read_var(orignc,'sea_surface_temperature')
    orignc.close()

    print "finished " + filename

    for j in points:
        
        sst[j].append(sst_val[j[0],j[1]])

        approx[j].append(approx_val[j[0],j[1]])
          
        smooth[j].append(smooth_val[j[0],j[1]])

        approx2[j].append(approx2_val[j[0],j[1]])

        collated[j].append(collated_vals[j[0],j[1]])
        
        collated2[j].append(collated2_vals[j[0],j[1]])

        smooth_collate[j].append(smooth_collate_val[j[0],j[1]])

        reinstated[j].append(reinstated_vals[j[0],j[1]])

        if(sst_clear_pass2[j[0],j[1]]):
            second_pass[j].append(sst[j][i])
        else:
        	second_pass[j].append(np.nan)

    

for j in points:
    for i in hour_ind:
        h[j].append(collated2[j][i])

x_inds = range(total_files)
for j,i in enumerate(points):
    
    fig = plt.figure(figsize=(12,10))
    fig.canvas.set_window_title(str(i))
    plt.title("2017 - 01 - 10   Location: " + str(i), fontsize=20)
    plt.grid()
    plt.plot(x_inds, approx[i], '-', c="#D98719", label="Approx 1")
    plt.plot(x_inds, approx2[i], '-m', label="Approx 2")
    plt.plot(x_inds, second_pass[i], 'o', c="#39FF14",fillstyle='none', markeredgewidth=2, markersize=12, label="Clear 1")
    plt.plot(x_inds, reinstated[i], '.r', markersize=15, label="Clear 2")
    plt.plot(x_inds,sst[i], 'b', label="Original SST")
    plt.plot(x_inds,smooth[i], 'g', label="Smooth")
    plt.plot(x_inds,collated[i], 'r', label="Collated 1")
    plt.plot(x_inds,collated2[i], '-k', lw=3, label="Collated 2")
    plt.plot(hour_ind,h[i], '.c', markersize=20, label="Hour")
    plt.legend(fontsize=15)
    plt.yticks(fontsize=15)
    plt.xticks(hour_ind, times, rotation=-45, fontsize=15)
    plt.xlim(x_lims)
    #plt.ylim(y_lims[j])
    plt.show()
