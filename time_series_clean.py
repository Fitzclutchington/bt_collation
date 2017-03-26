import netCDF4
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import glob
import scipy.io as sio

original_files = [line.rstrip('\n') for line in open('ahitest.txt')]
clear_files = sorted(glob.glob("data/clear/*.nc"))
smooth_files = sorted(glob.glob("data/smooth_test/*.nc"))
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


#points = [(3351, 3398),(3353, 3400),(2639, 4975),(544,3237),(524,3173),(573,3190),(2664,5209),(2692,5198),(2701,5174),(2705,5186),(3386,1678),(3386,1677)]
#points = [(1878,3085),(1879,3081)]
#points = [(1471, 875), (1477, 1647), (1539, 624), (1555, 1538), (1590, 640), (1596, 640), (1603, 638), (1706, 694), (1712, 693), (1723, 691), (1725, 524), (1726, 521)]
points = [(4451,3534),(3773,1502),(1577,2328),(1894,4455)]
times = []

hour_ind = []

for i,approx_file in enumerate(approx_filenames):
    if approx_file[10:12] == '00':
        times.append(approx_file[8:12])
        hour_ind.append(i)

 
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

dicts = [ sst, smooth, approx, approx2, collated, collated2, smooth_collate, reinstated, clear, second_pass, h ]

for d in dicts:
    for i in points:    
        d[i] = []

dicts = [ approx2, collated2, collated, reinstated, smooth_collate ]
data = [ approx2_files, collated_files2, collated_files, reinstated_files, smooth_collate_files ]
masks = [ pass2_files  ]
masks_data = [ pass2 ]

total_files = len(approx_files)
for i in range(total_files):

    base_file = approx_files[i].split('/')[2]
    filename = original_files[orig_filenames.index(base_file)]

    orignc=netCDF4.Dataset(filename)
    sst_val = read_var(orignc,'sea_surface_temperature')
    orignc.close()

    for j in points:            
        sst[j].append(sst_val[j[0],j[1]])


    for f,d in zip(data,dicts):
        cdf=netCDF4.Dataset(f)
        val = read_var(cdf,'sea_surface_temperature')
        cdf.close()

        for j in points:            
            d[j].append(val[j[0],j[1]])
        
    for f,d in zip(mask,masks_data):
        filename = pass2_files[pass2_filenames.index(base_file)]
        cdf = netCDF4.Dataset(filename)
        sst_mask = read_var(cdf,'brightness_temperature_11um2')
        cdf.close()

        for j in points:
            if(sst_mask[j[0],j[1]]):
                d[j].append(sst[j][i])
            else:
                d[j].append(np.nan)
    

for j in points:
    for i in hour_ind:
        h[j].append(collated2[j][i])

for i in points:
    outfile = str(i)
    
    #sio.savemat(outfile,{"sst":sst[i],'approx':approx[i],"second_pass":second_pass[i],"clear":clear[i],"smooth":smooth[i],"interpolated":interpolated[i],"collated":collated[i], "times":times})
    
    fig = plt.figure(figsize=(12,10))
    fig.canvas.set_window_title(str(i))
    plt.title(str(i) + " pchip smooth\n")
    #plt.plot(range(total_files),approx[i],'.',label="approx")
    plt.plot(range(total_files),approx2[i],'.',label="approx2")
    plt.plot(range(total_files),smooth_collate[i],label="smooth collate")
    plt.plot(range(total_files),second_pass[i],'o', markersize=15,label="pass2")
    plt.plot(range(total_files),reinstated[i],'o', markersize=15,label="reinstated")
    #plt.plot(range(total_files),clear[i],'.', markersize=15,label="clear")
    plt.plot(range(total_files),sst[i],label="orig")
    #plt.plot(range(total_files),smooth[i],label="smooth")
    
    #plt.plot(collated_inds,interpolated[i],'.',markersize=15,label="interp")
    plt.plot(range(total_files),collated2[i],'.',markersize=15,label="collated2")
    plt.plot(hour_ind,h[i],'.',markersize=15,label="hour")
    plt.legend()
    plt.xticks(hour_ind,times,rotation='vertical')
    plt.show()
