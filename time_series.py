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
interp_files = sorted(glob.glob("data/collated_interp/*.nc"))
pass2_files = sorted(glob.glob("data/pass2/*.nc"))
#smooth_collate_files = sorted(glob.glob("data/smooth_collate/*.nc"))

interp_files == [x for x in interp_files if x[31:33] == '00']
original_files = [line.rstrip('\n') for line in open('ahitest.txt')]

orig_filenames = [x.split('/')[6] for x in original_files]
pass2_filenames = [x.split('/')[2] for x in pass2_files]
approx_filenames = [x.split('/')[2] for x in approx_files]
clear_filenames = [x.split('/')[2] for x in clear_files]
smooth_filenames = [x.split('/')[2] for x in smooth_files]
#smooth_collate_filenames = [x.split('/')[2] for x in smooth_collate_files]

filenames = [x.split('/')[2] for x in smooth_files]
original_filenames = [x.split('/')[6] for x in original_files]

TL0_files = []
TL1_files = []
TL2_files = []
for i in range(len(pass2_files)):
    TL0_files.append("data/TL0/TL0_"+str(i)+".nc")
    TL1_files.append("data/TL1/TL1_"+str(i)+".nc")
    TL2_files.append("data/TL2/TL2_"+str(i)+".nc")

filter_window_lag = 1
second_pass_lag = 25/2
smooth_lag = 19/2
collated_inds = []

def read_var(cdf,variable):
    data = np.squeeze(cdf[variable][:])
    return data


points = [(3351, 3398),(3353, 3400),(2639, 4975),(544,3237),(524,3173),(573,3190),(2664,5209),(2692,5198),(2701,5174),(2705,5186),(3386,1678),(3386,1677)]
#points = [(1878,3085),(1879,3081)]
points = [(1471, 875), (1477, 1647), (1539, 624), (1555, 1538), (1590, 640), (1596, 640), (1603, 638), (1706, 694), (1712, 693), (1723, 691), (1725, 524), (1726, 521)]

times = []

for approx_file in approx_files:
    times.append(approx_file[20:24])

times = times[0:len(times):10]

approx = {}
sst = {}
second_pass = {}
clear = {}
smooth = {}
collated = {}
smooth_collate = {}

for i in points:
    approx[i] = []
    sst[i] = []
    clear[i] = []
    smooth[i] = []
    collated[i] = []
    second_pass[i] = []
    smooth_collate[i] = []

total_files = len(approx_files)
for i,f in enumerate(approx_files):

    approxnc=netCDF4.Dataset(f)
    approx_val = read_var(approxnc,'sea_surface_temperature')
    approxnc.close()

    filename = smooth_files[smooth_filenames.index(f.split('/')[2])]

    smoothnc=netCDF4.Dataset(filename)
    smooth_val = read_var(smoothnc,'sea_surface_temperature')
    smoothnc.close()

    
    #filename = smooth_collate_files[smooth_collate_filenames.index(f.split('/')[2])]
    #smoothnc=netCDF4.Dataset(filename)
    #smooth_collate_val = read_var(smoothnc,'sea_surface_temperature')
    #smoothnc.close()
    
    filename = pass2_files[pass2_filenames.index(f.split('/')[2])]

    pass2nc=netCDF4.Dataset(filename)
    sst_clear_pass2 = read_var(pass2nc,'brightness_temperature_11um2')
    pass2nc.close()

    filename = clear_files[clear_filenames.index(f.split('/')[2])]

    clearnc=netCDF4.Dataset(filename)
    sst_clear_val = read_var(clearnc,'brightness_temperature_11um2')
    clearnc.close()

    filename = original_files[orig_filenames.index(f.split('/')[2])]

    orignc=netCDF4.Dataset(filename)
    sst_val = read_var(orignc,'sea_surface_temperature')
    orignc.close()


    for j in points:
        
        sst[j].append(sst_val[j[0],j[1]])

        approx[j].append(approx_val[j[0],j[1]])
          
        smooth[j].append(smooth_val[j[0],j[1]])
        #smooth_collate[j].append(smooth_collate_val[j[0],j[1]])
        if(sst_clear_val[j[0],j[1]] == 255):
        	clear[j].append(sst[j][i])
        else:
        	clear[j].append(np.nan)

        if(sst_clear_pass2[j[0],j[1]]):
            second_pass[j].append(sst[j][i])
        else:
        	second_pass[j].append(np.nan)

    

collated_inds = []
for f in collated_files:
    
    gen=netCDF4.Dataset(f)
    original = read_var(gen,'sea_surface_temperature')
    gen.close()
    for j in points:
        collated[j].append(original[j[0],j[1]])

    #filename = interp_files[i]
    
    #gen=netCDF4.Dataset(filename)
    #original = read_var(gen,'sea_surface_temperature')
    #for j in points:
    #    interpolated[j].append(original[j[0],j[1]])
    collated_inds.append(approx_filenames.index(f.split('/')[2])) 
    
"""
smooth_collate_inds = []
for f in smooth_collate_files:
    
    gen=netCDF4.Dataset(f)
    original = read_var(gen,'sea_surface_temperature')
    gen.close()
    for j in points:
        smooth_collate[j].append(original[j[0],j[1]])
    smooth_collate_inds.append(approx_filenames.index(f.split('/')[2])) 
"""

    

for i in points:
    outfile = str(i)
    
    #sio.savemat(outfile,{"sst":sst[i],'approx':approx[i],"second_pass":second_pass[i],"clear":clear[i],"smooth":smooth[i],"interpolated":interpolated[i],"collated":collated[i], "times":times})
    
    fig = plt.figure(figsize=(12,10))
    fig.canvas.set_window_title(str(i))
    plt.title(str(i) + "pchip smooth\n")
    plt.plot(range(total_files),approx[i],'.',label="approx")
    #plt.plot(smooth_collate_inds,smooth_collate[i],label="smooth collate")
    plt.plot(range(total_files),second_pass[i],'o', markersize=15,label="pass2")
    plt.plot(range(total_files),clear[i],'.', markersize=15,label="clear")
    plt.plot(range(total_files),sst[i],label="orig")
    plt.plot(range(total_files),smooth[i],label="smooth")
    
    #plt.plot(collated_inds,interpolated[i],'.',markersize=15,label="interp")
    plt.plot(collated_inds,collated[i],'.',markersize=15,label="collated")
    plt.legend()
    #plt.xticks(range(times),times,rotation='vertical')
    plt.show()
