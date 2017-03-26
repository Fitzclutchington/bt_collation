import netCDF4
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import glob
import scipy.io as sio

original_files = [line.rstrip('\n') for line in open('ahitest.txt')]
#original_files = sorted(glob.glob("../data/2016-11-15/*.nc"))
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
points = [(3790,1683),(1564,1535),(1564,1536),(4428,3346),(4438,3347),(729,2917),(735,2917),(726,2870),(730,2870),(1487,1714),(1487,1720),(1494,1715),(437,4080),(442,4082),(802,4209),(3210,5058),(3218,5064)]
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
    """
    approxnc=netCDF4.Dataset(f)
    approx_val = read_var(approxnc,'sea_surface_temperature')
    approxnc.close()

    filename = smooth_files[smooth_filenames.index(f.split('/')[2])]

    smoothnc=netCDF4.Dataset(filename)
    smooth_val = read_var(smoothnc,'sea_surface_temperature')
    smoothnc.close()
    """
    approxnc=netCDF4.Dataset(approx2_files[i])
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

    """
    clearnc=netCDF4.Dataset(filename)
    sst_clear_val = read_var(clearnc,'brightness_temperature_11um2')
    clearnc.close()
    """
    filename = original_files[orig_filenames.index(base_file)]

    orignc=netCDF4.Dataset(filename)
    sst_val = read_var(orignc,'sea_surface_temperature')
    orignc.close()


    for j in points:
        
        sst[j].append(sst_val[j[0],j[1]])

        #approx[j].append(approx_val[j[0],j[1]])
          
        #smooth[j].append(smooth_val[j[0],j[1]])

        approx2[j].append(approx2_val[j[0],j[1]])

        collated[j].append(collated_vals[j[0],j[1]])
        
        collated2[j].append(collated2_vals[j[0],j[1]])

        smooth_collate[j].append(smooth_collate_val[j[0],j[1]])

        reinstated[j].append(reinstated_vals[j[0],j[1]])

        """
        if(sst_clear_val[j[0],j[1]] == 255):
        	clear[j].append(sst[j][i])
        else:
        	clear[j].append(np.nan)
        """
        if(sst_clear_pass2[j[0],j[1]]):
            second_pass[j].append(sst[j][i])
        else:
        	second_pass[j].append(np.nan)

    

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
    plt.plot(range(total_files),reinstated[i],'.', markersize=15,label="reinstated")
    #plt.plot(range(total_files),clear[i],'.', markersize=15,label="clear")
    plt.plot(range(total_files),sst[i],label="orig")
    #plt.plot(range(total_files),smooth[i],label="smooth")
    
    #plt.plot(collated_inds,interpolated[i],'.',markersize=15,label="interp")
    plt.plot(range(total_files),collated[i],label="collated")
    plt.plot(range(total_files),collated2[i],'.',markersize=15,label="collated2")
    plt.plot(hour_ind,h[i],'.',markersize=15,label="hour")
    plt.legend()
    plt.xticks(hour_ind,times,rotation='vertical')
    plt.show()
