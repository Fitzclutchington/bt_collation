#REFERENCE_SST IS DATE OF REFERENCE

import netCDF4
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import glob

def read_var(cdf,variable):
    data = np.squeeze(cdf[variable][:])
    return data

clear_files = sorted(glob.glob("data/clear/*.nc"))
approx_files = sorted(glob.glob("data/approx/*.nc"))
collated_files = sorted(glob.glob("data/collated_mat/*.nc"))
approx2_files = sorted(glob.glob("data/approx2/*.nc"))
collated2_files = sorted(glob.glob("data/collated_mat2/*.nc"))
pass2_files = sorted(glob.glob("data/pass2/*.nc"))
smooth_files = sorted(glob.glob("data/smooth_test/*.nc"))
scollate_files = sorted(glob.glob("data/smooth_collate/*.nc"))
reinstated_files = sorted(glob.glob("data/reinstated/*.nc"))

#smooth_collate_files = sorted(glob.glob("data/smooth_collate/*.nc"))
#collated_old = sorted(glob.glob("data/collated_old_eigen/*.nc"))
original_files = [line.rstrip('\n') for line in open('ahitest.txt')]
orig_filenames = [x.split('/')[6] for x in original_files]
pass2_filenames = [x.split('/')[2] for x in pass2_files]
approx_filenames = [x.split('/')[2] for x in approx_files]
clear_filenames = [x.split('/')[2] for x in clear_files]
smooth_filenames = [x.split('/')[2] for x in smooth_files]
#smooth_collate_filenames = [x.split('/')[2] for x in smooth_collate_files]

collated_min = []
collated2_min = []
sst_min = []
approx_min = []
approx2_min = []
pass2_min =[]
reinstated_min = []

collated_max = []
collated2_max = []
sst_max = []
approx_max = []
approx2_max = []
pass2_max =[]
reinstated_max = []


ref_file = "/home/fitz/sst_approximation/data/ACSPO_V2.41_H08_AHI_2017-01-09_0050-0100_20170110.033449.nc"
cdf = netCDF4.Dataset(ref_file)
ref_sst = read_var(cdf,"sst_reynolds")
sza = read_var(cdf,"satellite_zenith_angle")
mask = np.abs(sza) < 67

ref_sst[~mask] = np.nan

total_files = len(approx_files)

times =[]

hour_ind = []

for i,approx_file in enumerate(approx_filenames):
    if approx_file[10:12] == '00':
        times.append(approx_file[8:12])
        hour_ind.append(i)

for i in range(total_files):
    base_file = approx_files[i].split('/')[2]
    cdf = netCDF4.Dataset(approx_files[i])
    vals = read_var(cdf,"sea_surface_temperature")
   
    cdf.close()
    diff = vals-ref_sst
    approx_min.append(np.nanmin(diff))
    approx_max.append(np.nanmax(diff))

    cdf = netCDF4.Dataset(approx2_files[i])
    vals = read_var(cdf,"sea_surface_temperature")
    diff = vals-ref_sst
    approx2_min.append(np.nanmin(diff))
    approx2_max.append(np.nanmax(diff))
    cdf.close()

    pass2_filename = pass2_files[pass2_filenames.index(base_file)]
    cdf = netCDF4.Dataset(pass2_filename)
    pass2 = read_var(cdf,"brightness_temperature_11um2")
    cdf.close()


    sst_filename = original_files[orig_filenames.index(base_file)]
    cdf = netCDF4.Dataset(sst_filename)
    l2p_flags = read_var(cdf,'l2p_flags')
    sst_clear = read_var(cdf, 'sea_surface_temperature')
    
    sst_acspo = sst_clear.copy()
    sst_acspo[~mask] = np.nan
    cdf.close()

    cloud_mask = np.bitwise_and(l2p_flags,-16384).astype(bool)
    pass2_mask = pass2 == 0
    sst_acspo[cloud_mask] = np.nan

    diff = sst_acspo - ref_sst
    sst_min.append(np.nanmin(diff))
    sst_max.append(np.nanmax(diff))

    sst_pass2 = sst_clear.copy()
    sst_pass2[~mask] = np.nan
    sst_pass2[pass2_mask] = np.nan

    
    diff = sst_pass2-ref_sst
    pass2_min.append(np.nanmin(diff))
    pass2_max.append(np.nanmax(diff))


    print "finished iteration for " + approx_files[i]

    cdf = netCDF4.Dataset(collated_files[i])
    col = read_var(cdf,"sea_surface_temperature")
    col[~mask] = np.nan
    
    diff = col - ref_sst
    collated_min.append(np.nanmin(diff))
    collated_max.append(np.nanmax(diff))
   
    cdf.close()

    cdf = netCDF4.Dataset(collated2_files[i])
    col = read_var(cdf,"sea_surface_temperature")
    col[~mask] = np.nan
    
    diff = col - ref_sst
    collated2_min.append(np.nanmin(diff))
    collated2_max.append(np.nanmax(diff))
    cdf.close()

    
    cdf = netCDF4.Dataset(reinstated_files[i])
    col = read_var(cdf,"sea_surface_temperature")
    col[~mask] = np.nan
    
    diff = col - ref_sst
    reinstated_min.append(np.nanmin(diff))
    reinstated_max.append(np.nanmax(diff))
    cdf.close() 


plt.figure()
plt.plot(collated_min, label="collated")
plt.plot(collated2_min, label="collated2")
plt.plot(sst_min, label="acspo")
plt.plot(pass2_min, label="clear2")
plt.plot(reinstated_min, label="reinstated")
plt.plot(approx_min, label="approx")
plt.plot(approx2_min, label="approx2")
plt.legend()
plt.xticks(hour_ind,times,rotation='vertical')
plt.title("Minimums")
plt.show()


plt.figure()
plt.plot(collated_max, label="collated")
plt.plot(collated2_max, label="collated2")
plt.plot(sst_max, label="acspo")
plt.plot(pass2_max, label="clear2")
plt.plot(reinstated_max, label="reinstated")
plt.plot(approx_max, label="approx")
plt.plot(approx2_max, label="approx2")
plt.legend()
plt.xticks(hour_ind,times,rotation='vertical')
plt.title("Maximums")
plt.show()



