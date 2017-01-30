import netCDF4
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import glob

def read_var(cdf,variable):
    data = np.squeeze(cdf[variable][:])
    return data

collated_files = sorted(glob.glob("data/collated/*.nc"))

reference_file = "../data/ACSPO_V2.41_H08_AHI_2016-11-15_0300-0310_20161116.034934.nc"
ref_cdf = netCDF4.Dataset(reference_file)

sst_ref = read_var(ref_cdf,"sst_reynolds")
sza = np.squeeze(ref_cdf['solar_zenith_angle'])

lons = read_var(ref_cdf,"longitude")

valid_mask_67 = np.abs(sza) <= 67
valid_mask_60 = np.abs(sza) <= 60

stds_67 = {}
means_67 = {}
#stds_60 = []
#means_60 = []
#mins = []
#maxs = []

means = np.zeros((len(collated_files),len(zip(range(-185,181,5),range(-180,176,5)))))

#for i,j in zip(range(-185,181,5),range(-180,176,5)):
#    stds_67[(i,j)] = []
#    means_67[(i,j)] = []


for y,collated_file in enumerate(collated_files):
    collated_cdf = netCDF4.Dataset(collated_file)
    sst_collated = read_var(collated_cdf,"sea_surface_temperature")

    diffs = sst_collated - sst_ref

    for x,(i,j) in enumerate(zip(range(-185,181,5),range(-180,176,5))):
        mask1 = lons <= j
        mask2 = lons >= i
        total_mask = np.logical_and(valid_mask_67,np.logical_and(mask1,mask2))
        
        #stds_67[(i,j)].append(np.nanstd(diffs[total_mask]))
        means[y,x] = np.nanmean(diffs[total_mask])
    #stds_60.append(np.nanstd(diffs[valid_mask_60]))
    #means_60.append(np.nanmean(diffs[valid_mask_60]))
    #mins.append(np.nanmin(diffs[valid_mask]))
    #maxs.append(np.nanmax(diffs[valid_mask]))


plt.figure()
plt.imshow(means, vmin=-.5,vmax=1.5)
plt.colorbar()
for i,j in zip(range(-185,181,5),range(-180,176,5)):
    #plt.plot(stds_67[(i,j)],label="stds"+str((i,j)))
    plt.plot(means_67[(i,j)],label="means"+str((i,j)))

#plt.plot(maxs,label="maxs")
#plt.plot(mins,label="mins")
plt.legend()


collated_cdf = netCDF4.Dataset(collated_files[19])
sst_collated = read_var(collated_cdf,"sea_surface_temperature")
diffs = sst_collated - sst_ref
mask = diffs < -5
mask = np.logical_and(mask,valid_mask)
count = mask.astype(int).sum()
diffs[~mask] = NAN
plt.figure()
plt.imshow(diffs)
plt.colorbar()



plt.figure()
plt.imshow(diffs,vmin=-5,vmax=1)
plt.colorbar()