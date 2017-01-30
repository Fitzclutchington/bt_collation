import netCDF4
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import glob

def read_var(cdf,variable):
    data = np.squeeze(cdf[variable][:])
    return data


interp_files = sorted(glob.glob("data/collated_interp/*.nc"))
collated_files = sorted(glob.glob("data/collated/*.nc"))
original_files = [line.rstrip('\n') for line in open('ahitest.txt')]
filenames = [x.split('/')[6] for x in original_files]

collated = []
collated_interp = []
sst_std = []

clear_col = []
clear_interp = []
acspo = []

ref_file = "../data/ACSPO_V2.41_H08_AHI_2016-09-29_0010-0020_20160930.022950.nc"
cdf = netCDF4.Dataset(ref_file)
ref_sst = read_var(cdf,"sst_reynolds")
sza = read_var(cdf,"satellite_zenith_angle")
mask = np.abs(sza) < 67
ref_sst[~mask] = np.nan

for f1,f2 in zip(collated_files,interp_files):
	cdf = netCDF4.Dataset(f1)
	col = read_var(cdf,"sea_surface_temperature")
	col[~mask] = np.nan
	collated.append(np.nanstd(col - ref_sst))

	cdf = netCDF4.Dataset(f2)
	col_interp = read_var(cdf,"sea_surface_temperature")
	col_interp[~mask] = np.nan
	collated_interp.append(np.nanstd(col_interp - ref_sst))

	clear_col.append(np.isfinite(col).sum())
	clear_interp.append(np.isfinite(col_interp).sum())

	sst_filename = original_files[filenames.index(f1.split('/')[2])]
	cdf = netCDF4.Dataset(sst_filename)
	l2p_flags = read_var(cdf,'l2p_flags')
	sst_orig = read_var(cdf, 'sea_surface_temperature')
	cloud_mask = np.bitwise_and(l2p_flags,-16384).astype(bool)

	sst_orig[cloud_mask] = np.nan
	sst_orig[~mask] = np.nan
	sst_std.append(np.nanstd(sst_orig - ref_sst))

	acspo.append((~cloud_mask).sum())


plt.figure()
plt.plot(collated, label="collated")
plt.plot(collated_interp, label="interpolated")
plt.plot(sst_std, label="acspo")
plt.legend()
plt.title("Standard Deviation Interpolated Collated")
plt.show()

plt.figure()
plt.plot(clear_col, label="collated")
plt.plot(clear_interp, label="interpolated")
plt.plot(acspo,label="acspo")
plt.legend()
plt.title("Number Clear Interpolated Collated")
plt.show()


collated_filename = collated_files[15]
interp_filename = interp_files[15]
sst_filename = original_files[filenames.index(collated_files[15].split('/')[2])]

cdf = netCDF4.Dataset(sst_filename)

sst = read_var(cdf,"sea_surface_temperature")

l2p_flags = read_var(cdf,'l2p_flags')
land_mask = np.bitwise_and(l2p_flags,2).astype(bool)
land = np.zeros((5500,5500,4))
r = 146/256.0
g = 98/256.0
b = 57/256.0
land[land_mask] = [r,g,b,1]

cdf = netCDF4.Dataset(collated_filename)
collated = read_var(cdf, "sea_surface_temperature")

cdf = netCDF4.Dataset(interp_filename)
interp = read_var(cdf,"sea_surface_temperature")

plt.figure()


ax1 = plt.subplot(131)
img1 = ax1.imshow(collated,vmin=270,vmax=307)
ax1.imshow(land,interpolation='nearest')
div1 = make_axes_locatable(ax1)
cax1 = div1.append_axes("right", size="5%", pad=0.05)
cbar1 = plt.colorbar(img1, cax=cax1)

ax2 = plt.subplot(132,sharex=ax1,sharey=ax1)
img2 = ax2.imshow(interp,vmin=270,vmax=307)
ax2.imshow(land,interpolation='nearest')
div2 = make_axes_locatable(ax2)
cax2 = div2.append_axes("right", size="5%", pad=0.05)
cbar2 = plt.colorbar(img2, cax=cax2)

ax3 = plt.subplot(133,sharex=ax1,sharey=ax1)
img3 = ax3.imshow(sst,vmin=270,vmax=307)
ax3.imshow(land,interpolation='nearest')
div3 = make_axes_locatable(ax3)
cax3 = div3.append_axes("right", size="5%", pad=0.05)
cbar3 = plt.colorbar(img3, cax=cax3)

"""
ax3 = plt.subplot(133,sharex=ax1,sharey=ax1)
img3 = ax3.imshow(means,vmin=270,vmax=307)
ax3.imshow(land,interpolation='nearest')
div3 = make_axes_locatable(ax3)
cax3 = div3.append_axes("right", size="5%", pad=0.05)
cbar3 = plt.colorbar(img3, cax=cax3)
"""
"""
plt.figure()
ax1 = plt.subplot(121)
img1 = ax1.imshow(means-mins,vmin=0,vmax=5)
ax1.imshow(land,interpolation='nearest')
div1 = make_axes_locatable(ax1)
cax1 = div1.append_axes("right", size="5%", pad=0.05)
cbar1 = plt.colorbar(img1, cax=cax1)

ax2 = plt.subplot(122,sharex=ax1,sharey=ax1)
img2 = ax2.imshow(maxes-means,vmin=0,vmax=5)
ax2.imshow(land,interpolation='nearest')
div2 = make_axes_locatable(ax2)
cax2 = div2.append_axes("right", size="5%", pad=0.05)
cbar2 = plt.colorbar(img2, cax=cax2)
"""
