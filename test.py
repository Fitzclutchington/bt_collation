import netCDF4
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import glob
import scipy.io as sio
import cmocean

def read_var(cdf,variable):
    data = np.squeeze(cdf[variable][:])
    return data


original_file = '../sst_approximation/data/2017-01-08/20170108120000-STAR-L2P_GHRSST-SSTsubskin-AHI_H08-ACSPO_V2.41-v02.0-fv01.0.nc'
pass1_file = 'data/clear/20170108120000-STAR-L2P_GHRSST-SSTsubskin-AHI_H08-ACSPO_V2.41-v02.0-fv01.0.nc'
pass2_file = 'data/pass2/20170108120000-STAR-L2P_GHRSST-SSTsubskin-AHI_H08-ACSPO_V2.41-v02.0-fv01.0.nc'
smooth_file = 'data/smooth/20170108120000-STAR-L2P_GHRSST-SSTsubskin-AHI_H08-ACSPO_V2.41-v02.0-fv01.0.nc'
smooth_collate_file = 'data/smooth_collate/20170108120000-STAR-L2P_GHRSST-SSTsubskin-AHI_H08-ACSPO_V2.41-v02.0-fv01.0.nc'
approx_file = 'data/approx/20170108120000-STAR-L2P_GHRSST-SSTsubskin-AHI_H08-ACSPO_V2.41-v02.0-fv01.0.nc'
approx2_file = 'data/approx2/20170108120000-STAR-L2P_GHRSST-SSTsubskin-AHI_H08-ACSPO_V2.41-v02.0-fv01.0.nc'
collated_file = 'data/collated_mat/20170108120000-STAR-L2P_GHRSST-SSTsubskin-AHI_H08-ACSPO_V2.41-v02.0-fv01.0.nc'
collated2_file = 'data/collated_mat2/20170108120000-STAR-L2P_GHRSST-SSTsubskin-AHI_H08-ACSPO_V2.41-v02.0-fv01.0.nc'
reinstated_file = 'data/reinstated/20170108120000-STAR-L2P_GHRSST-SSTsubskin-AHI_H08-ACSPO_V2.41-v02.0-fv01.0.nc'

files = [original_file, smooth_file,smooth_collate_file,approx_file,approx2_file,collated_file,collated2_file, reinstated_file]
titles = ['Original SST', 'SST Moving Average 101 X 101 X 19', 'SST Moving Average (10 x 10 x 7)','SST Approximation 1st Iteration (Least Squares)','Approximation 2nd Iteration (Least Squares)','High Resolution Reference','SST Collation', "Reinstated Clear SST"]

mask_files = [pass1_file,pass2_file]
mask_titles = ['Masking Sequence (1st Pass)', 'Clear After Second Pass']
cdf = netCDF4.Dataset(original_file)
sst = read_var(cdf,'sea_surface_temperature')

l2p_flags = read_var(cdf,'l2p_flags')
land_mask = np.bitwise_and(l2p_flags,2).astype(bool)
land = np.zeros((5500,5500,4))
r = 146/256.0
g = 98/256.0
b = 57/256.0
land[land_mask] = [r,g,b,1]

cloud_mask = np.bitwise_and(l2p_flags,-16384).astype(bool)

date_time = "12:00 UTC 01-08-2017\n"
for f,title in zip(files,titles):
	print f
	cdf = netCDF4.Dataset(f)
	vals = read_var(cdf,'sea_surface_temperature')

	plt.figure()
	plt.title(date_time + title)
	plt.imshow(vals,vmin=270,vmax=307)
	plt.colorbar()
	plt.imshow(land)	
	plt.show()


for f,title in zip(mask_files, mask_titles):
	cdf = netCDF4.Dataset(f)
	mask = read_var(cdf,'brightness_temperature_11um2')
	mask = mask == 0
	clear_sst = sst.copy()
	clear_sst[mask] = np.nan

	plt.figure()
	plt.title(date_time + title)
	plt.imshow(clear_sst,vmin=270,vmax=307)
	plt.colorbar()
	plt.imshow(land)	
	plt.show()

clear_sst = sst.copy()
clear_sst[cloud_mask] = np.nan

plt.figure()
plt.title(date_time + "SST ACSPO Cloud Flags")
plt.imshow(clear_sst,vmin=270,vmax=307)
plt.colorbar()
plt.imshow(land)	
plt.show()
