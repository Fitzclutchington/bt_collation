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
interp_files = sorted(glob.glob("data/collated_interp/*.nc"))
collated_files = sorted(glob.glob("data/collated_mat/*.nc"))
pass2_files = sorted(glob.glob("data/pass2/*.nc"))
smooth_files = sorted(glob.glob("data/smooth/*.nc"))
#smooth_collate_files = sorted(glob.glob("data/smooth_collate/*.nc"))
#collated_old = sorted(glob.glob("data/collated_old_eigen/*.nc"))
original_files = [line.rstrip('\n') for line in open('ahitest.txt')]
orig_filenames = [x.split('/')[6] for x in original_files]
pass2_filenames = [x.split('/')[2] for x in pass2_files]
approx_filenames = [x.split('/')[2] for x in approx_files]
clear_filenames = [x.split('/')[2] for x in clear_files]
smooth_filenames = [x.split('/')[2] for x in smooth_files]
#smooth_collate_filenames = [x.split('/')[2] for x in smooth_collate_files]

collated = []
collated_interp = []
sst_std = []
approx_std = []
clear_std =[]
#smooth_collated =[]

num_clear = []
clear_pass2 = []
clear_approx = []
clear_col = []
clear_interp = []
#clear_smooth = []
acspo = []

ref_file = "/home/fitz/sst_approximation/data/ACSPO_V2.41_H08_AHI_2017-01-09_0050-0100_20170110.033449.nc"
cdf = netCDF4.Dataset(ref_file)
ref_sst = read_var(cdf,"sst_reynolds")
sza = read_var(cdf,"satellite_zenith_angle")
mask = np.abs(sza) < 67

ref_sst[~mask] = np.nan
#masked_sst = ref_sst.copy()
#masked_sst[~mask] = np.nan
"""
plt.figure()
#cmap = cmocean.cm.thermal

ax1 = plt.subplot(121)
img1 = ax1.imshow(masked_sst,vmin=270,vmax=307)
#ax1.imshow(cloud,interpolation='nearest')
#ax1.imshow(land,interpolation='nearest')
div1 = make_axes_locatable(ax1)
cax1 = div1.append_axes("right", size="5%", pad=0.05)
cbar1 = plt.colorbar(img1, cax=cax1)


ax2 = plt.subplot(122, sharex=ax1,sharey=ax1)
img2 = ax2.imshow(ref_sst,vmin=-6,vmax=6)
#ax2.imshow(cloud,interpolation='nearest')
#ax2.imshow(land,interpolation='nearest')
div2 = make_axes_locatable(ax2)
cax2 = div2.append_axes("right", size="5%", pad=0.05)
cbar2 = plt.colorbar(img2, cax=cax2)
"""
"""
o2ld = read_var(cdf,"ocean_to_land_dist")
good_mask = np.isfinite(o2ld)

filename = "data/collated_mat/20170108120000-STAR-L2P_GHRSST-SSTsubskin-AHI_H08-ACSPO_V2.41-v02.0-fv01.0.nc"
cdf = netCDF4.Dataset(filename)
sst = read_var(cdf,"sea_surface_temperature")

sst[~mask] = np.nan
vals = np.nanstd(sst - ref_sst)

#plt.figure()
#plt.hist(o2ld[np.logical_and(m1,good_mask)])
#plt.show()
#m2 = vals < -2
#full_mask = np.logical_or(m1,m2)
vals[~m1] = np.nan


plt.figure()
#cmap = cmocean.cm.thermal

ax1 = plt.subplot(121)
img1 = ax1.imshow(sst,vmin=270,vmax=307)
#ax1.imshow(cloud,interpolation='nearest')
ax1.imshow(land,interpolation='nearest')
div1 = make_axes_locatable(ax1)
cax1 = div1.append_axes("right", size="5%", pad=0.05)
cbar1 = plt.colorbar(img1, cax=cax1)


ax2 = plt.subplot(122, sharex=ax1,sharey=ax1)
img2 = ax2.imshow(vals,vmin=-6,vmax=6)
#ax2.imshow(cloud,interpolation='nearest')
ax2.imshow(land,interpolation='nearest')
div2 = make_axes_locatable(ax2)
cax2 = div2.append_axes("right", size="5%", pad=0.05)
cbar2 = plt.colorbar(img2, cax=cax2)
"""
for approx_file in approx_files:
	cdf = netCDF4.Dataset(approx_file)
	approx = read_var(cdf,"sea_surface_temperature")

	cdf.close()

	approx_std.append(np.nanstd(approx-ref_sst))
	clear_approx.append(np.isfinite(approx).sum())

	pass2_filename = pass2_files[pass2_filenames.index(approx_file.split('/')[2])]
	cdf = netCDF4.Dataset(pass2_filename)
	pass2 = read_var(cdf,"brightness_temperature_11um2")
	cdf.close()

	#clear_filename = clear_files[clear_filenames.index(approx_file.split('/')[2])]
	#cdf = netCDF4.Dataset(clear_filename)
	#clear = read_var(cdf,"brightness_temperature_11um2")
	#clear_mask = clear == 0
	#num_clear.append((~clear_mask).sum())

	sst_filename = original_files[orig_filenames.index(approx_file.split('/')[2])]
	print sst_filename
	cdf = netCDF4.Dataset(sst_filename)
	l2p_flags = read_var(cdf,'l2p_flags')
	sst_pass2 = read_var(cdf, 'sea_surface_temperature')
	
	sst_acspo = sst_pass2.copy()
	sst_acspo[~mask] = np.nan
	cdf.close()

	cloud_mask = np.bitwise_and(l2p_flags,-16384).astype(bool)
	pass2_mask = pass2 == 0
	sst_acspo[cloud_mask] = np.nan
	sst_std.append(np.nanstd(sst_acspo - ref_sst))

	sst_pass2[~mask] = np.nan
	sst_pass2[pass2_mask] = np.nan
	clear_std.append(np.nanstd(sst_pass2-ref_sst))

	clear_pass2.append(np.isfinite(sst_pass2).sum())
	acspo.append(np.isfinite(sst_acspo).sum())
	print "finished iteration for " + approx_file

collated_inds = []	
for f1 in collated_files:
	cdf = netCDF4.Dataset(f1)
	col = read_var(cdf,"sea_surface_temperature")
	col[~mask] = np.nan
	collated.append(np.nanstd(col - ref_sst))
	cdf.close()

	clear_col.append(np.isfinite(col).sum())
	collated_inds.append(approx_filenames.index(f1.split('/')[2]))

"""
smooth_collate_inds = []	
for f1 in smooth_collate_files:
	cdf = netCDF4.Dataset(f1)
	col = read_var(cdf,"sea_surface_temperature")
	col[~mask] = np.nan
	smooth_collated.append(np.nanstd(col - ref_sst))
	cdf.close()

	clear_smooth.append(np.isfinite(col).sum())
	smooth_collate_inds.append(approx_filenames.index(f1.split('/')[2]))
"""

plt.figure()
plt.plot(collated_inds,collated, label="collated")
#plt.plot(smooth_collate_inds,smooth_collated, label="smooth collated")
plt.plot(sst_std, label="acspo")
plt.plot(clear_std, label="clear")
plt.plot(approx_std, label="approx")
plt.legend()
plt.title("Standard Deviation Interpolated Collated")
plt.show()

"""
plt.plot(collated_inds,collated_interp, label="interpolated")
plt.plot(sst_std, label="acspo")
plt.plot(clear_std, label="clear")
plt.plot(approx_std, label="approx")
plt.legend()
plt.title("Standard Deviation Interpolated Collated")
plt.show()
"""
plt.figure()
plt.plot(collated_inds,clear_col, label="collated")
#plt.plot(smooth_collate_inds,clear_smooth, label="smooth collated")
plt.plot(acspo,label="acspo")
plt.plot(clear_pass2, label="pass2")
plt.plot(clear_approx, label="approx")
plt.legend()
plt.title("Number Clear Interpolated Collated")
plt.show()

outfile = "num_clear"
sio.savemat(outfile,{"clear":num_clear})
"""
plt.plot(collated_inds,clear_interp, label="interpolated")


"""

