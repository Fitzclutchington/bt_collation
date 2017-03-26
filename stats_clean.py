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
collated2 = []
scollated = []
sst_std = []
approx_std = []
approx2_std = []
clear_std =[]
#smooth_collated =[]

num_clear = []
clear_pass2 = []
clear_approx = []
clear_approx2 = []
clear_col = []
clear_col2 = []
clear_scol = []
clear_interp = []
#clear_smooth = []
acspo = []

ref_file = "/home/fitz/sst_approximation/data/ACSPO_V2.41_H08_AHI_2017-01-09_0050-0100_20170110.033449.nc"
cdf = netCDF4.Dataset(ref_file)
ref_sst = read_var(cdf,"sst_reynolds")
sza = read_var(cdf,"satellite_zenith_angle")
mask = np.abs(sza) < 67

ref_sst[~mask] = np.nan

total_files = len(approx_files)

for i in range(total_files):
	base_file = approx_files[i].split('/')[2]
	cdf = netCDF4.Dataset(approx_files[i])
	vals = read_var(cdf,"sea_surface_temperature")
	cdf.close()
	approx_std.append(np.nanstd(vals-ref_sst))
	clear_approx.append(np.isfinite(vals).sum())

	cdf = netCDF4.Dataset(approx2_files[i])
	vals = read_var(cdf,"sea_surface_temperature")
	approx2_std.append(np.nanstd(vals-ref_sst))
	clear_approx2.append(np.isfinite(vals).sum())
	cdf.close()

	#pass2_filename = pass2_files[pass2_filenames.index(approx_file.split('/')[2])]
	#cdf = netCDF4.Dataset(pass2_filename)
	#pass2 = read_var(cdf,"brightness_temperature_11um2")
	#cdf.close()

	#clear_filename = clear_files[clear_filenames.index(approx_file.split('/')[2])]
	#cdf = netCDF4.Dataset(clear_filename)
	#clear = read_var(cdf,"brightness_temperature_11um2")
	#clear_mask = clear == 0
	#num_clear.append((~clear_mask).sum())

	sst_filename = original_files[orig_filenames.index(base_file)]
	cdf = netCDF4.Dataset(sst_filename)
	l2p_flags = read_var(cdf,'l2p_flags')
	sst_acspo = read_var(cdf, 'sea_surface_temperature')
	
	sst_acspo[~mask] = np.nan
	cdf.close()

	cloud_mask = np.bitwise_and(l2p_flags,-16384).astype(bool)
	#pass2_mask = pass2 == 0
	sst_acspo[cloud_mask] = np.nan
	sst_std.append(np.nanstd(sst_acspo - ref_sst))

	#sst_pass2[~mask] = np.nan
	#sst_pass2[pass2_mask] = np.nan
	#clear_std.append(np.nanstd(sst_pass2-ref_sst))

	#clear_pass2.append(np.isfinite(sst_pass2).sum())
	acspo.append(np.isfinite(sst_acspo).sum())
	print "finished iteration for " + approx_files[i]

	cdf = netCDF4.Dataset(collated_files[i])
	col = read_var(cdf,"sea_surface_temperature")
	col[~mask] = np.nan
	collated.append(np.nanstd(col - ref_sst))
	clear_col.append(np.isfinite(col).sum())
	cdf.close()

	cdf = netCDF4.Dataset(collated2_files[i])
	col = read_var(cdf,"sea_surface_temperature")
	col[~mask] = np.nan
	collated2.append(np.nanstd(col - ref_sst))
	clear_col2.append(np.isfinite(col).sum())
	cdf.close()

	cdf = netCDF4.Dataset(scollate_files[i])
	col = read_var(cdf,"sea_surface_temperature")
	col[~mask] = np.nan
	scollated.append(np.nanstd(col - ref_sst))
	clear_scol.append(np.isfinite(col).sum())
	cdf.close()	
	


plt.figure()
plt.plot(collated_inds,collated, label="collated")
plt.plot(collated_inds,collated2, label="collated2")
plt.plot(collated_inds,scollated, label="scollated")
#plt.plot(smooth_collate_inds,smooth_collated, label="smooth collated")
plt.plot(sst_std, label="acspo")
#plt.plot(clear_std, label="clear")
plt.plot(approx_std, label="approx")
plt.plot(approx2_std, label="approx2")
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
plt.plot(collated_inds,clear_col2, label="collated2")
plt.plot(collated_inds,clear_scol, label="scollated")
#plt.plot(smooth_collate_inds,clear_smooth, label="smooth collated")
plt.plot(acspo,label="acspo")
#plt.plot(clear_pass2, label="pass2")
plt.plot(clear_approx, label="approx")
plt.plot(clear_approx2, label="approx2")
plt.legend()
plt.title("Number Clear Interpolated Collated")
plt.show()

outfile = "num_clear"
sio.savemat(outfile,{"clear":num_clear})
"""
plt.plot(collated_inds,clear_interp, label="interpolated")


"""

