import netCDF4
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import glob
import scipy.io as sio

def read_var(cdf,variable):
    data = np.squeeze(cdf[variable][:])
    return data

cold_files = sorted(glob.glob("data/cold_mask/*.nc"))
nn_files = sorted(glob.glob("data/nn_mask/*.nc"))
diag_files = sorted(glob.glob("data/diag_mask/*.nc"))
ratio_files = sorted(glob.glob("data/ratio_mask/*.nc"))
dt_files = sorted(glob.glob("data/dt_mask/*.nc"))
clear_files = sorted(glob.glob("data/clear/*.nc"))
approx_files = sorted(glob.glob("data/approx/*.nc"))
original_files = [line.rstrip('\n') for line in open('ahitest.txt')]
orig_filenames = [x.split('/')[6] for x in original_files]
clear_filenames = [x.split('/')[2] for x in clear_files]

clear_old = np.squeeze(sio.loadmat("num_clear.mat")['clear'])
num_clear = []

cold_mask_clear = []
nn_mask_clear = []
diag_mask_clear = []
ratio_mask_clear = []
dt_mask_clear = []

for i,f in enumerate(cold_files):
	clear_filename = cold_files[i]
	cdf = netCDF4.Dataset(clear_filename)
	clear = read_var(cdf,"brightness_temperature_11um2")
	mask = clear == 255
	cold_mask_clear.append(mask.sum())

	clear_filename = nn_files[i]
	cdf = netCDF4.Dataset(clear_filename)
	clear = read_var(cdf,"brightness_temperature_11um2")
	mask = clear == 255
	nn_mask_clear.append(mask.sum())

	clear_filename = diag_files[i]
	cdf = netCDF4.Dataset(clear_filename)
	clear = read_var(cdf,"brightness_temperature_11um2")
	mask = clear == 255
	diag_mask_clear.append(mask.sum())

	clear_filename = ratio_files[i]
	cdf = netCDF4.Dataset(clear_filename)
	clear = read_var(cdf,"brightness_temperature_11um2")
	mask = clear == 255
	ratio_mask_clear.append(mask.sum())

	clear_filename = dt_files[i]
	cdf = netCDF4.Dataset(clear_filename)
	clear = read_var(cdf,"brightness_temperature_11um2")
	mask = clear == 255
	dt_mask_clear.append(mask.sum())


plt.figure()
plt.title("Number of Observations per Mask")
plt.plot(diag_mask_clear,label="2nd Order Diff")
plt.plot(nn_mask_clear,label="Diff")
plt.plot(cold_mask_clear,label="BT08 BT10 Test")
plt.plot(ratio_mask_clear,label="BT08 BT11 BT12 Test")
plt.plot(dt_mask_clear,label="dt analysis")
plt.legend()
plt.show()

for f in approx_files:
	clear_filename = clear_files[clear_filenames.index(f.split('/')[2])]
	cdf = netCDF4.Dataset(clear_filename)
	clear = read_var(cdf,"brightness_temperature_11um2")
	clear_mask = clear == 255
	num_clear.append((clear_mask).sum())

plt.figure()
plt.title("Comparison of 1st Pass With and Without Noise Mask")
plt.plot(num_clear,label="clear w/o noise mask")
plt.plot(clear_old, label="clear w/ noise mask")
plt.legend()
plt.show()

plt.figure()
plt.title("Difference Between 1st Pass With and Without Noise Mask")
plt.plot(num_clear - clear_old)
plt.show()

clear = []
sst = []
for f in clear_files:
	orig_filename = original_files[orig_filenames.index(f.split('/')[2])]
	cdf = netCDF4.Dataset(orig_filename)
	sst_val = read_var(cdf,"sea_surface_temperature")[2639,4975]
	sst.append(sst_val)

	
	cdf = netCDF4.Dataset(f)
	clear_pix = read_var(cdf,"brightness_temperature_11um2")[2639,4975]
	if clear_pix == 255:
		clear.append(sst_val)
	else:
		clear.append(np.nan)
	
plt.figure()
plt.plot(sst,label='sst')
plt.plot(clear,'r.',markersize=15,label="clear")
plt.legend()
plt.show()
	