import netCDF4
import numpy as np
import glob
import sys
from matplotlib import pyplot as plt

def read_var(cdf,variable):
    data = np.squeeze(cdf[variable][:])
    return data

if len(sys.argv) < 8:
	print "usage: python save_crop.py <min_x> <max_x> <min_y> <max_y>"
	

max_x = int(sys.argv[2])
min_x = int(sys.argv[1])-1
max_y = int(sys.argv[4])
min_y = int(sys.argv[3])-1

variable_name = 'sea_surface_temperature'

original_filename = "/home/fitz/sst_approximation/data/2017-01-08/20170108120000-STAR-L2P_GHRSST-SSTsubskin-AHI_H08-ACSPO_V2.41-v02.0-fv01.0.nc"
base_filename = original_filename.split('/')[6]

folders = ["smooth_test","smooth_collate","approx","approx2","collated_mat","collated_mat2","reinstated"]
masks = ["clear","pass2"]

titles = ['SST Moving Average 101 X 101 X 19', 'SST Moving Average (11 x 11 x 7)','SST Approximation 1st Iteration (Least Squares)','Approximation 2nd Iteration (Least Squares)','High Resolution Reference','SST Collation', "Reinstated Clear SST"]
mask_titles = ['Masking Sequence (1st Pass)', 'Clear After Second Pass']

output_folder = str(min_y)+"_"+str(max_y)+"_"+str(min_x)+"_"+str(max_x)

cdf = netCDF4.Dataset(original_filename)

l2p_flags = read_var(cdf,'l2p_flags')
land_mask = np.bitwise_and(l2p_flags,2).astype(bool)
land = np.zeros((5500,5500,4))
r = 146/256.0
g = 98/256.0
b = 57/256.0
land[land_mask] = [r,g,b,1]

land = land[min_y:max_y,min_x:max_x]

cdf = netCDF4.Dataset(original_filename)
sst_data = read_var(cdf,variable_name)[min_y:max_y,min_x:max_x]
cdf.close()

plt.figure()
plt.imshow(sst_data,vmin=270,vmax=307)
plt.colorbar()
plt.imshow(land)
plt.title("Original SST")
plt.show()


for f,title in zip(folders,titles):

	filename = '/'.join(["../data",f,base_filename])

	cdf = netCDF4.Dataset(filename)
	data = read_var(cdf,variable_name)[min_y:max_y,min_x:max_x]
	cdf.close()
	
	plt.figure()
	plt.imshow(data,vmin=270,vmax=307)
	plt.colorbar()
	plt.imshow(land)
	plt.title(title)
	plt.show()

for f,mask_title in zip(masks,mask_titles):

	filename = '/'.join(["../data",f,base_filename])

	cdf = netCDF4.Dataset(filename)
	data = read_var(cdf,"brightness_temperature_11um2")[min_y:max_y,min_x:max_x]
	cdf.close()
	
	m = data == 0
	clear_sst = sst_data.copy()
	clear_sst[m] = np.nan

	plt.figure()
	plt.imshow(clear_sst,vmin=270,vmax=307)
	plt.colorbar()
	plt.imshow(land)
	plt.title(mask_title)
	plt.show()
