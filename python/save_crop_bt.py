import netCDF4
import numpy as np
import glob
import sys

def read_var(cdf,variable):
    data = np.squeeze(cdf[variable][:])
    return data

if len(sys.argv) < 3:
	print "usage: python save_crop.py <output_folder> <file_list>"
	print "date as yyyy-mm-dd"
	

min_ys = [4600, 4150, 3500, 3500, 1300, 500]
max_ys = [5101, 4651, 4501, 4001, 1801, 1001]
min_xs = [2700, 3250, 1200, 3000, 1200, 2300]
max_xs = [3701, 4251, 1701, 4001, 2201, 3301]

output_folder = sys.argv[1].strip('/')


file_list = sys.argv[2]
data_folder = '../data'
output_folder = '/'.join([data_folder,output_folder])

x_lims = []
original_files = [line.rstrip('\n') for line in open(file_list)]
orig_filenames = [x.split('/')[-1] for x in original_files]

for i,f in enumerate(orig_filenames):
    if f[8:12] == '0000':
        x_lims.append(i)

filenames = orig_filenames[x_lims[0]:x_lims[1]+1]
original_files = original_files[x_lims[0]:x_lims[1]+1]
orig_filenames = orig_filenames[x_lims[0]:x_lims[1]+1]

bt08_folder = 'bt08'
bt10_folder = 'bt10'
bt11_folder = 'bt11'
bt12_folder = 'bt12'

variables = ['brightness_temperature_08um6', 'brightness_temperature_10um4','brightness_temperature_11um2','brightness_temperature_12um3']

for j in range(len(max_xs)):
	max_y = max_ys[j]
	min_y = min_ys[j]
	max_x = max_xs[j]
	min_x = min_xs[j]

	crop_folder = '_'.join(map(str,[min_y,max_y,min_x,max_x]))
	print crop_folder
	w = max_x-min_x
	h = max_y-min_y

	# save sst and l2p flags
	variable_name = "sea_surface_temperature"
	
	for i,f in enumerate(original_files):

		f_name = orig_filenames[i]
		cdf = netCDF4.Dataset(f)

		for variable in variables:
			save_folder = '/'.join([output_folder, crop_folder])
			save_loc = '/'.join([save_folder , variable,f_name])
			print save_loc

			rootgrp = netCDF4.Dataset(save_loc, "w", format="NETCDF4")
			height = rootgrp.createDimension("height", h)
			width = rootgrp.createDimension("width", w)

		

			data = read_var(cdf,variable)[min_y:max_y,min_x:max_x]

			var = rootgrp.createVariable(variable,"f4",("height","width"))
			var[:] = data

		cdf.close()			
		rootgrp.close()

		print save_loc