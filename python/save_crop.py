import netCDF4
import numpy as np
import glob
import sys

def read_var(cdf,variable):
    data = np.squeeze(cdf[variable][:])
    return data

if len(sys.argv) < 8:
	print "usage: python save_crop.py <min_x> <max_x> <min_y> <max_y> <variable> <folder> <output_folder>"
	

max_x = int(sys.argv[2])
min_x = int(sys.argv[1])-1
max_y = int(sys.argv[4])
min_y = int(sys.argv[3])-1

folder_name = sys.argv[6]
output_folder = sys.argv[7]
variable_name = sys.argv[5]
	
print folder_name
files = sorted(glob.glob(folder_name+"/*.nc"))

filenames = [x.split('/')[-1] for x in files]

for i,f in enumerate(files):

	cdf = netCDF4.Dataset(f)
	data = read_var(cdf,variable_name)[min_x:max_x,min_y:max_y]
	cdf.close()
	
	save_loc = output_folder + '/' + filenames[i]
	rootgrp = netCDF4.Dataset(save_loc, "w", format="NETCDF4")
	height = rootgrp.createDimension("height", max_x-min_x)
	width = rootgrp.createDimension("width", max_y-min_y)

	var = rootgrp.createVariable(variable_name,"f4",("height","width"))
	var[:] = data
	rootgrp.close()
	print f