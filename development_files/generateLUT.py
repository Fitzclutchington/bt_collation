import numpy as np
import netCDF4
import sys

if  len(sys.argv) < 2:
	print "usage: python generateLUT.py <reference_file>"

#open reference file
reference_file = sys.argv[1]
ref_cdf = netCDF4.Dataset(reference_file)

sst_ref = np.squeeze(ref_cdf['sst_reynolds'])
sza = np.squeeze(ref_cdf['solar_zenith_angle'])
lons = np.squeeze(ref_cdf['longitude'])
negative_lon_mask = lons < 0
lons[negative_lon_mask] += 360

dt = 1
temp_min = 270
temp_max = 310
temp = np.arange(temp_min,temp_max+1,1)