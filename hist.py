import netCDF4
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import glob

def read_var(cdf,variable):
    data = np.squeeze(cdf[variable][:])
    return data

pass2_files = sorted(glob.glob("data/pass2/*.nc"))

pass2_filenames = [x.split('/')[2] for x in pass2_files]

filename = 'data/pass2/20170108070000-STAR-L2P_GHRSST-SSTsubskin-AHI_H08-ACSPO_V2.41-v02.0-fv01.0.nc'

cdf = netCDF4.Dataset(filename)
mask = read_var(cdf,"brightness_temperature_11um2")

filename = '../sst_approximation/data/2017-01-08/20170108070000-STAR-L2P_GHRSST-SSTsubskin-AHI_H08-ACSPO_V2.41-v02.0-fv01.0.nc'
cdf = netCDF4.Dataset(filename)
sst = read_var(cdf,"sea_surface_temperature").data

sst[mask == 0] = np.nan
sst[sst <  0] = np.nan

ref_file = "/home/fitz/sst_approximation/data/ACSPO_V2.41_H08_AHI_2017-01-09_0050-0100_20170110.033449.nc"
cdf = netCDF4.Dataset(ref_file)
ref_sst = read_var(cdf,"sst_reynolds")
sza = read_var(cdf,"satellite_zenith_angle")
mask = np.abs(sza) < 67

ref_sst[~mask] = np.nan

diffs = sst - ref_sst


vals = diffs[np.isfinite(diffs)]

plt.figure()
plt.hist(vals,bins=50)
plt.show()

