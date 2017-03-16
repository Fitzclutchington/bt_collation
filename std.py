import netCDF4
import numpy as np 

def read_var(cdf,variable):
    data = np.squeeze(cdf[variable][:])
    return data

ref_file = "/home/fitz/sst_approximation/data/ACSPO_V2.41_H08_AHI_2017-01-08_2350-2400_20170109.063811.nc"
cdf = netCDF4.Dataset(ref_file)
ref_sst = read_var(cdf,"sst_reynolds")
sza = read_var(cdf,"satellite_zenith_angle")
mask = np.abs(sza) <= 67

ref_sst[~mask] = np.nan

filename = "../sst_approximation/data/2017-01-08/20170108120000-STAR-L2P_GHRSST-SSTsubskin-AHI_H08-ACSPO_V2.41-v02.0-fv01.0.nc"
cdf = netCDF4.Dataset(filename)
sst = read_var(cdf,"sea_surface_temperature")
l2p_flags = read_var(cdf,'l2p_flags')
ql = read_var(cdf,"quality_level")

ql_mask = ql == 5

cloud_mask = np.bitwise_and(l2p_flags,-16384).astype(bool)

sst[~ql_mask] = np.nan

sst[~mask] = np.nan
plt.figure()
plt.imshow(sst)
val = np.nanstd(ref_sst - sst)
print val