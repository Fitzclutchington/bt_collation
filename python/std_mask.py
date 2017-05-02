

"""
File to test the brightness temperature masks for area with cloud leakage
"""
from scipy.ndimage.filters import uniform_filter
import netCDF4
import numpy as np
import glob
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

def read_var(cdf,variable):
    data = np.squeeze(cdf[variable][:])
    return data

cdf = netCDF4.Dataset("../../sst_approximation/data/2017-04-01/20170401191000-STAR-L2P_GHRSST-SSTsubskin-AHI_H08-ACSPO_V2.41-v02.0-fv01.0.nc")

bt08 = read_var(cdf, "brightness_temperature_08um6")
bt10 = read_var(cdf, "brightness_temperature_10um4")
bt11 = read_var(cdf, "brightness_temperature_11um2")
bt12 = read_var(cdf, "brightness_temperature_12um3")
sst_test =  read_var(cdf, "sea_surface_temperature")
sza = read_var(cdf,"satellite_zenith_angle")



l2p_flags = read_var(cdf,'l2p_flags')
land_mask = np.bitwise_and(l2p_flags,2).astype(bool)
land = np.zeros((5500,5500,4))
r = 146/256.0
g = 98/256.0
b = 57/256.0
land[land_mask] = [r,g,b,1]


c1 = uniform_filter(sst_test, (3,3))
c2 = uniform_filter(sst_test*sst_test, (3,3))
std = ((c2 - c1*c1)**.5)

cdf = netCDF4.Dataset("../res_mean.nc")
eigen = read_var(cdf,"data");

m1 = std > .6
m2 = eigen < .06

mask = np.logical_and(m1,m2)


plt.figure()
#cmap = cmocean.cm.thermal


ax1 = plt.subplot(221)
img1 = ax1.imshow(std,vmin=0.4,vmax=1)
#ax1.imshow(cloud,interpolation='nearest')
ax1.imshow(land,interpolation='nearest')
div1 = make_axes_locatable(ax1)
cax1 = div1.append_axes("right", size="5%", pad=0.05)
cbar1 = plt.colorbar(img1, cax=cax1)

ax2 = plt.subplot(222, sharex=ax1,sharey=ax1)
img2 = ax2.imshow(eigen,vmin=0.03,vmax=0.1)
#ax2.imshow(cloud,interpolation='nearest')
ax2.imshow(land,interpolation='nearest')
div2 = make_axes_locatable(ax2)
cax2 = div2.append_axes("right", size="5%", pad=0.05)
cbar2 = plt.colorbar(img2, cax=cax2)

ax2 = plt.subplot(223, sharex=ax1,sharey=ax1)
img2 = ax2.imshow(sst_test,vmin=270,vmax=307)
#ax2.imshow(cloud,interpolation='nearest')
ax2.imshow(land,interpolation='nearest')
div2 = make_axes_locatable(ax2)
cax2 = div2.append_axes("right", size="5%", pad=0.05)
cbar2 = plt.colorbar(img2, cax=cax2)

ax2 = plt.subplot(224, sharex=ax1,sharey=ax1)
img2 = ax2.imshow(mask)
#ax2.imshow(cloud,interpolation='nearest')
ax2.imshow(land,interpolation='nearest')
div2 = make_axes_locatable(ax2)
cax2 = div2.append_axes("right", size="5%", pad=0.05)
cbar2 = plt.colorbar(img2, cax=cax2)


plt.show()

plt.figure()
#cmap = cmocean.cm.thermal

