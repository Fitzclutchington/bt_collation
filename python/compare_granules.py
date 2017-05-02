import netCDF4
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import glob
import scipy.io as sio
import cmocean
import sys

if len(sys.argv) < 3:
	print "usage: python acspo_mask.py <sst_file> <folder>"
	sys.exit()


def read_var(cdf,variable):
    data = np.squeeze(cdf[variable][:])
    return data


original_file = sys.argv[1]
folder = sys.argv[2]

cdf = netCDF4.Dataset(original_file)

l2p_flags = read_var(cdf,'l2p_flags')
land_mask = np.bitwise_and(l2p_flags,2).astype(bool)
land = np.zeros((5500,5500,4))
r = 146/256.0
g = 98/256.0
b = 57/256.0
land[land_mask] = [r,g,b,1]

cloud_mask = np.bitwise_and(l2p_flags,-16384).astype(bool)
cloud = np.zeros((5500,5500,4))
r = 173/256.0
g = 216/256.0
b = 230/256.0
cloud[cloud_mask] = [r,g,b,1]

sst = read_var(cdf,"sea_surface_temperature")
cdf.close()

compare_file = folder + '/' + original_file.split('/')[-1]
cdf = netCDF4.Dataset(compare_file)
compare_data = read_var(cdf,"sea_surface_temperature")
cdf.close()

reinstated_file = '../data/reinstated' + '/' + original_file.split('/')[-1]
cdf = netCDF4.Dataset(reinstated_file)
reinstated_data = read_var(cdf,"sea_surface_temperature")

acspo = sst.copy()
acspo[cloud_mask] = np.nan

plt.figure()
cmap = cmocean.cm.thermal

ax1 = plt.subplot(221)
img1 = ax1.imshow(compare_data,vmin=270,vmax=307)
ax1.set_title("Collated 2")
ax1.imshow(land,interpolation='nearest')
div1 = make_axes_locatable(ax1)
cax1 = div1.append_axes("right", size="5%", pad=0.05)
cbar1 = plt.colorbar(img1, cax=cax1)


ax2 = plt.subplot(222, sharex=ax1,sharey=ax1)
img2 = ax2.imshow(sst,vmin=270,vmax=307)
ax2.set_title("SST")
#ax2.imshow(cloud,interpolation='nearest')
ax2.imshow(land,interpolation='nearest')
div2 = make_axes_locatable(ax2)
cax2 = div2.append_axes("right", size="5%", pad=0.05)
cbar2 = plt.colorbar(img2, cax=cax2)


ax2 = plt.subplot(223, sharex=ax1,sharey=ax1)
img2 = ax2.imshow(reinstated_data,vmin=270,vmax=307)
ax2.set_title("Reinstated")
#ax2.imshow(cloud,interpolation='nearest')
ax2.imshow(land,interpolation='nearest')
div2 = make_axes_locatable(ax2)
cax2 = div2.append_axes("right", size="5%", pad=0.05)
cbar2 = plt.colorbar(img2, cax=cax2)

ax2 = plt.subplot(224, sharex=ax1,sharey=ax1)
img2 = ax2.imshow(acspo,vmin=270,vmax=307)
ax2.set_title("ACSPO Mask")
#ax2.imshow(cloud,interpolation='nearest')
ax2.imshow(land,interpolation='nearest')
div2 = make_axes_locatable(ax2)
cax2 = div2.append_axes("right", size="5%", pad=0.05)
cbar2 = plt.colorbar(img2, cax=cax2)

plt.show()

