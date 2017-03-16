import netCDF4
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import glob

def read_var(cdf,variable):
    data = np.squeeze(cdf[variable][:])
    return data

collated_files = sorted(glob.glob("data/collated_mat/*.nc"))
original_files = [line.rstrip('\n') for line in open('ahitest.txt')]

cdf = netCDF4.Dataset(original_files[0])
l2p_flags = read_var(cdf,'l2p_flags')
land_mask = np.bitwise_and(l2p_flags,2).astype(bool)
land = np.zeros((5500,5500,4))
r = 146/256.0
g = 98/256.0
b = 57/256.0
land[land_mask] = [r,g,b,1]

cdf = netCDF4.Dataset(collated_files[71])
collated_03 = read_var(cdf,"sea_surface_temperature")

#print "number finite at 0100 = " + str(np.isfinite(collated).sum())

cdf = netCDF4.Dataset(collated_files[77])
collated_04 = read_var(cdf,"sea_surface_temperature")

#print "number finite at 0200 = " + str(np.isfinite(collated).sum())

cdf = netCDF4.Dataset(collated_files[75])
collated_0330 = read_var(cdf,"sea_surface_temperature")

#print "number finite at 0110 = " + str(np.isfinite(collated).sum())

plt.figure()

ax1 = plt.subplot(131)
img1 = ax1.imshow(collated_03,vmin=270,vmax=307)
ax1.imshow(land,interpolation='nearest')
#ax1.imshow(cloud,interpolation='nearest')
div1 = make_axes_locatable(ax1)
cax1 = div1.append_axes("right", size="5%", pad=0.05)
cbar1 = plt.colorbar(img1, cax=cax1)



ax2 = plt.subplot(132,sharex=ax1,sharey=ax1)
img2 = ax2.imshow(collated_04,vmin=270,vmax=307)
ax2.imshow(land,interpolation='nearest')
div2 = make_axes_locatable(ax2)
cax2 = div2.append_axes("right", size="5%", pad=0.05)
cbar2 = plt.colorbar(img2, cax=cax2)

ax3 = plt.subplot(133,sharex=ax1,sharey=ax1)
img3 = ax3.imshow(collated_0330,vmin=270,vmax=307)
ax3.imshow(land,interpolation='nearest')
div3 = make_axes_locatable(ax3)
cax3 = div3.append_axes("right", size="5%", pad=0.05)
cbar3 = plt.colorbar(img3, cax=cax3)

m1 = np.isfinite(collated_04)
m2 = np.isnan(collated_03)

mask1 = np.logical_and(m1,m2)
print mask1.sum()

m1 = np.isfinite(collated_03)
m2 = np.isnan(collated_04)

mask2 = np.logical_and(m1,m2)
print mask2.sum()

plt.figure()
plt.imshow(mask)
plt.show()