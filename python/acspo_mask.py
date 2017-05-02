import netCDF4
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import glob
import scipy.io as sio
import cmocean
import sys

if len(sys.argv) < 2:
	print "usage: python acspo_mask.py <sst_file>"
	sys.exit()


def read_var(cdf,variable):
    data = np.squeeze(cdf[variable][:])
    return data


y_min = 4250
y_max = 4601
x_min = 3900
x_max = 4301

inds = np.index_exp[y_min:y_max, x_min:x_max]

original_file = sys.argv[1]

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
bt08 = read_var(cdf, "brightness_temperature_08um6")
bt10 = read_var(cdf, "brightness_temperature_10um4")
bt11 = read_var(cdf, "brightness_temperature_11um2")
bt12 = read_var(cdf, "brightness_temperature_12um3")

bt_3 = bt08 - 0.5373*bt10 - 0.45*bt12
bt_3_mask = bt_3 < 2

sst_mask = sst.copy()

sst[cloud_mask] = np.nan
sst_mask[bt_3_mask] = np.nan

im = np.zeros((5500,5500))
m = np.logical_and(~cloud_mask,bt_3_mask)

im[m] = 1

plt.figure()
#cmap = cmocean.cm.thermal

plt.imshow(im)
plt.colorbar()
plt.title("acspo clear bt3 cloud")
plt.show()

filename = "../data/pass2/"+original_file.split('/')[-1]
cdf = netCDF4.Dataset(filename)

pass2 = read_var(cdf,"brightness_temperature_11um2")
pass2_mask = pass2 == 255


m = np.logical_and(pass2,bt_3_mask)

im2 = np.zeros((5500,5500))
im2[m] = 1

plt.figure()
#cmap = cmocean.cm.thermal

plt.imshow(im2)
plt.colorbar()
plt.title("pass2 clear bt3 cloud")
plt.show()
"""
plt.figure()
#cmap = cmocean.cm.thermal

ax1 = plt.subplot(121)
img1 = ax1.imshow(sst,vmin=270,vmax=307)
ax1.set_title("ACSPO")
#ax1.imshow(cloud,interpolation='nearest')
ax1.imshow(land,interpolation='nearest')
div1 = make_axes_locatable(ax1)
cax1 = div1.append_axes("right", size="5%", pad=0.05)
cbar1 = plt.colorbar(img1, cax=cax1)


ax2 = plt.subplot(122, sharex=ax1,sharey=ax1)
img2 = ax2.imshow(sst_mask,vmin=270,vmax=307)
ax2.set_title("BT 3 Mask")
#ax2.imshow(cloud,interpolation='nearest')
ax2.imshow(land,interpolation='nearest')
div2 = make_axes_locatable(ax2)
cax2 = div2.append_axes("right", size="5%", pad=0.05)
cbar2 = plt.colorbar(img2, cax=cax2)

plt.show()
"""

"""
y_nocloud = bt08[~cloud_mask] - bt11[~cloud_mask]
x_nocloud = bt11[~cloud_mask] - bt12[~cloud_mask]

y_cloud = bt08[bt_3_mask] - bt11[bt_3_mask]
x_cloud = bt11[bt_3_mask] - bt12[bt_3_mask]

plt.figure()
plt.scatter(x_cloud,y_cloud,label='cloud',c='b',alpha=0.5)
plt.scatter(x_nocloud,y_nocloud,label='no cloud',c='r',alpha=0.5)
plt.ylabel('bt08 - bt11')
plt.xlabel('bt11-bt12')
plt.grid()
plt.legend()

plt.show()
"""