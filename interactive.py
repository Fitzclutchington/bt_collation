import netCDF4
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import glob
import scipy.io as sio


collated_files = sorted(glob.glob("data/collated_mat2/*.nc"))

def onclick(event):
    print 'button={}, x={}, y={}, xdata={}, ydata={}'.format(event.button, event.x, event.y, event.xdata, event.ydata)

def read_var(cdf,variable):
    data = np.squeeze(cdf[variable][:])
    return data


original_file = '../sst_approximation/data/2017-01-08/20170108120000-STAR-L2P_GHRSST-SSTsubskin-AHI_H08-ACSPO_V2.41-v02.0-fv01.0.nc'

cdf = netCDF4.Dataset(original_file)

l2p_flags = read_var(cdf,'l2p_flags')
land_mask = np.bitwise_and(l2p_flags,2).astype(bool)
land = np.zeros((5500,5500,4))
r = 146/256.0
g = 98/256.0
b = 57/256.0
land[land_mask] = [r,g,b,1]

filename = collated_files[146]
cdf = netCDF4.Dataset(filename) 
collated = read_var(cdf,"sea_surface_temperature")

fig = plt.figure()

ax1 = plt.subplot(111)
img1 = ax1.imshow(collated,vmin=270,vmax=307)
ax1.imshow(land,interpolation='nearest')
div1 = make_axes_locatable(ax1)
cax1 = div1.append_axes("right", size="5%", pad=0.05)
cbar1 = plt.colorbar(img1, cax=cax1)

plt.show()

cid = fig.canvas.mpl_connect('button_press_event', onclick)