"""
File to test the brightness temperature masks for area with cloud leakage
"""

import netCDF4
import numpy as np
import glob
import sys
import matplotlib.pyplot as plt

def read_var(cdf,variable):
    data = np.squeeze(cdf[variable][:])
    return data

cdf = netCDF4.Dataset("../../sst_approximation/data/2017-01-09/20170109160000-STAR-L2P_GHRSST-SSTsubskin-AHI_H08-ACSPO_V2.41-v02.0-fv01.0.nc")

bt08 = read_var(cdf, "brightness_temperature_08um6")
bt10 = read_var(cdf, "brightness_temperature_10um4")
bt11 = read_var(cdf, "brightness_temperature_11um2")
bt12 = read_var(cdf, "brightness_temperature_12um3")

bt_2_mask = (100 * (bt10 - bt08)) / (bt10 + bt08)
bt_3_mask = bt08 + 0.8*bt11 - 1.8*bt12

plt.figure()
plt.imshow(bt_2_mask,vmin=0.2,vmax=1)
plt.colorbar()
plt.show()

plt.figure()
plt.imshow(bt_3_mask,vmin=0.3,vmax=1)
plt.colorbar()
plt.show()

print "2 bt mask = " + str(bt_2_mask[3935,1305])
print "3 bt mask = " + str(bt_3_mask[3935,1305])