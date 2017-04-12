import netCDF4
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import glob
import scipy.io as sio
import cmocean

original_files = [line.rstrip('\n') for line in open('ahitest.txt')]
#original_files = sorted(glob.glob("../data/2016-11-15/*.nc"))
clear_mask_files = sorted(glob.glob("data/clear/*.nc"))
smooth_files = sorted(glob.glob("data/smooth/*.nc"))
approx_files = sorted(glob.glob("data/approx/*.nc"))
collated_files = sorted(glob.glob("data/collated_mat/*.nc"))
interp_files = sorted(glob.glob("data/collated_interp/*.nc"))
interp_files = [x for x in interp_files if x[31:33] == '00']
pass2_files = sorted(glob.glob("data/pass2/*.nc"))
#pass2_files_old = sorted(glob.glob("data/pass2_old/*.nc"))

filenames = [x.split('/')[2] for x in smooth_files]
original_filenames = [x.split('/')[6] for x in original_files]

TL0_files = []
TL1_files = []
TL2_files = []
for i in range(len(pass2_files)):
    TL0_files.append("data/TL0/TL0_"+str(i)+".nc")
    TL1_files.append("data/TL1/TL1_"+str(i)+".nc")
    TL2_files.append("data/TL2/TL2_"+str(i)+".nc")

filter_window_lag = 1
second_pass_lag = 25/2
smooth_lag = 19/2
collated_inds = []

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

cloud_mask = np.bitwise_and(l2p_flags,-16384).astype(bool)
cloud = np.zeros((5500,5500,4))
r = 173/256.0
g = 216/256.0
b = 230/256.0
cloud[cloud_mask] = [r,g,b,1]

filename = 'data/approx/20170108120000-STAR-L2P_GHRSST-SSTsubskin-AHI_H08-ACSPO_V2.41-v02.0-fv01.0.nc'
cdf = netCDF4.Dataset(filename) 
approx = read_var(cdf,"sea_surface_temperature")

filename = 'data/collated_mat/20170108120000-STAR-L2P_GHRSST-SSTsubskin-AHI_H08-ACSPO_V2.41-v02.0-fv01.0.nc'
cdf = netCDF4.Dataset(filename) 
collated = read_var(cdf,"sea_surface_temperature")

filename = 'data/smooth/20170108120000-STAR-L2P_GHRSST-SSTsubskin-AHI_H08-ACSPO_V2.41-v02.0-fv01.0.nc'
cdf = netCDF4.Dataset(filename) 
smooth = read_var(cdf,"sea_surface_temperature")

filename = '../sst_approximation/data/2017-01-08/20170108120000-STAR-L2P_GHRSST-SSTsubskin-AHI_H08-ACSPO_V2.41-v02.0-fv01.0.nc'
cdf = netCDF4.Dataset(filename) 
sst = read_var(cdf,"sea_surface_temperature")


sst_acspo = sst.copy()
sst_acspo[cloud_mask] = np.nan

filename = 'data/clear/20170108120000-STAR-L2P_GHRSST-SSTsubskin-AHI_H08-ACSPO_V2.41-v02.0-fv01.0.nc'
cdf = netCDF4.Dataset(filename) 
clear = read_var(cdf,"brightness_temperature_11um2")
clear_mask = clear == 255

sst_clear = sst.copy()
sst_clear[~clear_mask] = np.nan 


filename = 'data/pass2/20170108120000-STAR-L2P_GHRSST-SSTsubskin-AHI_H08-ACSPO_V2.41-v02.0-fv01.0.nc'
cdf = netCDF4.Dataset(filename) 
clear = read_var(cdf,"brightness_temperature_11um2")
clear_mask = clear == 255

sst_pass2 = sst.copy()
sst_pass2[~clear_mask] = np.nan 

ref_file = "/home/fitz/sst_approximation/data/ACSPO_V2.41_H08_AHI_2017-01-09_0050-0100_20170110.033449.nc"
cdf = netCDF4.Dataset(ref_file)
ref_sst = read_var(cdf,"sst_reynolds")
sza = read_var(cdf,"satellite_zenith_angle")
mask = np.abs(sza) < 67

ref_sst[~mask] = np.nan

plt.figure()
#cmap = cmocean.cm.thermal

ax1 = plt.subplot(221)
img1 = ax1.imshow(collated,vmin=270,vmax=307)
#ax1.imshow(cloud,interpolation='nearest')
ax1.imshow(land,interpolation='nearest')
div1 = make_axes_locatable(ax1)
cax1 = div1.append_axes("right", size="5%", pad=0.05)
cbar1 = plt.colorbar(img1, cax=cax1)


ax2 = plt.subplot(222, sharex=ax1,sharey=ax1)
img2 = ax2.imshow(sst_pass2,vmin=270,vmax=307)
#ax2.imshow(cloud,interpolation='nearest')
ax2.imshow(land,interpolation='nearest')
div2 = make_axes_locatable(ax2)
cax2 = div2.append_axes("right", size="5%", pad=0.05)
cbar2 = plt.colorbar(img2, cax=cax2)

ax3 = plt.subplot(223, sharex=ax1,sharey=ax1)
img3 = ax3.imshow(sst_acspo,vmin=270,vmax=307)
#ax2.imshow(cloud,interpolation='nearest')
ax3.imshow(land,interpolation='nearest')
div3 = make_axes_locatable(ax3)
cax3 = div3.append_axes("right", size="5%", pad=0.05)
cbar3 = plt.colorbar(img3, cax=cax3)

ax4 = plt.subplot(224, sharex=ax1,sharey=ax1)
img4 = ax4.imshow(approx,vmin=270,vmax=307)
#ax2.imshow(cloud,interpolation='nearest')
ax4.imshow(land,interpolation='nearest')
div4 = make_axes_locatable(ax4)
cax4 = div4.append_axes("right", size="5%", pad=0.05)
cbar4 = plt.colorbar(img4, cax=cax4)

mask = sst >= collated - 0.2
new_pass2 = sst_pass2.copy()
new_pass2[mask] = sst[mask]

plt.figure()
#cmap = cmocean.cm.thermal

ax1 = plt.subplot(121)
img1 = ax1.imshow(sst_pass2,vmin=270,vmax=307)
#ax1.imshow(cloud,interpolation='nearest')
ax1.imshow(land,interpolation='nearest')
div1 = make_axes_locatable(ax1)
cax1 = div1.append_axes("right", size="5%", pad=0.05)
cbar1 = plt.colorbar(img1, cax=cax1)

ax2 = plt.subplot(122, sharex=ax1,sharey=ax1)
img2 = ax2.imshow(new_pass2,vmin=270,vmax=307)
#ax2.imshow(cloud,interpolation='nearest')
ax2.imshow(land,interpolation='nearest')
div2 = make_axes_locatable(ax2)
cax2 = div2.append_axes("right", size="5%", pad=0.05)
cbar2 = plt.colorbar(img2, cax=cax2)
"""
plt.figure()
#cmap = cmocean.cm.thermal

ax1 = plt.subplot(221)
img1 = ax1.imshow(sst,vmin=270,vmax=307)
#ax1.imshow(cloud,interpolation='nearest')
ax1.imshow(land,interpolation='nearest')
div1 = make_axes_locatable(ax1)
cax1 = div1.append_axes("right", size="5%", pad=0.05)
cbar1 = plt.colorbar(img1, cax=cax1)


ax2 = plt.subplot(222, sharex=ax1,sharey=ax1)
img2 = ax2.imshow(smooth,vmin=270,vmax=307)
#ax2.imshow(cloud,interpolation='nearest')
ax2.imshow(land,interpolation='nearest')
div2 = make_axes_locatable(ax2)
cax2 = div2.append_axes("right", size="5%", pad=0.05)
cbar2 = plt.colorbar(img2, cax=cax2)

ax3 = plt.subplot(223, sharex=ax1,sharey=ax1)
img3 = ax3.imshow(approx,vmin=270,vmax=307)
#ax2.imshow(cloud,interpolation='nearest')
ax3.imshow(land,interpolation='nearest')
div3 = make_axes_locatable(ax3)
cax3 = div3.append_axes("right", size="5%", pad=0.05)
cbar3 = plt.colorbar(img3, cax=cax3)

ax4 = plt.subplot(224, sharex=ax1,sharey=ax1)
img4 = ax4.imshow(collated,vmin=270,vmax=307)
#ax2.imshow(cloud,interpolation='nearest')
ax4.imshow(land,interpolation='nearest')
div4 = make_axes_locatable(ax4)
cax4 = div4.append_axes("right", size="5%", pad=0.05)
cbar4 = plt.colorbar(img4, cax=cax4)
"""
ref_file = "../data/ACSPO_V2.41_H08_AHI_2016-09-29_0010-0020_20160930.022950.nc"
cdf = netCDF4.Dataset(ref_file)
ref_sst = read_var(cdf,"sst_reynolds")

print "std collated0 - ref = {}".format(np.nanstd(col - ref_sst)) 
print "std collated1 - ref = {}".format(np.nanstd(col_old - ref_sst))
print "std collated2 - ref = {}".format(np.nanstd(col_26 - ref_sst))
#points = [(2705,5186),(2701,5174),(2692,5198),(2664,5209),(2697,5249),(3621,4180),(3625,4202),(3663,4135),(3774,4567),(3800,4766)]
#points = [(2713,5040),(2709,5061),(2691,5056)]
#point = [(1340,2350),(1329,2370),(1241,2561),(1246,2547)]
points = [(3312,1888),(3258,1952),(3229,1819),(3416,1777),(3173,1667),(544,3237),(524,3173),(573,3190)]
points = [(3351, 3398),(3353, 3400),(2639, 4975),(544,3237),(524,3173),(573,3190),(2664,5209),(2692,5198),(2701,5174),(2705,5186)]
#points = [(2568,3836),(3691,1716),(1519,1880),(4486,1815)]
#points = [(2661,3310),(3306,2087),(2833,3344)]
#points = [(360,3498),(301,3510),(3457,5123),(3473,5113),(2345,4602),(2427,4663),(5102,2775),(5116,2796)]
#points = [(2450,4450),(2449,4424),(2426,4440),(1903,1322),(4693,3096),(4683,3051),(1182,2539),(1189,2497)]
#points = [(1331,2959),(1239,4582),(1420,1840),(2659,415),(3737,550),(3765,523),(3512,1621),(4703,1482),(4712,1470)]
#points = [(2397,3564),(4388,1383)]
times = []

for approx_file in approx_files:
    times.append(approx_file[20:24])

approx = {}
sst = {}
second_pass = {}
clear = {}
smooth = {}
collated = {}
interpolated = {}

for i in points:
    approx[i] = []
    sst[i] = []
    clear[i] = []
    smooth[i] = []
    collated[i] = []
    second_pass[i] = []
    interpolated[i] = []

total_files = len(approx_files)
for i in range(total_files):
    filename = approx_files[i]
    approxnc=netCDF4.Dataset(filename)
    approx_val = read_var(approxnc,'sea_surface_temperature')
   
    filename = smooth_files[i]
    smoothnc=netCDF4.Dataset(filename)
    smooth_val = read_var(smoothnc,'sea_surface_temperature')

    filename = pass2_files[i+smooth_lag]
    pass2nc=netCDF4.Dataset(filename)
    sst_clear_pass2 = read_var(pass2nc,'brightness_temperature_11um2')

    filename = clear_mask_files[i+second_pass_lag+smooth_lag]
    clearnc=netCDF4.Dataset(filename)
    sst_clear_val = read_var(clearnc,'brightness_temperature_11um2')

    filename = original_files[i+filter_window_lag+smooth_lag+second_pass_lag]
    orignc=netCDF4.Dataset(filename)
    sst_val = read_var(orignc,'sea_surface_temperature')
 


    for j in points:
        
        sst[j].append(sst_val[j[0],j[1]])
        
        approx[j].append(approx_val[j[0],j[1]])
          
        smooth[j].append(smooth_val[j[0],j[1]])
        
        if(sst_clear_val[j[0],j[1]] == 255):
        	clear[j].append(sst[j][i])
        else:
        	clear[j].append(np.nan)

        if(sst_clear_pass2[j[0],j[1]]):
            second_pass[j].append(sst[j][i])
        else:
        	second_pass[j].append(np.nan)

    


for i in range(len(collated_files)):
    filename = collated_files[i]
    gen=netCDF4.Dataset(filename)
    original = read_var(gen,'sea_surface_temperature')
    for j in points:
        collated[j].append(original[j[0],j[1]])

    #filename = interp_files[i]
    
    #gen=netCDF4.Dataset(filename)
    #original = read_var(gen,'sea_surface_temperature')
    #for j in points:
    #    interpolated[j].append(original[j[0],j[1]])

collated_inds = []

for f in collated_files:
    collated_inds.append(filenames.index(f.split('/')[2])) 


for i in points:
    outfile = str(i)
    
    #sio.savemat(outfile,{"sst":sst[i],'approx':approx[i],"second_pass":second_pass[i],"clear":clear[i],"smooth":smooth[i],"interpolated":interpolated[i],"collated":collated[i], "times":times})
    
    fig = plt.figure(figsize=(12,10))
    fig.canvas.set_window_title(str(i))
    plt.title(str(i) + "\n")
    plt.plot(range(total_files),approx[i],'.',label="approx")
    
    plt.plot(range(total_files),second_pass[i],'o', markersize=15,label="pass2")
    plt.plot(range(total_files),clear[i],'.', markersize=15,label="clear")
    plt.plot(range(total_files),sst[i],label="orig")
    plt.plot(range(total_files),smooth[i],label="smooth")
    
    #plt.plot(collated_inds,interpolated[i],'.',markersize=15,label="interp")
    plt.plot(collated_inds,collated[i],'.',markersize=15,label="collated")
    plt.legend()
    plt.xticks(range(total_files),times,rotation='vertical')
    #plt.savefig("Time_at_" + str(i[0]) + "_" + str(i[1]) + ".png")
    #plt.show()
    
########################################################################################################################################
points = [(544,3237),(524,3173),(573,3190)]
points = [(2701,5174),(2705,5186)]
points = [(3351, 3398),(3353, 3400),(2639, 4975),(544,3237),(524,3173),(573,3190),(2664,5209),(2692,5198),(2701,5174),(2705,5186)]
#TL0 = {}
#TL1 = {}
#TL2 = {}
sst = {}
second_pass = {}
clear = {}
times = []

for approx_file in approx_files:
    times.append(approx_file[20:24])


for i in points:
    sst[i] = []
    clear[i] = []
    second_pass[i] = []
    #TL0[i] = []
    #TL1[i] = []
    #TL2[i] = []
    
total_files = len(approx_files)
for i in range(total_files):


    filename = pass2_files[i]
    pass2nc=netCDF4.Dataset(filename)
    sst_clear_pass2 = read_var(pass2nc,'brightness_temperature_11um2')

    filename = clear_mask_files[i+second_pass_lag]
    clearnc=netCDF4.Dataset(filename)
    sst_clear_val = read_var(clearnc,'brightness_temperature_11um2')

    filename = original_files[i+filter_window_lag+second_pass_lag]
    orignc=netCDF4.Dataset(filename)
    sst_val = read_var(orignc,'sea_surface_temperature')
    """
    filename = TL0_files[i]
    orignc = netCDF4.Dataset(filename)
    TL0_val = read_var(orignc, "brightness_temperature_11um2")

    filename = TL1_files[i]
    orignc = netCDF4.Dataset(filename)
    TL1_val = read_var(orignc, "brightness_temperature_11um2")

    filename = TL2_files[i]
    orignc = netCDF4.Dataset(filename)
    TL2_val = read_var(orignc, "brightness_temperature_11um2")
    """
    for j in points:
        
        sst[j].append(sst_val[j[0],j[1]])
        
        #TL0[j].append(TL0_val[j[0],j[1]])
        #TL1[j].append(TL1_val[j[0],j[1]])
        #TL2[j].append(TL2_val[j[0],j[1]])
        
        if(sst_clear_val[j[0],j[1]] == 255):
            clear[j].append(sst[j][i])
        else:
            clear[j].append(np.nan)

        if(sst_clear_pass2[j[0],j[1]]):
            second_pass[j].append(sst[j][i])
        else:
            second_pass[j].append(np.nan)

        


for i in points:
    fig = plt.figure(figsize=(12,10))
    fig.canvas.set_window_title(str(i))
    plt.title(str(i) + "\n")
    
    plt.plot(range(total_files),second_pass[i],'o', markersize=15,label="pass2")
    plt.plot(range(total_files),clear[i],'.', markersize=15,label="clear")
    plt.plot(range(total_files),sst[i],'k',label="orig")
    #plt.plot(range(total_files),TL0[i],'r',label="TL0")
    #plt.plot(range(total_files),TL1[i],'g',label="TL1")
    #plt.plot(range(total_files),TL2[i],'b', label="TL2")
    plt.xticks(range(total_files),times,rotation='vertical')

    plt.legend()

######################################################################################################################
clear_count = []
pass2_count = []


filenames = [x.split('/')[2] for x in pass2_files]
clear_mask_files = sorted(glob.glob("data/clear/*.nc"))
pass2_files = sorted(glob.glob("data/pass2_old/*.nc"))
for clear_file, pass2_file in zip(clear_mask_files[12:12+len(pass2_files)],pass2_files):
    clear_cdf = netCDF4.Dataset(clear_file)
    clear = read_var(clear_cdf,"brightness_temperature_11um2")
    
    pass2_cdf = netCDF4.Dataset(pass2_file)
    pass2 = read_var(pass2_cdf,"brightness_temperature_11um2")
    
    
    #orig_cdf = netCDF4.Dataset(original_file)
    #l2p_flags = read_var(orig_cdf,"l2p_flags")

    #acspo_count.append(np.bitwise_and(l2p_flags,-16384).astype(bool).astype(int).sum())
    clear_count.append((clear == 255).sum())
    pass2_count.append((pass2 == 255).sum())




plt.figure()
plt.plot(clear_count,label="first pass")
plt.plot(pass2_count,label="second pass")

#plt.plot(acspo_count,label="acspo")
plt.legend()

points = [(3370,5263),(3426,5215),(2505,295),(3984,3312),(1458,2631),(3540,1824),(3511,1767),(4460,3920),(3682,1725)]

sst = {}
clear = {}
pass2 = {}


for i in points:
    sst[i] = []
    clear[i] = []
    pass2[i] = []
    
clear_mask_files = sorted(glob.glob("data/clear/*.nc"))[12:128]

for i in range(total_files):

    filename = clear_mask_files[i]
    clearnc=netCDF4.Dataset(filename)
    sst_clear_val = read_var(clearnc,'brightness_temperature_11um2')

    filename = pass2_files[i]
    pass2nc=netCDF4.Dataset(filename)
    pass2_val = read_var(pass2nc,'brightness_temperature_11um2')

    filename = original_files[i+1]
    orignc=netCDF4.Dataset(filename)
    sst_val = read_var(orignc,'sea_surface_temperature')


    for j in points:
        
        sst[j].append(sst_val[j[0],j[1]])
        
                
        if(sst_clear_val[j[0],j[1]] == 255):
            clear[j].append(sst[j][i])
        else:
            clear[j].append(np.nan)

        if(pass2_val[j[0],j[1]]):
            pass2[j].append(sst[j][i])
        else:
            pass2[j].append(np.nan)


for i in points:
    plt.figure(figsize=(12,10))
   
    plt.plot(pass2[i],'.', markersize=15,label="second pass")
    plt.plot(clear[i],'r.', markersize=10,label="first pass")
    
    plt.plot(sst[i],label="sst")
    plt.legend()
    plt.title("Time series at point {}".format(str(i)))
    plt.savefig("Time_at_" + str(i[0]) + "_" + str(i[1]) + ".png")
    #plt.show()

lag0_files = []
for i in range(0,126):
    lag0_files.append("data/lag0/"+str(i)+".nc")

totals =[]
for lag_file in lag0_files:
    lag_cdf = netCDF4.Dataset(lag_file)
    lag = read_var(lag_cdf,"brightness_temperature_11um2")

    total = np.nansum(lag)
    totals.append(total)


plt.figure()
plt.plot(totals)


reference_file = "../data/ACSPO_V2.41_H08_AHI_2016-11-15_0300-0310_20161116.034934.nc"
ref_cdf = netCDF4.Dataset(reference_file)
sst_ref = read_var(ref_cdf,"sst_reynolds")

collated_cdf = netCDF4.Dataset(collated_files[0])
sst_collated = read_var(collated_cdf,"sea_surface_temperature")

sst_cdf = netCDF4.Dataset(original_files[filenames.index(collated_files[0].split("/")[2]) + filter_window_lag + smooth_lag + second_pass_lag])
sst = read_var(sst_cdf,"sea_surface_temperature")

l2p_flags = read_var(sst_cdf,'l2p_flags')
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

"""
plt.figure()


ax1 = plt.subplot(131)
img1 = ax1.imshow(sst,vmin=270,vmax=307)
ax1.imshow(cloud,interpolation='nearest')
ax1.imshow(land,interpolation='nearest')
div1 = make_axes_locatable(ax1)
cax1 = div1.append_axes("right", size="5%", pad=0.05)
cbar1 = plt.colorbar(img1, cax=cax1)


ax2 = plt.subplot(132, sharex=ax1,sharey=ax1)
img2 = ax2.imshow(sst_collated,vmin=270,vmax=307)
#ax2.imshow(cloud,interpolation='nearest')
ax2.imshow(land,interpolation='nearest')
div2 = make_axes_locatable(ax2)
cax2 = div2.append_axes("right", size="5%", pad=0.05)
cbar2 = plt.colorbar(img2, cax=cax2)

ax3 = plt.subplot(133, sharex=ax1, sharey=ax1)
img3 = ax3.imshow(sst_collated - sst_ref,vmin=-2,vmax=2)
ax3.imshow(land,interpolation='nearest')
div3 = make_axes_locatable(ax3)
cax3 = div3.append_axes("right", size="5%", pad=0.05)
cbar3 = plt.colorbar(img3, cax=cax3)
"""



points = [(3370,5263),(3426,5215),(2505,295),(3984,3312),(1458,2631),(3540,1824),(3511,1767),(4460,3920),(3682,1725),(726,2379),(720,2427)]
points = [(3984,3312),(4460,3920),(3682,1725),(726,2379),(720,2427)]
points = [(1449,2721),(1419,2758),(3856,3276),(3800,3161),(1772,1121)]
points = [(3984,3312),(720,2427),(726,2379)]
#approx = {}
points = [(3370,5263),(3426,5215),(2505,295),(3984,3312),(1458,2631),(3540,1824),(3511,1767),(4460,3920),(3682,1725),(726,2379),(720,2427),(4180,454),(4048,372),(2171,178),(2178,169)]
points = [(1080,2160),(1090,2160),(1100,2160),(1110,2160),(1085,2160),(1095,2160),(1105,2160)]

points = [(1110,2160),(1095,2160),(1105,2160)]
approx = {}
sst = {}
second_pass = {}
clear = {}
smooth = {}
collated = {}
second_pass_old = {}


for i in points:
    approx[i] = []
    sst[i] = []
    clear[i] = []
    smooth[i] = []
    collated[i] = []
    second_pass[i] = []
    second_pass_old[i] = []

total_files = len(approx_files)
for i in range(total_files):
    
    filename = approx_files[i]
    approxnc=netCDF4.Dataset(filename)
    approx_val = read_var(approxnc,'brightness_temperature_11um2')
    
    filename = smooth_files[i]
    smoothnc=netCDF4.Dataset(filename)
    smooth_val = read_var(smoothnc,'brightness_temperature_11um2')
    
    filename = pass2_files[i+smooth_lag]
    pass2nc=netCDF4.Dataset(filename)
    sst_clear_pass2 = read_var(pass2nc,'brightness_temperature_11um2')
    
    filename = pass2_files_old[i+smooth_lag]
    pass2nc_old=netCDF4.Dataset(filename)
    sst_clear_pass2_old = read_var(pass2nc_old,'brightness_temperature_11um2')
    
    """
    filename = TL1_files[i]
    TL1nc_new = netCDF4.Dataset(filename)
    TL1_new = read_var(TL1nc_new,'brightness_temperature_11um2')

    filename = TL2_files[i]
    TL2nc_new = netCDF4.Dataset(filename)
    TL2_new = read_var(TL2nc_new,'brightness_temperature_11um2')
    """

    filename = clear_mask_files[i+second_pass_lag+smooth_lag]
    clearnc=netCDF4.Dataset(filename)
    sst_clear_val = read_var(clearnc,'brightness_temperature_11um2')

    filename = original_files[i+filter_window_lag+second_pass_lag+smooth_lag]
    orignc=netCDF4.Dataset(filename)
    sst_val = read_var(orignc,'brightness_temperature_11um2')
 

    for j in points:
        
        sst[j].append(sst_val[j[0],j[1]])
        
        approx[j].append(approx_val[j[0],j[1]])
          
        smooth[j].append(smooth_val[j[0],j[1]])


        if(sst_clear_val[j[0],j[1]] == 255):
            clear[j].append(sst[j][i])
        else:
            clear[j].append(np.nan)

        if(sst_clear_pass2[j[0],j[1]]):
            second_pass[j].append(sst[j][i])
        else:
            second_pass[j].append(np.nan)

        if(sst_clear_pass2_old[j[0],j[1]]):
            second_pass_old[j].append(sst[j][i])
        else:
            second_pass_old[j].append(np.nan)

        


for i in range(len(collated_files)):
    filename = collated_files[i]
    gen=netCDF4.Dataset(filename)
    original = read_var(gen,'brightness_temperature_11um2')
    for j in points:
        collated[j].append(original[j[0],j[1]])

collated_inds = []

for f in collated_files:
    collated_inds.append(filenames.index(f.split('/')[2]))



for i in points:
    plt.figure(figsize=(12,10))
    plt.title(str(i) + "\n")
    #plt.plot(range(total_files),approx[i],label="approx")
    
    
    
    plt.plot(range(total_files),clear[i],'.', markersize=10,label="clear")
    
    plt.plot(range(total_files),second_pass[i],'o', markersize=15,label="pass2")
    plt.plot(range(total_files),second_pass_old[i],'o', markersize=10,label="pass2_old")
    plt.plot(range(total_files),sst[i],label="orig")
    plt.plot(range(total_files),smooth[i],label="smooth")
    plt.plot(range(total_files),approx[i],label="approx")
    plt.plot(collated_inds,collated[i],'.',markersize=15,label="collated")
    plt.legend()

clear_count = np.zeros((5500,5500))
for clear_file in clear_mask_files:
    clearnc=netCDF4.Dataset(clear_file)
    sst_clear_val = read_var(clearnc,'brightness_temperature_11um2')
    mask = sst_clear_val == 255
    clear_count[mask] += 1

plt.figure()
plt.imshow(clear_count,vmin=0,vmax=len(clear_mask_files))
plt.colorbar()

original_file = original_files[13+80]
cdf = netCDF4.Dataset(original_file)
sst = read_var(cdf,"brightness_temperature_11um2")
plt.figure()
plt.imshow(sst,vmin=263,vmax=301)
plt.colorbar()

diffs = []
for i in range(len(smooth_files)-1):
    smoothcdf = netCDF4.Dataset(smooth_files[i])
    smooth1 = read_var(smoothcdf,"brightness_temperature_11um2")

    smoothcdf = netCDF4.Dataset(smooth_files[i+1])
    smooth2 = read_var(smoothcdf,"brightness_temperature_11um2")
    diffs.append(smooth1[726,2379] - smooth2[726,2379])

plt.figure()
plt.plot(diffs)
#plt.plot(smooth[(726,2379)])