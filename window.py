#make movie to show shifrt in time

import netCDF4
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import animation
import glob

sst_files = [line.rstrip('\n') for line in open('ahinov.txt')]
clear_mask_files = sorted(glob.glob("data/clear/*.nc"))

fig = plt.figure()
sst_ax = plt.subplot(121)
clear_ax = plt.subplot(122)
sst_im, = sst_ax.imshow([],vmin=292,vmax=307, animated=True, cmap='cool')
clear_im = clear_ax.imshow([],vmin=0,vmax=255, animated=True)

# initialization function: plot the background of each frame
def init():
    sst_im.set_array([])
    clear_im.set_array([])
    return sst_im, , clear_im

# animation function.  This is called sequentially
def animate(i):

    line.set_data(x, y)
    return line,

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=200, interval=20, blit=True)

# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html
anim.save('basic_animation.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

plt.show()