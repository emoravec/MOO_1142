"""
Makes plots of minkasi maps.

Created in 2024-07
"""
import aplpy

from astropy.nddata import Cutout2D
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from astropy.nddata import Cutout2D
from astropy import units as u

import matplotlib
import matplotlib.pyplot as pyplot
import numpy as np
import os
import time
# -------------------------------------------------------------------------------------------- #
# colormap parameters blue-orange colormap
import matplotlib.cm as cm

bottom = cm.get_cmap('Oranges', 128)
top = cm.get_cmap('Blues_r', 128)
newcolors = np.vstack((top(np.linspace(0, 1, 128)),bottom(np.linspace(0, 1, 128))))

from matplotlib.colors import ListedColormap
cm.register_cmap('OrangeBlue', cmap = ListedColormap(newcolors))
cmap = 'OrangeBlue'

moo1142_cen = SkyCoord(175.6909647,15.4551255,unit='deg')
# ------------------ #
def do_plotting(filename,fits_path,save_path):
    print(fits_path)
    minkasi_image_hdu = fits.open(fits_path)
    minkasi_image_data = minkasi_image_hdu[0].data
    minkasi_image_wcs = WCS(minkasi_image_hdu[0].header)
    # make cutout for scaling purposes
    img_size = u.Quantity((5.0,5.0), u.arcmin)
    cutout = Cutout2D(minkasi_image_data, position=moo1142_cen, size=img_size, wcs=minkasi_image_wcs)
    cutout_max = np.max(cutout.data)
    print('The max for the cutout is:',cutout_max)
    cutout_min = np.min(cutout.data)
    print('The min for the cutout is:',cutout_min)
    # ------------------ #
    # basic figure properties
    fig = pyplot.figure(figsize=(5, 5))
    img = aplpy.FITSFigure(fits_path, figure=fig, downsample=1, smooth=False, convention='calabretta')
    img.set_theme('publication')
    img.set_title(filename)

    # show image
    img.show_colorscale(cmap=cmap,stretch='linear',smooth=1,vmin=cutout_min,vmax=cutout_max) # pixels in the minkasi maps are 2-3"/pixels so a smoothing of smooth=2-4 times that roughly gives the M2 beam of 9" - 3 seems like a lot and 1 seems like not enough?

    # recenter and add center marker
    img.recenter(moo1142_cen.ra.degree, moo1142_cen.dec.degree, radius=2/60.0) #radius=2.0/60.0
    img.show_markers(moo1142_cen.ra.degree,moo1142_cen.dec.degree,facecolor='black',edgecolor=None,marker='+',s=50,linewidths=2,alpha = 0.5)
    
    # other plotting params
    img.ax.tick_params(axis='both',which='both',direction='in')
    img.add_scalebar(0.5/60.0, '30\"', color='black')

    # add M2 beam
    img.add_beam(major=9.0/3600.0, minor=9.0/3600.0, angle=0)
    img.beam.set_color('white')
    img.beam.set_edgecolor('green')
    img.beam.set_facecolor('white')
    img.beam.set_corner('bottom left')

    # add colorbar
    img.add_colorbar('right')
    img.colorbar.set_width(0.12)
    img.colorbar.set_axis_label_text('K$_{rj}$') # minkasi maps are usually brightness temperature, K_rj (Tony M July 3rd, 2024 Slack)

    img.save(save_path + 'pdfs/' + filename +'.pdf', format='pdf', dpi=300)

def plot_minkasi_maps(N,M,tod_class,save_path):
    # make initial image
    initial_filename = 'MOO_1142_'+tod_class+'_initial'
    image_fits_path = tod_class_path + initial_filename +'.fits'
    do_plotting(initial_filename,image_fits_path,save_path)

    # plot iterations
    for n in N:
        for m in M:
            # read in file and image properties
            n_m_file_name = 'MOO_1142_'+tod_class+'_niter_'+str(n)+'_'+str(m)
            image_fits_path = tod_class_path + n_m_file_name +'.fits'
            do_plotting(n_m_file_name,image_fits_path,save_path)

# ------------------ #
location = '/Users/emoravec/Documents/Research/merging_clusters/analysis/MOO_1142/M2/images/minkasi/'
tod_class = 'good'
tod_class_path = location + tod_class + '/'
N_array = [1,2,3,4,5]
M_array = [1,5,15,25,50]
# ------------------ #
# Record the start time
start_time = time.time()

# do plotting
plot_minkasi_maps(N_array,M_array,tod_class,tod_class_path)

# Record the end time
end_time = time.time()

# Calculate the elapsed time
elapsed_time = end_time - start_time

# Convert the elapsed time to hours and minutes
hours = int(elapsed_time // 3600)
minutes = int((elapsed_time % 3600) // 60)
seconds = int(elapsed_time % 60)

# Print the elapsed time in hours and minutes
print(f"Elapsed time: {hours} hours, {minutes} minutes, and {seconds} seconds")
