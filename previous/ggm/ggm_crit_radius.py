"""
2024-07
Script to plot a ggm image and plot various features of interest for comparison.
"""
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy import units as u
import aplpy
import matplotlib.pyplot as plt
import numpy as np
# ------------------------------------------------------------
plt.rcParams.update({'font.size': 14})

cen = SkyCoord(175.6909647,15.4551255,unit='deg')
location = '/Users/emoravec/Documents/Research/merging_clusters/analysis/MOO_1142/'
# ------------------------------------------------------------
def plot_ggm(input_image):
    hdu = fits.open(input_image)
    fig = plt.figure(figsize=(8,8))
    # plot the new cutout
    img = aplpy.FITSFigure(data=hdu,figure=fig,convention='wells')
    img.show_colorscale(cmap='RdBu_r')
    img.recenter(x=cen.ra,y=cen.dec,radius=2.50/60)
    # circles
    img.show_circles(cen.ra,cen.dec,radius=100./3600.,color='y',linewidth=2)
    img.show_circles(cen.ra,cen.dec,radius=83./3600.,color='mediumpurple',linewidth=2)
    img.show_markers(cen.ra,cen.dec,s=100,facecolor='k',marker='+')
    # plot properties
    img.ticks.set_tick_direction('in')
    img.ticks.set_color('black')
    # save close image
    plt.savefig(location+'ggm/M2/2024-07/pdfs/ggm_rs/MOO_1142_ggm4_rs.pdf',dpi=300)
    plt.close()
# ------------------------------------------------------------
# M2 images
mooinp = location+'ggm/M2/2024-07/fits/pt_sub/moo1142_M2_only_good.pt_sub.ggm.4.cutout.fits'

plot_ggm(mooinp)

  