"""
Plots M2 with XMM contours and spec or photzs.

Created in 2024-07
"""
import aplpy

from astropy.nddata import Cutout2D
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.table import Table,vstack
from astropy import units as u
from astropy.wcs import WCS

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
# -------------------------------------------------------------------------------------------- #
# colormap parameters blue-orange colormap
import matplotlib.cm as cm

bottom = cm.get_cmap('Oranges', 128)
top = cm.get_cmap('Blues_r', 128)
newcolors = np.vstack((top(np.linspace(0, 1, 128)),bottom(np.linspace(0, 1, 128))))

from matplotlib.colors import ListedColormap
cm.register_cmap('OrangeBlue', cmap = ListedColormap(newcolors))
cmap = 'OrangeBlue'
# -------------------------------------------------------------------------------------------- #
# parameter definitions
moo1142_cen = SkyCoord(175.6909647,15.4551255,unit='deg')
# -------------------------------------------------------------------------------------------- #
# files
location = Path('/Users/emoravec/Documents/Research/merging_clusters/analysis/MOO_1142/')
M2_image_path = location / 'M2/images/midas/2024-07-03_Kelvin_MOO_1142_2asp_pca0_qm2_fitel_0f070-to-49f9Hz_1p0rr_L_dt20_snr_iter1.fits'
XMM_image_path = location / 'xray/images/XMM/XMM_adapt-400-7200.fits'
photz_member_smoothed_distribution_path = location / 'members/photz/MOO_1142.photz_wide_IntPz0.3_galdens.fits'
photz_members_path = location / 'members/photz/MOO_1142.photz_wide_IntPz0.3.fits'
specz_members_path = location / 'members/specz/MOO_1142+1527.specz_members.28aug23.fits'
# -------------------------------------------------------------------------------------------- #
# get data
photz_members_smoothed_hdu = fits.open(photz_member_smoothed_distribution_path)

photz_members = Table.read(photz_members_path)
photz_member_coords = np.array([photz_members['ra'], photz_members['dec']]).T

specz_members = Table.read(specz_members_path)
specz_members_zrange = specz_members[(specz_members['z'] > 1.0) & (specz_members['z'] < 1.4)]
specz_member_coords = np.array([specz_members_zrange['RA'], specz_members_zrange['Dec']]).T

def plot_data(base_image,show_specz=False):
    base_image_hdu = fits.open(base_image)
    base_image_data = base_image_hdu[0].data

    fig = plt.figure(figsize=(8,8))
    img = aplpy.FITSFigure(base_image_hdu,figure=fig)
    img.show_colorscale(cmap='RdBu_r')

    # plot properties
    img.recenter(x=moo1142_cen.ra.deg,y=moo1142_cen.dec.deg,radius=2.50/60)
    img.set_title('MOO 1142 M2 SNR')
    img.ticks.set_tick_direction('in')
    img.ticks.set_color('black')

    # additions
    clevels_neg=[-25,-20,-15,-10,-5]; clevels_pos=[5,10,15,20]
    img.show_contour(hdu=base_image_hdu,colors="white",linestyles="solid",levels=clevels_neg,returnlevels=True,zorder=1) #,smooth=3, convention='calabretta',
    img.show_contour(hdu=base_image_hdu,colors="gray",linestyles="dashed",levels=clevels_pos,returnlevels=True,zorder=2) #,smooth=3, convention='calabretta',

    # speczs
    if show_specz == True:
        img.show_markers(specz_member_coords[:,0],specz_member_coords[:,1],coords_frame='world',s=150,edgecolor='k',facecolor='yellow',marker='*',zorder=3)

    # plt.savefig(location / 'multi-wavelength_figures/M2_base/MOO_1142_baseM2.pdf', format='pdf', dpi=300)
    # plt.close()
    plt.show()
    
# -------------------------------------------------------------------------------------------- #
plot_data(base_image = M2_image_path,
          show_specz = True)
