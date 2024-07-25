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
specz_members_path = '/Users/emoravec/Documents/Research/merging_clusters/analysis/MOO_1142/members/MOO_1142+1527.specz_members.28aug23.fits'
# -------------------------------------------------------------------------------------------- #
# get data
photz_members_smoothed_hdu = fits.open(photz_member_smoothed_distribution_path)

photz_members = Table.read(photz_members_path)
photz_member_coords = np.array([photz_members['ra'], photz_members['dec']]).T

specz_members = Table.read(specz_members_path)
specz_member_coords = np.array([specz_members['RA'], specz_members['Dec']]).T