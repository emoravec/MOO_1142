import matplotlib
matplotlib.rc('font',**{'family':'serif','sans-serif':['Times'],'size':11})
matplotlib.rc('text',usetex=True)

import matplotlib.pyplot as plt
plt.rc('text.latex', preamble=r'\usepackage{newtxtext}')
plt.rc('text.latex', preamble=r'\usepackage{newtxmath}')

from palettable.cartocolors.sequential import Mint_7 as pal4
cmap = pal4.mpl_colormap
#######################################################################################################################
import numpy as np
import pandas as pd
from pathlib import Path

from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table

import aplpy
#######################################################################################################################
dpi = 72.27*390.00/504.00

factorx = 1.0 # 0.60
factory = 0.8 # 0.50
figsize = (factorx*504.00/dpi,factory*504.00/dpi)
#######################################################################################################################
def plot_galdens_basic(hdu,coord_array,image_size):
    fig = plt.figure(figsize=figsize)

    # units are already in #/arcmin^2
    density_data, density_head = hdu[0].data, hdu[0].header

    # image parameters
    density_img = aplpy.FITSFigure(data=fits.PrimaryHDU(data=density_data,header=density_head),figure=fig)
    density_img.recenter(x=density_head['CRVAL1'],y=density_head['CRVAL2'],radius=image_size/60.00)
    
    # # contour levels - 0.5
    # density_max = density_data.max()
    # density_dlt = 0.5
    # density_lev = np.arange(1,density_max+density_dlt,density_dlt)

    # contour levels - 1.0
    density_max = density_data.max()
    density_dlt = 1.00
    density_lev = np.arange(1,density_max+density_dlt,density_dlt)

    # show contours
    density_img.show_contour(data=fits.PrimaryHDU(data=density_data,header=density_head),cmap=cmap,filled=True,levels=density_lev,extend='max')
    density_con = density_img._layers['contour_set_1']
    density_con.cmap.set_under('white')

    # plot parameters
    density_img.set_title('MOO 1142+1527 Pyrostat Specz Members')
    density_img.ax.tick_params(axis='both',which='both',direction='in',color='k')

    # colorbar
    density_bar = plt.colorbar(density_con)
    density_bar.set_label(r'Galaxy density [$\mathrm{arcmin^{-2}}$]',labelpad=10)

    # show data
    density_img.show_markers(coord_array[:,0],coord_array[:,1],coords_frame='world',s=20,alpha=0.50,edgecolor='none',facecolor='black')

    return density_img
#######################################################################################################################
# read in the galaxy density map
location = Path('/Users/emoravec/Documents/Research/merging_clusters/analysis/MOO_1142/members')
galaxy_density_path = location / 'photz_MC1_cats/MOO_1142+1527_R2_galdens_scaled_IntPz_gt0.3.fits'
galaxy_density_hdu = fits.open(galaxy_density_path)
members = Table.read(location / 'photz_MC1_cats/MOO_1142_R2_photz_members_IntPzgt03.fits')
# make coordinates array
coords = np.array([members['RA'], members['DEC']]).T
# Make galdens basic plot
density_img = plot_galdens_basic(hdu=galaxy_density_hdu,
                                 coord_array=coords,
                                 image_size=4)
# Add open face star for BCG location
density_img.show_markers(175.69784, 15.45318, coords_frame='world', s=100, marker='*', facecolor='none', edgecolor='yellow', linewidth=1)
# # Save the plot of just galaxy density map with photz members
# plt.savefig(location / 'photz_MC1_cats/MOO_1142+1527_R2_galdens_contours0p5.pdf',format='pdf',dpi=300)

### Add spec-zs
# # all good quality speczs that we have
# all_speczs = Table.read(location / 'specz/MOO_1142+1527.speczs.30jul25.goodQ.fits')
# density_img.show_markers(all_speczs['RA'], all_speczs['Dec'], coords_frame='world', s=100, marker='*', facecolor='none', edgecolor='magenta', linewidth=1)
# # Save the plot of photz members and spec-z members
# plt.savefig(location / 'plots/photz_specz/MOO_1142+1527_R2_galdens_all_good_specz.pdf',format='pdf',dpi=300)

# only members
with open(location / 'specz/MOO_1142+1527.speczs.goodQ.pyrostat.members.14aug25.txt', 'r') as f:
    first_line = f.readline().strip()
    column_names = first_line.replace('#', '').strip().split()
    
# Read the data skipping the comment line
specz_members = pd.read_csv(location / 'specz/MOO_1142+1527.speczs.goodQ.pyrostat.members.14aug25.txt', 
                           delim_whitespace=True, 
                           comment='#',
                           names=column_names)
density_img.show_markers(specz_members['RA'], specz_members['Dec'], coords_frame='world', s=100, marker='*', facecolor='none', edgecolor='black', linewidth=1)
plt.savefig(location / 'plots/photz_specz/MOO_1142+1527_R2_galdens_specz_members.pdf',format='pdf',dpi=300)

