"""
Plots the galaxy density, overlaid with X-ray and SZ contours.

Created in 2025-10
"""
import matplotlib.pyplot as plt
plt.rc('font',**{'family':'serif','sans-serif':['Times'],'size':11})
plt.rc('text',usetex=True)

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
location = Path(__file__).resolve().parent.parent
#######################################################################################################################
def plot_galdens_basic(hdu,coord_array,image_size):
    fig = plt.figure(figsize=figsize)

    # units are already in #/arcmin^2
    density_data, density_head = hdu[0].data, hdu[0].header

    # image parameters
    density_img = aplpy.FITSFigure(data=fits.PrimaryHDU(data=density_data,header=density_head),figure=fig)
    density_img.recenter(x=density_head['CRVAL1'],y=density_head['CRVAL2'],radius=image_size/60.00)
    
    # contour levels - 1.0
    density_max = density_data.max()
    density_dlt = 1.00
    density_lev = np.arange(1,density_max+density_dlt,density_dlt)

    # show contours
    density_img.show_contour(data=fits.PrimaryHDU(data=density_data,header=density_head),cmap=cmap,filled=True,levels=density_lev,extend='max')
    density_con = density_img._layers['contour_set_1']
    density_con.cmap.set_under('white')

    # plot parameters
    density_img.set_title('MOO 1142+1527')
    density_img.ax.tick_params(axis='both',which='both',direction='in',color='k')

    # colorbar
    density_bar = plt.colorbar(density_con)
    density_bar.set_label(r'Galaxy density [$\mathrm{arcmin^{-2}}$]',labelpad=10)

    # show data
    density_img.show_markers(coord_array[:,0],coord_array[:,1],coords_frame='world',s=20,alpha=0.50,edgecolor='none',facecolor='black')

    return density_img

def make_xray_contour_levels_log(min,max):
    num_contours = 10
    contours = np.logspace(np.log10(min),np.log10(max),num_contours)
    return contours
#######################################################################################################################
### Galaxy Density
galaxy_density_path = location / 'members/photz_MC1_cats/MOO_1142+1527_R2_galdens_scaled_IntPz_gt0.3.fits'
galaxy_density_hdu = fits.open(galaxy_density_path)
members = Table.read(location / 'members/photz_MC1_cats/MOO_1142_R2_photz_members_IntPzgt03.fits')
# make coordinates array
coords = np.array([members['RA'], members['DEC']]).T
# Make galdens basic plot
density_img = plot_galdens_basic(hdu=galaxy_density_hdu,
                                 coord_array=coords,
                                 image_size=4)
### BCG
density_img.show_markers(175.69784, 15.45318, coords_frame='world', s=100, marker='*', facecolor='none', edgecolor='yellow', linewidth=1)

### X-ray Contours
xmm_img = location / 'xray/images/XMM/XMM_adapt-400-7200.fits'
xmm_contour_levels = make_xray_contour_levels_log(10,320)
xmm_levels = density_img.show_contour(str(xmm_img),returnlevels=True,levels=xmm_contour_levels,smooth=1,colors='black',linestyles='-',linewidths=1,zorder=1)

### SZ
m2_img = location / 'M2/images/midas/2024-07-03_Kelvin_MOO_1142_2asp_pca0_qm2_fitel_0f070-to-49f9Hz_1p0rr_L_dt20_snr_iter1.fits'
m2_contour_levels = np.array([-25,-20,-15,-10,-5,-3,3,5,10])
m2_levels = density_img.show_contour(str(m2_img),returnlevels=True,levels=m2_contour_levels,smooth=1,colors='red',linestyles='--',linewidths=1,zorder=1)

plt.savefig(location / 'multi-wavelength_figures/galdens_base/MOO_1142+1527_galdens_xmm_m2-pt-src.png',format='png',dpi=300)
