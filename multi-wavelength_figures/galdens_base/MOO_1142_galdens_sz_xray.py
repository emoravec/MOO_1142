"""
Plots the galaxy density, overlaid with X-ray (XMM net-comb-center) and SZ (minkasi signal-smoothed-10arcsec) contours.

Created in 2026-04
"""
import matplotlib.pyplot as plt
plt.rc('font', **{'family': 'serif', 'sans-serif': ['Times'], 'size': 11})
plt.rc('text', usetex=True)

from palettable.cartocolors.sequential import Mint_7 as pal4
cmap = pal4.mpl_colormap
#######################################################################################################################
import numpy as np
from pathlib import Path

from astropy.io import fits
from astropy.table import Table
from astropy import units as u
from astropy.coordinates import SkyCoord

import aplpy
#######################################################################################################################
dpi = 72.27 * 390.00 / 504.00

factorx = 1.0
factory = 0.8
figsize = (factorx * 504.00 / dpi, factory * 504.00 / dpi)
location = Path(__file__).resolve().parents[2]
by_eye_center = SkyCoord('11:42:46.45 15:27:12.38', unit=(u.hourangle,u.deg))

# Choose which M2 product to contour: 'signal' or 'snr'
M2_IMAGE_OPTION = 'snr'
#######################################################################################################################
def plot_galdens_basic(hdu, coord_array, image_size):
    fig = plt.figure(figsize=figsize)

    # units are already in #/arcmin^2
    density_data, density_head = hdu[0].data, hdu[0].header

    density_img = aplpy.FITSFigure(
        data=fits.PrimaryHDU(data=density_data, header=density_head),
        figure=fig,
        convention='calabretta',
    )
    density_img.recenter(by_eye_center.ra, by_eye_center.dec, radius=image_size / 60.00)

    density_max = density_data.max()
    density_dlt = 1.00
    density_lev = np.arange(1, density_max + density_dlt, density_dlt)

    density_img.show_contour(
        data=fits.PrimaryHDU(data=density_data, header=density_head),
        cmap=cmap,
        filled=True,
        levels=density_lev,
        extend='max',
        convention='calabretta',
    )
    density_con = density_img._layers['contour_set_1']
    density_con.cmap.set_under('white')

    density_img.set_title('MOO 1142+1527')
    density_img.ax.tick_params(axis='both', which='both', direction='in', color='k')

    density_bar = plt.colorbar(density_con)
    density_bar.set_label(r'Galaxy density [$\mathrm{arcmin^{-2}}$]', labelpad=10)

    density_img.show_markers(
        coord_array[:, 0],
        coord_array[:, 1],
        coords_frame='world',
        s=20,
        alpha=0.50,
        edgecolor='none',
        facecolor='black',
    )

    return density_img

def make_contour_levels_sqrt(min_value, max_value, num_contours):
    return np.square(np.linspace(np.sqrt(min_value), np.sqrt(max_value), num_contours))

def make_contour_levels_linear(min_value, max_value, num_contours):
    return np.linspace(min_value, max_value, num_contours)
#######################################################################################################################
# Galaxy Density
galaxy_density_path = location / 'members/photz_MC1_cats/MOO_1142+1527_R2_galdens_scaled_IntPz_gt0.3.fits'
galaxy_density_hdu = fits.open(galaxy_density_path)
members = Table.read(location / 'members/photz_MC1_cats/MOO_1142_R2_photz_members_IntPzgt03.fits')
coords = np.array([members['RA'], members['DEC']]).T

# Base Plot
density_img = plot_galdens_basic(
    hdu=galaxy_density_hdu,
    coord_array=coords,
    image_size=2,
)

# BCG
density_img.show_markers(
    175.69784,
    15.45318,
    coords_frame='world',
    s=100,
    marker='*',
    facecolor='none',
    edgecolor='yellow',
    linewidth=1,
)

# X-ray Contours
xmm_img = location / 'xray/images/XMM/smoothed/XMM_comb-net-center_smoothed_10arcsec.fits'
xmm_contour_levels = make_contour_levels_sqrt(7.5e-6, 1.2e-4, 12)
density_img.show_contour(
    str(xmm_img),
    returnlevels=True,
    levels=xmm_contour_levels,
    smooth=1,
    convention='calabretta',
    colors='black',
    linestyles='-',
    linewidths=1,
    zorder=1,
)

# SZ Contours
if M2_IMAGE_OPTION == 'signal':
    m2_img = location / 'M2/images/minkasi/2025-11_JackOS/smoothed/MOO_1142_signal_smoothed_10arcsec.fits'
    m2_contour_levels = make_contour_levels_linear(-5.0e-4, 2.5e-4, 13)
    output_name = 'MOO_1142+1527_base=galdens_contours=xmm_m2_minkasi_zoom.png'
elif M2_IMAGE_OPTION == 'snr':
    m2_img = location / 'M2/images/minkasi/2025-11_JackOS/cutouts/MOO_1142_snrmap_4p5_by_eye_ctr.fits'
    m2_contour_levels = np.arange(0, 12, 1)
    output_name = 'MOO_1142+1527_base=galdens_contours=xmm_m2_snr_zoom.png'
else:
    raise ValueError("M2_IMAGE_OPTION must be either 'signal' or 'snr'.")

density_img.show_contour(
    str(m2_img),
    returnlevels=True,
    levels=m2_contour_levels,
    smooth=1,
    convention='calabretta',
    colors='red',
    linestyles='-',
    linewidths=1,
    zorder=1,
)

plt.savefig(location / 'multi-wavelength_figures/galdens_base/pngs' / output_name, format='png', dpi=300)
#plt.show()