"""
Plots the galdensity map with spec-z members overlaid. The color and size of the dot corresponds to the delta velocity of each member.

Created in 2024-10 by E. Moravec
"""
import matplotlib
matplotlib.rc('font',**{'family':'serif','sans-serif':['Times'],'size':11})
matplotlib.rc('text',usetex=True)

import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm

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
def plot_galdens_basic(hdu,coord_array,image_size):
    fig = plt.figure(figsize=(12,8))

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
    density_img.set_title('MOO 1142+1527 Members')
    density_img.ax.tick_params(axis='both',which='both',direction='in',color='k')

    # colorbar
    density_bar = plt.colorbar(density_con, pad=0.02)  # Reduce from default (~0.05) to 0.02
    density_bar.set_label(r'Galaxy density [$\mathrm{arcmin^{-2}}$]', labelpad=10)

    # show data
    density_img.show_markers(coord_array[:,0],coord_array[:,1],coords_frame='world',s=20,alpha=0.50,edgecolor='none',facecolor='black')

    return density_img

def plot_members_galdens_velocity(density_img, specz_members, velocity_range=(-3000, 3000)):
    """
    Add spectroscopic members to galaxy density plot with size and color based on Delta_v.
    
    Parameters:
    -----------
    density_img : aplpy.FITSFigure
        The galaxy density plot
    specz_members : pandas.DataFrame
        DataFrame with RA, Dec, and Delta_v columns
    velocity_range : tuple
        Min and max velocity for color scaling
    """
    
    # Calculate sizes based on absolute Delta_v (larger for higher velocities)
    abs_delta_v = np.abs(specz_members['Delta_v'])
    # Scale sizes from 50 to 200 based on velocity range
    min_size, max_size = 50, 200
    if abs_delta_v.max() > abs_delta_v.min():
        sizes = min_size + (max_size - min_size) * (abs_delta_v - abs_delta_v.min()) / (abs_delta_v.max() - abs_delta_v.min())
    else:
        sizes = np.full(len(specz_members), min_size)
    
    # Create velocity-based colormap plot with SymLogNorm
    vmin, vmax = velocity_range
    norm = SymLogNorm(linthresh=200, linscale=1, vmin=vmin, vmax=vmax)
    cmap_vel = plt.cm.RdBu_r  # Red for positive, Blue for negative (add _r back)
    
    # Plot each point individually to control size and color
    for i, row in specz_members.iterrows():
        color = cmap_vel(norm(row['Delta_v']))
        density_img.show_markers(row['RA'], row['Dec'], 
                                coords_frame='world', 
                                s=sizes[i], 
                                facecolor=color,
                                edgecolor='black', 
                                linewidth=0.5,
                                alpha=0.8,
                                marker='o')
    
    # Add velocity colorbar
    # Create a ScalarMappable for the colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap_vel, norm=norm)
    sm.set_array([])
    
    # Add colorbar to the right side
    cbar = plt.colorbar(sm, ax=density_img.ax, shrink=0.6, pad=0.02)  # Reduce from 0.12 to 0.02
    cbar.set_label(r'$\Delta_v$ [km/s]', labelpad=10)
    
    # Customize the colorbar tick labels to show regular numbers instead of scientific notation
    cbar.ax.tick_params(labelsize=10)
    # Set custom tick locations and labels
    tick_locations = [-3000, -1000, -100, 0, 100, 1000, 2000]
    tick_labels = ['-3000', '-1000', '-100', '0', '100', '1000', '2000']
    cbar.set_ticks(tick_locations)
    cbar.set_ticklabels(tick_labels)
    
    return density_img
#######################################################################################################################
# read in the galaxy density map
location = Path(__file__).resolve().parent
galaxy_density_path = location / 'photz_MC1_cats/MOO_1142+1527_R2_galdens_scaled_IntPz_gt0.3.fits'
galaxy_density_hdu = fits.open(galaxy_density_path)
members = Table.read(location / 'photz_MC1_cats/MOO_1142_R2_photz_members_IntPzgt03.fits')

# make coordinates array for photz members
coords = np.array([members['RA'], members['DEC']]).T

# Make galdens basic plot with photz members
density_img = plot_galdens_basic(hdu=galaxy_density_hdu,
                                 coord_array=coords,
                                 image_size=4)

# Add open face star for BCG location
density_img.show_markers(175.69784, 15.45318, coords_frame='world', 
                        s=100, marker='*', facecolor='none', 
                        edgecolor='yellow', linewidth=2)

# Read in the specz members with delta velocity information
specz_members = pd.read_csv(location / 'specz/MOO_1142+1527.speczs.goodQ.pyrostat.members.with_delta_v.25oct25.csv')

# Add specz members with velocity-based size and color
density_img = plot_members_galdens_velocity(density_img, specz_members, velocity_range=(-3000, 2000))

# Save the plot
plt.savefig(location / 'plots/photz_specz/MOO_1142+1527_R2_galdens_specz_members_velocity.pdf',
           format='pdf', dpi=300, bbox_inches='tight')
# plt.show()