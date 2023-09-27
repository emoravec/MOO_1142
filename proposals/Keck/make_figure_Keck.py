"""

Created in 2023-09
"""
import matplotlib.pyplot as plt
import matplotlib; matplotlib.use('TkAgg')
matplotlib.rc('font',**{'family':'serif','sans-serif':['Times'],'size':11})
matplotlib.rc('text',usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{newtxtext}')
plt.rc('text.latex', preamble=r'\usepackage{newtxmath}')
plt.rc('font',family='serif',size=14)
plt.rc('text',usetex=True)

from palettable.cartocolors.sequential import Mint_7 as pal4
cmap = pal4.mpl_colormap

import numpy as np
import pandas as pd
from pathlib import Path

import aplpy
from astropy.io import fits
from astropy.table import Table,vstack
from astropy.coordinates import SkyCoord
from astropy import units as u
# ------------------------------------------------------------------------------------------- #
dpi = 72.27*390.00/504.00

factorx = 1.0 # 0.60
factory = 0.8 # 0.50
figsize = (factorx*504.00/dpi,factory*504.00/dpi)

location = Path('/Users/emoravec/Documents/Research/merging_clusters/')
chandra_img_values = Table.read(location / 'tables/Chandra_sample/Chandra_min_max.csv',format='ascii.csv')
chandra_img_values.add_index('Name')

location = Path('/Users/emoravec/Documents/Research/merging_clusters/')
# ------------------------------------------------------------------------------------------- #
def get_opt_ir_data(cluster):
    # get smoothed density distribution
    smoothed_distribution_filename = location / f'analysis/vla21b/density_dist/galdens_output/IntPz/{cluster}_galdens_IntPz_gt0.3.fits'
    hdu = fits.open(smoothed_distribution_filename)
    # get members
    table_path = location / f'catalogs/vla21b/cluster_members/IntPz/{cluster}_IntPz0.3.fits'
    members = Table.read(table_path)
    member_coords = np.array([members['ra'], members['dec']]).T
    return hdu,member_coords

def make_xray_contour_levels_lin(min,max):
    num_contours = 10
    contours = np.linspace(min,max,num_contours)
    return contours

def make_xray_contour_levels_log(min,max):
    num_contours = 10
    contours = np.logspace(np.log10(min),np.log10(max),num_contours)
    return contours

def make_keck_plot(hdu,coord_array,specz_coord_array,cluster,image_size):
    fig = plt.figure(figsize=figsize)

    # make units of data to be #/arcmin^2
    density_data, density_head = hdu[0].data, hdu[0].header
    # want to convert density_data to arcmin^2 so convert data from #/pix to #/arcmin
    # cdelt = deg/sqrt(pix)
    # #/pix * (sqrt(pix)/arcmin)^2 = #/arcmin^2 = data * (1/(cdelt1 * 60 * cdelt2 * 60))
    density_data *= 1.0/np.abs(density_head['CDELT1']*60.00)/np.abs(density_head['CDELT2']*60.00)

    # image parameters
    density_img = aplpy.FITSFigure(data=fits.PrimaryHDU(data=density_data,header=density_head),figure=fig)
    center = SkyCoord('11:42:47','15:27:00',unit=(u.hourangle,u.deg))
    density_img.recenter(center.ra.deg,center.dec.deg,radius=image_size/60.0)
    
    # contour levels
    density_min = np.round(density_data.min()-0.50) # adding a shift of 0.5 up or down to ensure encapsulates range -> round(1.1+0.5) -> round(1.6) -> 2
    density_max = np.round(density_data.max()+0.50)
    density_dlt=1.00
    density_lev = np.linspace(1,density_max,1+int((density_max-1)/density_dlt)) # have +1 to make sure it includes the last level

    # show contours
    density_img.show_contour(data=fits.PrimaryHDU(data=density_data,header=density_head),cmap=cmap,filled=True,levels=density_lev,extend='max')
    density_con = density_img._layers['contour_set_1']
    density_con.cmap.set_under('white')

    # show cluster members
    density_img.show_markers(coord_array[:,0],coord_array[:,1],coords_frame='world',s=40,alpha=0.50,edgecolor='none',facecolor='black')

    # show SPECTROSCOPIC cluster members
    density_img.show_markers(specz_coord_array[:,0],specz_coord_array[:,1],coords_frame='world',s=150,edgecolor='k',facecolor='salmon',marker='*')

    # plot parameters
    density_img.set_title(cluster.replace('_',' '))
    density_img.ax.tick_params(axis='both',which='both',direction='in',color='k',length=7)
    density_img.ax.tick_params(which='minor', length=4)

    # colorbar
    density_bar = plt.colorbar(density_con)
    density_bar.set_label(r'Galaxy density [$\mathrm{arcmin^{-2}}$]',labelpad=10)
    #density_bar.ax.set_yticklabels(['{0:2d}'.format(int(np.round(float(tick.get_text()[1:-1])))) for tick in density_bar.ax.yaxis.get_ticklabels()])

    # X-ray contours
    xray_img = f'/Users/emoravec/Documents/Research/merging_clusters/images/Chandra/smooth_masked/{cluster}_smooth_masked_cutout_5am.fits'
    xray_contour_levels = make_xray_contour_levels_log(chandra_img_values.loc[c]['min'],chandra_img_values.loc[c]['max'])
    levels = density_img.show_contour(xray_img,returnlevels=True,levels=xray_contour_levels,smooth=1,colors='black',linestyles='-',linewidths=1,zorder=1)

    # save plot
    plt.savefig(f'/Users/emoravec/Documents/Research/merging_clusters/analysis/MOO_1142/proposals/Keck/{cluster}_opt-ir-dist_Chandra_contours.pdf',format='pdf',dpi=300)
# ------------------------------------------------------------------------------------------- #
cluster_info = Table.read(location / 'tables/high_z/vla21b_sample_MC.csv',format='ascii.csv')
cluster_info.add_index('Name')
c = 'MOO_1142+1527'

hdu,member_coords= get_opt_ir_data(c)
specz_table_path = '/Users/emoravec/Documents/Research/merging_clusters/analysis/MOO_1142/members/MOO_1142+1527.specz_members.28aug23.fits'
specz_members = Table.read(specz_table_path)
specz_member_coords = np.array([specz_members['RA'], specz_members['Dec']]).T

make_keck_plot(hdu=hdu,
                coord_array=member_coords,
                specz_coord_array=specz_member_coords,
                cluster=c,
                image_size=3.0)