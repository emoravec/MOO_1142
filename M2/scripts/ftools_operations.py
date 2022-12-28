"""
Calls various FTOOLS routines to do unsharp masking.
Created in 2022-12
"""
import subprocess
import matplotlib.pyplot as plt
import numpy as np

from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
# -------------------------------------------------------------------------------------------- #
location = '/Users/emoravec/Documents/Research/merging_clusters/'
M2_loc = 'analysis/MOO_1142/M2/'
# # determining pixel scale
# hdulist = fits.open(location+'images/M2/moo1142_snr_pt_source_sub.fits')
# header = hdulist[0].header
# cdelt2 = header['CDELT2']*u.deg
# cdelt2_arcs = np.round(cdelt2.to(u.arcsec).value,2) #-> 1

# smooth the SNR point source subtracted image on the order of the M2 beam
# FWHM = 2.3548*sigma where FWHM of M2 = 9"
# sigma_arcsec = 9/2.3548
# smooth_scale_M2_pix = np.rint(sigma_arcsec) # arcsec * pixel/(deg->arcsec) -> np.round((sigma_arcsec/cdelt2_arcs),2) -> but pixelscale is 1"/pix 
smooth_scale_M2_pix = 5
# image paths
input_image = location + 'images/M2/moo1142_snr_pt_source_sub.fits'
M2_smoothed_image = location + M2_loc + 'smoothed/MOO_1142_smoothed_M2Beam_'+str(smooth_scale_M2_pix)+'pix.fits'
# use FTOOLS fgauss to do the smoothing
subprocess.call(["fgauss",input_image,M2_smoothed_image,str(smooth_scale_M2_pix)])

# smooth image to many different large scale sizes
scale_array = [25,50,100] # roughly 5,10,15 x PSF which is 9"/2.35~4

print('Smoothing and subtracting.')
for s in scale_array:
    print(s)
    # use FTOOLS fgauss to do the smoothing
    smoothed_image = location + M2_loc + 'smoothed/MOO_1142_smoothed_scale_'+str(s)+'pix.fits'
    subprocess.call(["fgauss",input_image,smoothed_image,str(s)])
    # subtract the M2 beam image from the large-scale image
    subtracted_image = location + M2_loc + 'subtracted/MOO_1142_sub_scale'+str(s)+'pix.fits'
    subprocess.call(["farith", smoothed_image, M2_smoothed_image, subtracted_image, "SUB"])