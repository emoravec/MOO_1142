"""
Makes cutouts from M2 images.

@author: emoravec; created 2026-04-14
"""
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy import units as u
from astropy.wcs import WCS
import matplotlib.pyplot as plt
from pathlib import Path

location = Path(__file__).resolve().parent

# image center - from /Users/emoravec/Documents/Research/MOO_1142/profile_fitting/XMM/XMM_fit_coords_EBarbavara.txt
# which is Eleonora's fit to the main and west subcluster
main = SkyCoord(175.69766238,15.45364413,unit='deg')
west = SkyCoord(175.68837814,15.45566318,unit='deg')
by_eye_center = SkyCoord('11:42:46.45 15:27:12.38', unit=(u.hourangle,u.deg))
img_size = u.Quantity((4.5, 4.5), u.arcmin)

hdu = fits.open(location / '2025-11_JackOS/MOO_1142_snrmap.fits')[0]
wcs = WCS(hdu.header)

# Make the cutout, including the WCS
cutout = Cutout2D(hdu.data, position=by_eye_center, size=img_size, wcs=wcs)
#plt.imshow(cutout.data, origin='lower')   

# Put the cutout image in the FITS HDU
hdu.data = cutout.data

# Update the FITS header with the cutout WCS
hdu.header.update(cutout.wcs.to_header())

# Write the cutout to a new FITS file
cutout_filename = location / '2025-11_JackOS/cutouts/MOO_1142_snrmap_4p5_by_eye_ctr.fits'
hdu.writeto(cutout_filename, overwrite=True)