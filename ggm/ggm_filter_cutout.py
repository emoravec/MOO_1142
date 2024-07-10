"""
2024-07
Script to make a ggm image and create a cutout of that image for plotting purposes.
"""
import matplotlib; matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from astropy.nddata import Cutout2D
from astropy import units as u

import aplpy

import numpy as np
import scipy.ndimage

# ------------------------------------------------------------
dpi = 72.27*390.00/504.00

factorx = 0.60
factory = 0.50
figsize = (factorx*504.00/dpi,factory*504.00/dpi)
# ------------------------------------------------------------
location = '/Users/emoravec/Documents/Research/merging_clusters/analysis/MOO_1142/'

# M2 images
mooinp = location+'M2/images/midas/2024-07-03_Kelvin_MOO_1142_2asp_pca0_qm2_fitel_0f070-to-49f9Hz_1p0rr_L_dt20_map_iter1.fits'

moohdu = fits.open(mooinp)
moodat = moohdu[0].data

for kernel in range(1,9):
  print('kernel size: {0} pix'.format(kernel))
  newfig = plt.figure(figsize=figsize)

  newggm = scipy.ndimage.gaussian_gradient_magnitude(moodat,kernel)

  newhdu = fits.PrimaryHDU(newggm,moohdu[0].header)
  # M2 images
  newhdu.writeto(location+'ggm/M2/2024-07/fits/moo1142_M2_only_good.ggm.{0}.fits'.format(kernel),overwrite=True)
  # # if already made cutout open it instead.
  # newhdu = fits.open(location+'ggm/M2/2024-07/moo1142_M2_only_good.ggm.{0}.fits'.format(kernel))[0]

  # create cutout for scaling purposes
  cen = SkyCoord(175.6909647,15.4551255,unit='deg')
  img_size = u.Quantity((6.0, 6.0), u.arcmin)
  new_wcs = WCS(newhdu.header)
  # Make the cutout, including the WCS
  cutout = Cutout2D(newhdu.data, position=cen, size=img_size, wcs=new_wcs)
  # Put the cutout image in the FITS HDU
  newhdu.data = cutout.data
  # Update the FITS header with the cutout WCS
  newhdu.header.update(cutout.wcs.to_header())
  # Write the cutout to a new FITS file
  newhdu.writeto(location+'ggm/M2/2024-07/fits/moo1142_M2_only_good.ggm.{0}.cutout.fits'.format(kernel),overwrite=True)

  # plot the new cutout
  newimg = aplpy.FITSFigure(data=newhdu,figure=newfig,convention='wells')
  newimg.show_colorscale(cmap='RdBu_r')

  newimg.recenter(x=175.6909647,y=15.4551255,radius=2.50/60)
  newimg.add_label(0.9,0.9,relative=True,text='$\sigma$='+str(kernel),color='white',size=18)

  newimg.ticks.set_tick_direction('in')
  newimg.ticks.set_color('black')

  # M2
  plt.savefig(location+'ggm/M2/2024-07/pdfs/moo1142_M2_only_good.ggm.{0}.cutout.pdf'.format(kernel))

  #plt.show()
  plt.close()

  del newhdu