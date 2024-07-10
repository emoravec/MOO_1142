import matplotlib; matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.wcs import WCS

import aplpy

import numpy as np
import scipy.ndimage

# ------------------------------------------------------------

dpi = 72.27*390.00/504.00

factorx = 0.60
factory = 0.50
figsize = (factorx*504.00/dpi,factory*504.00/dpi)

# ------------------------------------------------------------

mspan = True
clean = True

# ------------------------------------------------------------

mooinp = 'data/MOOJ1142_snr_no_pnt_src.fits' if clean else \
         'data/MOOJ1142_first_map_precon_mpi_EigSmooth_rebinpowers_model_files0-38_bunch3.fits'

moohdu = fits.open(mooinp)

moodat = moohdu[0].data

for kernel in range(1,10):
  print('kernel size: {0} pix'.format(kernel))
  newfig = plt.figure(figsize=figsize)

  newggm = scipy.ndimage.gaussian_gradient_magnitude(moodat,kernel)

  newhdu = fits.PrimaryHDU(newggm,moohdu[0].header)
  newhdu.writeto('output/moo1142_01_ggm.{0}.fits'.format(kernel),overwrite=True)

  newmin = np.nanmin(newhdu.data) if mspan else -np.nanmax(np.abs(newhdu.data))
  newmax = np.nanmax(newhdu.data) if mspan else  np.nanmax(np.abs(newhdu.data))

  newimg = aplpy.FITSFigure(data=newhdu,figure=newfig,convention='wells')
  newimg.show_colorscale(cmap='RdBu_r',vmin=newmin,vmax=newmax)

  newimg.recenter(x=175.6909647,y=15.4551255,radius=2.50/60)

  newimg.ticks.set_tick_direction('in')
  newimg.ticks.set_color('black')

  plt.savefig('images/moo1142_ggm.{0}.pdf'.format(kernel))
  plt.close()

  del newhdu

