import matplotlib; matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.wcs import WCS

import aplpy

import numpy as np
import scipy.ndimage
from pathlib import Path

# ------------------------------------------------------------

dpi = 72.27*390.00/504.00

factorx = 0.60
factory = 0.50
figsize = (factorx*504.00/dpi,factory*504.00/dpi)

# ------------------------------------------------------------

mspan = True
clean = True

# ------------------------------------------------------------
location = f"{Path(__file__).resolve().parents[2]}/"
mooinp = location+'M2/images/MOO_1142_M2_self-cal_map.fits'

moohdu = fits.open(mooinp)

moodat = moohdu[0].data

for kernel in range(1,9):
  print('kernel size: {0} pix'.format(kernel))
  newfig = plt.figure(figsize=figsize)

  newggm = scipy.ndimage.gaussian_gradient_magnitude(moodat,kernel)

  newhdu = fits.PrimaryHDU(newggm,moohdu[0].header)
  newhdu.writeto(location+'ggm/M2/moo1142_M2.self-cal.ggm.{0}.fits'.format(kernel),overwrite=True)

  newmin = np.nanmin(newhdu.data) if mspan else -np.nanmax(np.abs(newhdu.data))
  newmax = np.nanmax(newhdu.data) if mspan else  np.nanmax(np.abs(newhdu.data))

  newimg = aplpy.FITSFigure(data=newhdu,figure=newfig,convention='wells')
  newimg.show_colorscale(cmap='RdBu_r',vmin=newmin,vmax=newmax)

  newimg.recenter(x=175.6909647,y=15.4551255,radius=2.50/60)
  newimg.ticks.set_tick_direction('in')
  newimg.ticks.set_color('black')

  plt.savefig(location+'ggm/M2/pdf/moo1142_M2.self-cal.ggm.{0}.pdf'.format(kernel))

  plt.close()

  del newhdu