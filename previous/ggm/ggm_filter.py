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

# M2 images
#mooinp = location+'ggm/M2_luca/data/MOOJ1142_snr_no_pnt_src.fits'
#mooinp = location+'M2/images/MOO_1142_M2_self-cal_map.fits'
mooinp = location+'M2/images/midas/2024-07-03_Kelvin_MOO_1142_2asp_pca0_qm2_fitel_0f070-to-49f9Hz_1p0rr_L_dt20_map_iter1.fits'

# Xray images
#mooinp = location+'xray/images/XMM/XMM_comb-net-image.fits'
#mooinp = location+'xray/images/XMM/XMM_comb-net-center.fits' - use this one!
#mooinp = location+'xray/images/XMM/XMM_comb-obj-im-400-7200.fits'
#mooinp = location+'xray/images/Chandra/Chandra_0.5-4_flux.img'

moohdu = fits.open(mooinp)
moodat = moohdu[0].data

for kernel in range(1,9):
  print('kernel size: {0} pix'.format(kernel))
  newfig = plt.figure(figsize=figsize)

  newggm = scipy.ndimage.gaussian_gradient_magnitude(moodat,kernel)

  newhdu = fits.PrimaryHDU(newggm,moohdu[0].header)
  # M2 images
  #newhdu.writeto(location+'ggm/M2/moo1142_M2.self-cal.ggm.{0}.fits'.format(kernel),overwrite=True)
  newhdu.writeto(location+'ggm/M2/2024-07/moo1142_M2_only_good.ggm.{0}.fits'.format(kernel),overwrite=True)
  # X-ray images
  #newhdu.writeto(location+'ggm/XMM/scipy/moo1142_xxm.combNetIm.ggm.{0}.fits'.format(kernel),overwrite=True)
  #newhdu.writeto(location+'ggm/XMM/scipy/moo1142_xxm.combNetCenter.ggm.{0}.fits'.format(kernel),overwrite=True)
  #newhdu.writeto(location+'ggm/XMM/scipy/moo1142_xxm.combObj.ggm.{0}.fits'.format(kernel),overwrite=True)
  #newhdu.writeto(location+'ggm/Chandra/moo1142_chandra_ggm.{0}.fits'.format(kernel),overwrite=True)

  newmin = np.nanmin(newhdu.data) if mspan else -np.nanmax(np.abs(newhdu.data))
  print(newmin)
  newmax = np.nanmax(newhdu.data) if mspan else  np.nanmax(np.abs(newhdu.data))
  print(newmax)

  newimg = aplpy.FITSFigure(data=newhdu,figure=newfig,convention='wells')
  newimg.show_colorscale(cmap='RdBu_r',vmin=newmin,vmax=newmax)

  newimg.recenter(x=175.6909647,y=15.4551255,radius=2.50/60)
  newimg.add_label(0.9,0.9,relative=True,text='$\sigma$='+str(kernel),color='white',size=18)

  newimg.ticks.set_tick_direction('in')
  newimg.ticks.set_color('black')

  # M2
  #plt.savefig(location+'ggm/M2/pngs/moo1142_M2.ggm.{0}.png'.format(kernel))
  #plt.savefig(location+'ggm/M2/pdf/moo1142_M2.self-cal.ggm.{0}.pdf'.format(kernel))
  plt.savefig(location+'ggm/M2/2024-07/pdfs/moo1142_M2_only_good.ggm.{0}.pdf'.format(kernel))
  # Xray
  #plt.savefig(location+'ggm/XMM/scipy/pdf/moo1142_xxm.combNetIm.ggm.{0}.pdf'.format(kernel))
  #plt.savefig(location+'ggm/XMM/scipy/pdf/moo1142_xxm.combNetCenter.ggm.{0}.pdf'.format(kernel))
  #plt.savefig(location+'ggm/XMM/scipy/pdf/moo1142_xxm.combNetObj.ggm.{0}.pdf'.format(kernel))
  #plt.savefig(location+'ggm/Chandra/pdf/moo1142_chandra_ggm.{0}.pdf'.format(kernel))

  #plt.show()
  plt.close()

  del newhdu