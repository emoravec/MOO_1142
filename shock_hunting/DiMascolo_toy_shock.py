import matplotlib; matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70.00,Om0=0.30)

from astropy import units
from astropy import constants as const

import numpy as np
import scipy.integrate

# Support function
# ======================================================================
def hypot3d(x,y,z):
  return np.sqrt(x*x+y*y+z*z)

def jitgrid(gridx,gridy,cntrx,cntry,sint,cost,cosy,e,rs):
  modgridx = (-(gridx-cntrx)*cosy*sint-(gridy-cntry)*cost)
  modgridy = ( (gridx-cntrx)*cosy*cost-(gridy-cntry)*sint)/np.sqrt(1.00-e**2)
  return modgridx*(np.pi/1.8e2)/rs, modgridy*(np.pi/1.8e2)/rs, np.sqrt((modgridx**2)+(modgridy**2))*(np.pi/1.8e2)/rs


# Models
# ======================================================================
# Line-of-sight eccentricity
# ----------------------------------------------------------------------
def elos(e): return e/np.sqrt(2.00-e*e)

# A10 radial profile
# ----------------------------------------------------------------------
def a10ProfileRadial(x,alpha,beta,gamma,ap,c500,mass): 
  return (x**(-gamma))*((1.00+(x**alpha))**((gamma-beta)/alpha))*(mass**((ap+0.10)/(1.00+(2.00*x/c500)**3.00)))

# 3D A10 model profile
# ----------------------------------------------------------------------
def a10ProfileIntegrand(x,xi,alpha,beta,gamma,ap,c500,mass): 
  return (x**(-gamma))*((1.00+(x**alpha))**((gamma-beta)/alpha))*(x/((x*x-xi*xi)**0.50))*(mass**((ap+0.10)/(1.00+(2.00*x/c500)**3.00)))

# A10 model integral
# ----------------------------------------------------------------------
def _a10ProfileIntegral(x,alpha,beta,gamma,ap,c500,mass,limdist=np.inf,epsrel=1.00E-06): 
  return 2.00*scipy.integrate.quad(a10ProfileIntegrand,x,limdist,args=(x,alpha,beta,gamma,ap,c500,mass),epsrel=epsrel)[0]
a10ProfileIntegral = np.vectorize(_a10ProfileIntegral)

# Integrated elliptical A10 model profile
# ----------------------------------------------------------------------
def a10Profile(grid,offset,amp,major,e,alpha,beta,gamma,ap,c500,mass,limdist=np.inf,epsrel=1.00E-06,freeLS=None): 
  integral = np.zeros_like(grid,dtype=np.float64)
  if isinstance(grid,np.ndarray):
    integral[grid<=limdist] = a10ProfileIntegral(grid[grid<=limdist],alpha,beta,gamma,ap,c500,mass,limdist,epsrel)
  else:
    integral = a10ProfileIntegral(grid,alpha,beta,gamma,ap,c500,mass,limdist,epsrel)
  ellipse = np.sqrt(1.00-elos(e)**2) if freeLS is None else freeLS
  return offset+amp*major*ellipse*integral


# Constants
# ======================================================================
fb = 0.175; mu = 0.590; mue = 1.140

gadab = 5.00/3.00

ysznorm = const.sigma_T/const.m_e/const.c**2
ysznorm = ysznorm.to(units.cm**3/units.keV/units.Mpc)


# Main
# ======================================================================

modx = 11.00
mody = 20.00

wline = np.linspace(-2.00/60.00,2.00/60.00,256)
zline = np.linspace( 0.00/60.00,2.00/60.00,512)
wgrid = np.meshgrid(wline,wline)

wgrid[0] = modx+wgrid[0]/np.cos(np.deg2rad(mody))
wgrid[1] = mody+wgrid[1]
wcube = np.array([np.broadcast_to(wgrid[0][:,:,None],(*wgrid[0].shape,zline.shape[0])),
                  np.broadcast_to(wgrid[1][:,:,None],(*wgrid[0].shape,zline.shape[0])),
                  np.broadcast_to(zline[None,None,:],(*wgrid[0].shape,zline.shape[0]))])

alpha = 1.0510E+00; beta = 5.4905E+00; gamma =  3.0810E-01
pnorm = 8.4030E+00; c500 = 1.1770E+00; ap    =  1.2000E-01
pnorm = pnorm*((cosmo.H0.value/70.00)**(-3.00/2.00))

e = 0.75
bias = 0.00E+00
zeta = 1.00E+00
m500 = np.log10(3.00E+14)

theta = 80.00

Hz = cosmo.H(zeta)
dz = cosmo.angular_diameter_distance(zeta)
         
r500 = ((3.00/4.00/np.pi/500.00/cosmo.critical_density(zeta))*(1.00-bias)*(10**m500)*units.solMass)**(1.00/3.00)

rs = r500.to(units.Mpc).value/c500; del r500

ps  = pnorm*(3.00/8.00/np.pi)*(fb*mu/mue)
ps *= (((((2.5e2*Hz*Hz)**2.00)*((1.00-bias)*(10**(m500-15.00))*units.solMass)/(const.G**(0.50)))**(2.00/3.00)).to(units.keV/units.cm**3)).value
ps *= 1e10*ysznorm.value

mass = ((1.00-bias)*(10**(m500-14.00))*(cosmo.H0.value/70.00)/3.00)

limepsr = 1.00E-06
limdist = 5.00E+00*c500

# ----------------------------------------------------------------------
    
rline = np.logspace(-5,2,100)
mline = a10Profile(rline,0.00,ps,rs,0.00,alpha,beta,gamma,ap,c500,mass,limdist,limepsr,None)
mnorm = a10Profile(0.00,0.00,ps,rs,0.00,alpha,beta,gamma,ap,c500,mass,limdist,limepsr,None)

sintheta = np.sin(np.deg2rad(theta)) 
costheta = np.cos(np.deg2rad(theta)) 
cosctrDec = np.cos(np.deg2rad(mody))

xgrid, ygrid, rgrid = jitgrid(wgrid[0],wgrid[1],modx,mody,sintheta,costheta,np.cos(np.deg2rad(mody)),e,rs/dz.value)
mgrid = np.interp(rgrid,rline,mline)

rcube = np.sqrt(xgrid[:,:,None]**2+ygrid[:,:,None]**2+(wcube[2]*(np.pi/1.80e2)/np.sqrt(1.00-elos(e)**2)/(rs/dz.value))**2)
mcube = a10ProfileRadial(rcube,alpha,beta,gamma,ap,c500,mass)*2.00*ps
mcube = mcube*np.deg2rad(np.abs(wcube[2][0,0,1]-wcube[2][0,0,2]))*dz.value/np.sqrt(1.00-elos(e)**2)

ogrid = mgrid.copy()

# ----------------------------------------------------------------------

offx = 0#0.20/60.00
offy = 0.00

frnr = 1.00/60.00

beta =  4.00
mach =  2.00
phi0 = 60.00
theta = 80.00

sdeca = 2.00*(1.00-(1.00+(rcube*np.rad2deg(rs/dz.value)/frnr)**beta)**(np.log(0.50)/np.log(2.00)))

sintheta = np.sin(np.deg2rad(theta)) 
costheta = np.cos(np.deg2rad(theta)) 
xgrid, _, rgrid = jitgrid(wgrid[0],wgrid[1],modx,mody,sintheta,costheta,np.cos(np.deg2rad(mody)),e,rs/dz.value)

sazim = np.sqrt(rgrid[:,:,None]**2+(wcube[2]/rs)**2)
sazim = np.rad2deg(np.arccos(-xgrid[:,:,None]/sazim))
sazim = np.where(sazim<=phi0,(mach-1.00)*np.cos(0.50*np.pi*sazim/phi0)**2+1.00,1.00)
sazim = (sazim-1.00)*sdeca+1.00

sazim = 2.00*gadab*sazim*sazim-(gadab-1.00)
sazim = sazim/(gadab+1.00)

xazim, yazim, _ = jitgrid(wgrid[0],wgrid[1],modx,mody,sintheta,costheta,np.cos(np.deg2rad(mody)),0.00,np.deg2rad(frnr))
scube = np.sqrt((xazim[:,:,None]+offx/frnr)**2+(yazim[:,:,None]+offy/frnr)**2+(wcube[2]/np.sqrt(1.00-elos(e)**2)*(np.pi/1.80e2)/np.deg2rad(frnr))**2)

sazim[scube>=1.00] = np.nan

#plt.subplot(221); plt.imshow(sazim[:,:,100],origin='lower')
#plt.subplot(222); plt.imshow(sazim[:,:,0],origin='lower')
#plt.subplot(223); plt.imshow(mgrid,origin='lower')
#plt.subplot(224); plt.imshow(mgrid+np.nansum(mcube*sazim-mcube,axis=-1),origin='lower')
#plt.subplot(224); plt.imshow(np.nansum(mcube*sazim,axis=-1),origin='lower')
#plt.show(); plt.close()

ogrid = ogrid+np.nansum(mcube*sazim-mcube,axis=-1)

# ----------------------------------------------------------------------

offx = -0.10/60.00
offy = 0.00

frnr =  1.20/60.00

beta = 10.00
mach =  2.00
phi0 = 40.00
theta = -110.00

sdeca = 2.00*(1.00-(1.00+(rcube*np.rad2deg(rs/dz.value)/frnr)**beta)**(np.log(0.50)/np.log(2.00)))

sintheta = np.sin(np.deg2rad(theta)) 
costheta = np.cos(np.deg2rad(theta)) 
xgrid, _, rgrid = jitgrid(wgrid[0],wgrid[1],modx,mody,sintheta,costheta,np.cos(np.deg2rad(mody)),e,rs/dz.value)

sazim = np.sqrt(rgrid[:,:,None]**2+(wcube[2]/rs)**2)
sazim = np.rad2deg(np.arccos(-xgrid[:,:,None]/sazim))
sazim = np.where(sazim<=phi0,(mach-1.00)*np.cos(0.50*np.pi*sazim/phi0)**2+1.00,1.00)
sazim = (sazim-1.00)*sdeca+1.00

sazim = 2.00*gadab*sazim*sazim-(gadab-1.00)
sazim = sazim/(gadab+1.00)

xazim, yazim, _ = jitgrid(wgrid[0],wgrid[1],modx,mody,sintheta,costheta,np.cos(np.deg2rad(mody)),0.00,np.deg2rad(frnr))
scube = np.sqrt((xazim[:,:,None]+offx/frnr)**2+(yazim[:,:,None]+offy/frnr)**2+(wcube[2]/np.sqrt(1.00-elos(e)**2)*(np.pi/1.80e2)/np.deg2rad(frnr))**2)

sazim[scube>=1.00] = np.nan

ogrid = ogrid+np.nansum(mcube*sazim-mcube,axis=-1)

fits.writeto('test.fits',ogrid,overwrite=True)

