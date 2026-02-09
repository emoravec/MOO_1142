"""
Compares the X-ray and SZ profiles of MOO 1142.

Created in 2026-02
"""
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table,vstack
from astropy import units as u

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path
# -------------------------------------------------------------------------------------------- #
location = Path('/Users/emoravec/Documents/Research/merging_clusters/analysis/MOO_1142/profile_fitting/compare_xray_sz')

def SZ_gNFW(r, P_0, r_s, alpha, beta, gamma):
    denominator = ((r / r_s) ** gamma) * (1 + (r / r_s) ** alpha) ** ((beta - gamma) / alpha)
    return P_0 / denominator

def Xray_gNFW(r, P_0, r_s, alpha, beta, gamma):
    expression = ((r / r_s) ** -gamma) * (1 + (r / r_s) ** alpha) ** ((gamma - beta) / alpha)
    return P_0 * expression

r = np.arange(0, 30, 0.1)  # radius in arcsec
# -------------------------------------------------------------------------------------------- #
# Create SZ profile
sz_params = pd.read_csv(location / 'parameter_files/main_cluster_gNFW_sz.csv')
main_cluster_sz_profile = SZ_gNFW(r = r,
                                  P_0=sz_params['P_0'].values[0],
                                  r_s = sz_params['r_s_arcsec'].values[0], 
                                  alpha = sz_params['alpha'].values[0], 
                                  beta = sz_params['beta'].values[0], 
                                  gamma = sz_params['gamma'].values[0])
# Create X-ray profile
xray_params = pd.read_csv(location / 'parameter_files/main_cluster_gNFW_xray.csv')
main_cluster_xray_profile = Xray_gNFW(r = r,
                                  P_0=xray_params['P_0'].values[0],
                                  r_s = xray_params['r_s_arcsec'].values[0], 
                                  alpha = xray_params['alpha'].values[0], 
                                  beta = xray_params['beta'].values[0], 
                                  gamma = xray_params['gamma'].values[0])
# Plot main cluster profiles
plt.title('Main Cluster Pressure Profiles')
plt.plot(r, main_cluster_sz_profile,c='blue',label='SZ Profile')
plt.plot(r, main_cluster_xray_profile,c='red',label='X-ray Profile')
plt.xlabel('r (")')
plt.ylabel('P (keV cm$^{-3}$)')
plt.legend()
#plt.show()
plt.savefig(location / 'plots/MOO_1142_main_cluster_sz_xray_profiles.png',dpi=300)