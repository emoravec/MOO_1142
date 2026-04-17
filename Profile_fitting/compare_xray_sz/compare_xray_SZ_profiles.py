"""
Compares the X-ray and SZ profiles of MOO 1142.

Created in 2026-02
"""
from astropy.coordinates import SkyCoord
from astropy.cosmology import Planck18
from astropy.io import fits
from astropy.table import Table,vstack
from astropy import units as u

import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import numpy as np
import pandas as pd
from pathlib import Path
# -------------------------------------------------------------------------------------------- #
location = Path(__file__).resolve().parent

def SZ_gNFW(r, P_0, r_s, alpha, beta, gamma):
    denominator = ((r / r_s) ** gamma) * (1 + (r / r_s) ** alpha) ** ((beta - gamma) / alpha)
    return P_0 / denominator

def Xray_gNFW(r, P_0, r_s, alpha, beta, gamma):
    expression = ((r / r_s) ** -gamma) * (1 + (r / r_s) ** alpha) ** ((gamma - beta) / alpha)
    return P_0 * expression

def iso_beta(r, P_0, r_c, beta):
    expression = (1 + (r / r_c) ** 2) ** (-1.5*beta)
    return P_0 * expression


r = np.arange(0, 100, 0.1)  # radius in arcsec
# -------------------------------------------------------------------------------------------- #
### MAIN CLUSTER ###
# Create SZ profile
main_sz_params = pd.read_csv(location / 'parameter_files/main_cluster_gNFW_sz.csv')
main_cluster_sz_profile = SZ_gNFW(r = r,
                                  P_0= main_sz_params['P_0'].values[0],
                                  r_s = main_sz_params['r_s_arcsec'].values[0], 
                                  alpha = main_sz_params['alpha'].values[0], 
                                  beta = main_sz_params['beta'].values[0], 
                                  gamma = main_sz_params['gamma'].values[0])
# Create X-ray profile
main_xray_params = pd.read_csv(location / 'parameter_files/main_cluster_gNFW_xray.csv')
main_cluster_xray_profile = Xray_gNFW(r = r,
                                  P_0= main_xray_params['P_0'].values[0],
                                  r_s = main_xray_params['r_s_arcsec'].values[0], 
                                  alpha = main_xray_params['alpha'].values[0], 
                                  beta = main_xray_params['beta'].values[0], 
                                  gamma = main_xray_params['gamma'].values[0])

### Subcluster ###
# Create SZ profile
sub_sz_params = pd.read_csv(location / 'parameter_files/subcluster_sph_isobeta_sz.csv')
subcluster_sz_profile = iso_beta(r = r,
                                 P_0= sub_sz_params['P_0'].values[0],
                                 r_c = sub_sz_params['r_c_arcsec'].values[0], 
                                 beta = sub_sz_params['beta'].values[0])
# Create X-ray profile
sub_xray_params = pd.read_csv(location / 'parameter_files/subcluster_sph_isobeta_xray.csv')
subcluster_xray_profile = iso_beta(r = r,
                                 P_0= sub_xray_params['P_0'].values[0],
                                 r_c = sub_xray_params['r_c_arcsec'].values[0], 
                                 beta = sub_xray_params['beta'].values[0])

# -------------------------------------------------------------------------------------------- #
# Create figure with 2 horizontal subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Cosmology for kpc conversion
z = 1.189  # MOO 1142 biweight redshift
kpc_per_arcsec = Planck18.kpc_proper_per_arcmin(z).to(u.kpc / u.arcsec).value

# Plot Main Cluster
ax1.set_title('Main Cluster Pressure Profile')
ax1.loglog(r, main_cluster_sz_profile, c='blue', label='M2')
ax1.loglog(r, main_cluster_xray_profile, c='red', label='XMM')
ax1.set_xlim(10, 9e1)
ax1.set_xlabel('r (")')
ax1.set_ylabel('P (keV cm$^{-3}$)')
ax1.legend()

# Format x-axis labels as regular numbers
ax1.set_xticks([10, 20, 30, 40, 50, 60, 70, 80, 90])
ax1.set_xticklabels(['10', '20', '30', '40', '50', '60', '70', '80', '90'])

# Add secondary x-axis for main cluster
ax1_top = ax1.secondary_xaxis('top', functions=(lambda x: x * kpc_per_arcsec, 
                                                 lambda x: x / kpc_per_arcsec))
ax1_top.set_xlabel('r (kpc)')
ax1_top.set_xticks([100, 200, 300, 400, 500, 600, 700])
ax1_top.set_xticklabels(['100', '200', '300', '400', '500', '600', '700'])

# Plot Subcluster
ax2.set_title('Sub-cluster Pressure Profile')
ax2.loglog(r, subcluster_sz_profile, c='blue', label='M2')
ax2.loglog(r, subcluster_xray_profile, c='red', label='XMM')
ax2.set_xlim(10, 9e1)
ax2.set_xlabel('r (")')
ax2.set_ylabel('P (keV cm$^{-3}$)')
ax2.legend()

# Format x-axis labels as regular numbers
ax2.set_xticks([10, 20, 30, 40, 50, 60, 70, 80, 90])
ax2.set_xticklabels(['10', '20', '30', '40', '50', '60', '70', '80', '90'])

# Add secondary x-axis for subcluster
ax2_top = ax2.secondary_xaxis('top', functions=(lambda x: x * kpc_per_arcsec, 
                                                 lambda x: x / kpc_per_arcsec))
ax2_top.set_xlabel('r (kpc)')
ax2_top.set_xticks([100, 200, 300, 400, 500, 600, 700])
ax2_top.set_xticklabels(['100', '200', '300', '400', '500', '600', '700'])

plt.tight_layout()
plt.savefig(location / 'plots/MOO_1142_main+west_sz_xray_profiles.png', dpi=300)
#plt.show()