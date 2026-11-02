"""
Combines the Xray (density - n_e) and SZ (pressure - P_e) profiles to calculate the temperature and entropy profiles of MOO 1142.

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

try:
    from .profile_models import gNFW, iso_beta
except ImportError:
    from profile_models import gNFW, iso_beta

# -------------------------------------------------------------------------------------------- #
location = Path(__file__).resolve().parent

def kb_temperature_profile(pressure_profile, density_profile):
    """Calculate the temperature profile from the pressure and density profiles."""
    return pressure_profile / density_profile

def entropy_profile(pressure_profile, density_profile):
    """Calculate the entropy profile from the pressure and density profiles."""
    return pressure_profile / (density_profile ** (5/3))

r = np.arange(0, 100, 0.1)  # radius in arcsec
main_avg_spec_temp_keV = 6.76
pos_avg_spec_temp_keV = main_avg_spec_temp_keV + 0.75
neg_avg_spec_temp_keV = main_avg_spec_temp_keV - 0.54

binned_temp_radii = np.array([6.76, 15.35, 26.69, 53.12])
binned_temp_radii_bins = np.array([0, 11, 21, 36, 94])
binned_temp_kt_deproj = np.array([5.75, 8.37, 6.98, 5.25])
binned_temp_kt_lower = np.array([-0.96, -1.45, -1.33, -0.85])
binned_temp_kt_upper = np.array([1.12, 1.90, 1.93, 1.35])
# -------------------------------------------------------------------------------------------- #
### MAIN CLUSTER ###
# Create SZ pressure profile
main_sz_params = pd.read_csv(location / 'parameter_files/main_cluster_gNFW_sz.csv')
main_cluster_sz_profile = gNFW(r = r,
                               p_0= main_sz_params['P_0'].values[0],
                               r_s = main_sz_params['r_s_arcsec'].values[0], 
                               alpha = main_sz_params['alpha'].values[0], 
                               beta = main_sz_params['beta'].values[0], 
                               gamma = main_sz_params['gamma'].values[0])
# Create X-ray density profile
main_xray_params = pd.read_csv(location / 'parameter_files/main_cluster_gNFW_xray.csv')
main_cluster_xray_profile = gNFW(r = r,
                                 p_0= main_xray_params['n_s'].values[0],
                                 r_s = main_xray_params['r_s_arcsec'].values[0], 
                                 alpha = main_xray_params['alpha'].values[0], 
                                 beta = main_xray_params['beta'].values[0], 
                                 gamma = main_xray_params['gamma'].values[0])
# Thermodynamic profiles for main cluster
main_cluster_temperature_profile = kb_temperature_profile(main_cluster_sz_profile, main_cluster_xray_profile)
main_cluster_entropy_profile = entropy_profile(main_cluster_sz_profile, main_cluster_xray_profile)

### Subcluster ###
# Create SZ pressure profile
sub_sz_params = pd.read_csv(location / 'parameter_files/subcluster_sph_isobeta_sz.csv')
subcluster_sz_profile = iso_beta(r = r,
                                                                 p_0= sub_sz_params['P_0'].values[0],
                                 r_c = sub_sz_params['r_c_arcsec'].values[0], 
                                 beta = sub_sz_params['beta'].values[0])
# Create X-ray density profile
sub_xray_params = pd.read_csv(location / 'parameter_files/subcluster_sph_isobeta_xray.csv')
subcluster_xray_profile = iso_beta(r = r,
                                                                     p_0= sub_xray_params['n_0'].values[0],
                                                                     r_c = sub_xray_params['r_c_arcsec'].values[0], 
                                                                     beta = sub_xray_params['beta'].values[0])

# Thermodynamic profiles for subcluster
subcluster_temperature_profile = kb_temperature_profile(subcluster_sz_profile, subcluster_xray_profile)
subcluster_entropy_profile = entropy_profile(subcluster_sz_profile, subcluster_xray_profile)

# -------------------------------------------------------------------------------------------- #
# Create figure with 2 horizontal subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Cosmology for kpc conversion
z = 1.189  # MOO 1142 biweight redshift
kpc_per_arcsec = Planck18.kpc_proper_per_arcmin(z).to(u.kpc / u.arcsec).value

# Plot kB Temperature profile
ax1.loglog(r, main_cluster_temperature_profile, c='k')
# add average spectral temperature
ax1.axhline(y=main_avg_spec_temp_keV, color='red', linestyle='--')
ax1.axhspan(neg_avg_spec_temp_keV, pos_avg_spec_temp_keV, color='red', alpha=0.3)

# add binned spectral temperatures (x bars span the bin edges)
ax1.errorbar(binned_temp_radii, binned_temp_kt_deproj,
             xerr=[binned_temp_radii - binned_temp_radii_bins[:-1],
                   binned_temp_radii_bins[1:] - binned_temp_radii],
             yerr=[np.abs(binned_temp_kt_lower), binned_temp_kt_upper],
             fmt='o', color='blue', ecolor='blue', capsize=3, zorder=3)

# Format x-axis
ax1.set_xlim(5, 9e1)
ax1.set_xlabel('r (")')
ax1.set_ylim(10**0, 10**1)
ax1.set_yticks(np.arange(1, 11, 1))
ax1.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f'{value:g}'))
ax1.set_ylabel('k$_{B}$ T$_{e}$ (keV)')
ax1.tick_params(axis='y', which='both', right=True, labelright=False)

# misc plot properties
ax1.set_title('Main Subcluster Temperature Profile')
ax1.legend().remove()

# Format x-axis labels as regular numbers
ax1.set_xticks([5, 10, 20, 30, 40, 50, 60, 70, 80, 90])
ax1.set_xticklabels(['5', '10', '20', '30', '40', '50', '60', '70', '80', '90'])

# Add secondary x-axis for main cluster
ax1_top = ax1.secondary_xaxis('top', functions=(lambda x: x * kpc_per_arcsec, 
                                                 lambda x: x / kpc_per_arcsec))
ax1_top.set_xlabel('r (kpc)')
ax1_top.set_xticks([100, 200, 300, 400, 500, 600, 700])
ax1_top.set_xticklabels(['100', '200', '300', '400', '500', '600', '700'])

# Plot entropy profile
ax2.loglog(r, main_cluster_entropy_profile, c='k')

# Format x-axis
ax2.set_xlim(5, 9e1)
ax2.set_xlabel('r (")')
ax2.set_ylim(10**1, 10**3)
ax2.set_ylabel('K (keV cm$^{-2}$)')
ax2.tick_params(axis='y', which='both', right=True, labelright=False)

# misc plot properties
ax2.legend().remove()
ax2.set_title('Main Subcluster Entropy Profile')

# Format x-axis labels as regular numbers
ax2.set_xticks([5, 10, 20, 30, 40, 50, 60, 70, 80, 90])
ax2.set_xticklabels(['5', '10', '20', '30', '40', '50', '60', '70', '80', '90'])

# Add secondary x-axis for subcluster
ax2_top = ax2.secondary_xaxis('top', functions=(lambda x: x * kpc_per_arcsec, 
                                                 lambda x: x / kpc_per_arcsec))
ax2_top.set_xlabel('r (kpc)')
ax2_top.set_xticks([100, 200, 300, 400, 500, 600, 700])
ax2_top.set_xticklabels(['100', '200', '300', '400', '500', '600', '700'])

plt.tight_layout()
plt.savefig(location / 'plots/MOO_1142_temp_entropy_profiles_main.png', dpi=300)
#plt.show()