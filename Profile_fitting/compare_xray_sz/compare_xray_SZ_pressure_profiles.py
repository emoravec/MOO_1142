"""
Compares the X-ray and SZ profiles of MOO 1142.

Created in 2026-02
"""
from astropy.cosmology import Planck18
from astropy import units as u

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path

try:
    from .profile_models import gNFW, iso_beta
except ImportError:
    from profile_models import gNFW, iso_beta

# -------------------------------------------------------------------------------------------- #
location = Path(__file__).resolve().parent

def _get_param_value(params, key):
    value = params[key]
    if hasattr(value, 'values'):
        return value.values[0]
    return value[0]

def load_profile_parameters(base_path=location):
    return {
        'main_sz': pd.read_csv(base_path / 'parameter_files/main_cluster_gNFW_sz.csv'),
        'main_xray': pd.read_csv(base_path / 'parameter_files/main_cluster_gNFW_xray.csv'),
        'sub_sz': pd.read_csv(base_path / 'parameter_files/subcluster_sph_isobeta_sz.csv'),
        'sub_xray': pd.read_csv(base_path / 'parameter_files/subcluster_sph_isobeta_xray.csv'),
    }

def build_pressure_profiles(r, params):
    main_cluster_sz_profile = gNFW(
        r=r,
        p_0=_get_param_value(params['main_sz'], 'P_0'),
        r_s=_get_param_value(params['main_sz'], 'r_s_arcsec'),
        alpha=_get_param_value(params['main_sz'], 'alpha'),
        beta=_get_param_value(params['main_sz'], 'beta'),
        gamma=_get_param_value(params['main_sz'], 'gamma'),
    )
    main_cluster_xray_profile = gNFW(
        r=r,
        p_0=_get_param_value(params['main_xray'], 'P_0'),
        r_s=_get_param_value(params['main_xray'], 'r_s_arcsec'),
        alpha=_get_param_value(params['main_xray'], 'alpha'),
        beta=_get_param_value(params['main_xray'], 'beta'),
        gamma=_get_param_value(params['main_xray'], 'gamma'),
    )
    subcluster_sz_profile = iso_beta(
        r=r,
        p_0=_get_param_value(params['sub_sz'], 'P_0'),
        r_c=_get_param_value(params['sub_sz'], 'r_c_arcsec'),
        beta=_get_param_value(params['sub_sz'], 'beta'),
    )
    subcluster_xray_profile = iso_beta(
        r=r,
        p_0=_get_param_value(params['sub_xray'], 'P_0'),
        r_c=_get_param_value(params['sub_xray'], 'r_c_arcsec'),
        beta=_get_param_value(params['sub_xray'], 'beta'),
    )

    return {
        'main_sz': main_cluster_sz_profile,
        'main_xray': main_cluster_xray_profile,
        'sub_sz': subcluster_sz_profile,
        'sub_xray': subcluster_xray_profile,
    }

def make_pressure_profile_figure(r, profiles, z=1.189):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    kpc_per_arcsec = Planck18.kpc_proper_per_arcmin(z).to(u.kpc / u.arcsec).value

    ax1.set_title('Main Cluster Pressure Profile')
    ax1.loglog(r, profiles['main_sz'], c='blue', label='M2')
    ax1.loglog(r, profiles['main_xray'], c='red', label='XMM')
    ax1.set_xlim(10, 9e1)
    ax1.set_xlabel('r (")')
    ax1.set_ylabel('P (keV cm$^{-3}$)')
    ax1.legend()
    ax1.set_xticks([10, 20, 30, 40, 50, 60, 70, 80, 90])
    ax1.set_xticklabels(['10', '20', '30', '40', '50', '60', '70', '80', '90'])

    ax1_top = ax1.secondary_xaxis('top', functions=(lambda x: x * kpc_per_arcsec,
                                                    lambda x: x / kpc_per_arcsec))
    ax1_top.set_xlabel('r (kpc)')
    ax1_top.set_xticks([100, 200, 300, 400, 500, 600, 700])
    ax1_top.set_xticklabels(['100', '200', '300', '400', '500', '600', '700'])

    ax2.set_title('Sub-cluster Pressure Profile')
    ax2.loglog(r, profiles['sub_sz'], c='blue', label='M2')
    ax2.loglog(r, profiles['sub_xray'], c='red', label='XMM')
    ax2.set_xlim(10, 9e1)
    ax2.set_xlabel('r (")')
    ax2.set_ylabel('P (keV cm$^{-3}$)')
    ax2.legend()
    ax2.set_xticks([10, 20, 30, 40, 50, 60, 70, 80, 90])
    ax2.set_xticklabels(['10', '20', '30', '40', '50', '60', '70', '80', '90'])

    ax2_top = ax2.secondary_xaxis('top', functions=(lambda x: x * kpc_per_arcsec,
                                                    lambda x: x / kpc_per_arcsec))
    ax2_top.set_xlabel('r (kpc)')
    ax2_top.set_xticks([100, 200, 300, 400, 500, 600, 700])
    ax2_top.set_xticklabels(['100', '200', '300', '400', '500', '600', '700'])

    fig.tight_layout()
    return fig

def main():
    r = np.arange(0, 100, 0.1)
    params = load_profile_parameters()
    profiles = build_pressure_profiles(r, params)
    fig = make_pressure_profile_figure(r, profiles)
    fig.savefig(location / 'plots/pressure/MOO_1142_main+west_pressure_profiles.png', dpi=300)

if __name__ == '__main__':
    main()