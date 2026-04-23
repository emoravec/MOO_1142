"""
Compare the XMM number-density and SZ pressure profiles of MOO 1142.

The left-hand panel shows the XMM number-density profiles for the main cluster and west subcluster,
including the 90% confidence intervals from the saved XMM posterior.
The right-hand panel shows the SZ pressure profiles for the same two subclusters.

Created in 2026-04.
"""
import builtins
import sys
from pathlib import Path
from typing import Dict, List

from astropy import units as u
from astropy.cosmology import Planck18

import dill
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.special
import numpy.random._pickle as numpy_random_pickle
# -------------------------------------------------------------------------------------------- #
location = Path(__file__).resolve().parent
xmm_fit_dir = location.parent / "XMM/Barbavara_fit_2026-04"
xmm_state_path = xmm_fit_dir / "gnfw_circ+beta_circ_acfixed_NESTED/pmc_final.state"

PROFILE_SAMPLE_COUNT = 1000
R_SZ = np.arange(0.1, 100.0, 0.1)
R_XMM = np.arange(0.1, 100.0, 0.1)

ALPHA_MAIN = 2.26
GAMMA_MAIN = 0.465
F_S_MAIN = 1.632e4
F_S_SUB = 2.597e5
KT_MAIN_KEV = 6.76
KT_SUB_KEV = 7.13
XMM_PIXEL_SCALE_ARCSEC = 2.5
# -------------------------------------------------------------------------------------------- #
def SZ_gNFW(r: np.ndarray, p_0: float, r_s: float, alpha: float, beta: float, gamma: float) -> np.ndarray:
	denominator = ((r / r_s) ** gamma) * (1 + (r / r_s) ** alpha) ** ((beta - gamma) / alpha)
	return p_0 / denominator


def Xray_gNFW(r: np.ndarray, n_0: float, r_s: float, alpha: float, beta: float, gamma: float) -> np.ndarray:
	expression = ((r / r_s) ** -gamma) * (1 + (r / r_s) ** alpha) ** ((gamma - beta) / alpha)
	return n_0 * expression


def iso_beta(r: np.ndarray, c_0: float, r_c: float, beta: float) -> np.ndarray:
	expression = (1 + (r / r_c) ** 2) ** (-1.5 * beta)
	return c_0 * expression


def load_xmm_state(state_path: Path) -> Dict:
	"""Load the saved pocoMC state with compatibility shims for this environment."""
	sys.modules["__builtin__"] = builtins
	original_randomstate_ctor = numpy_random_pickle.__randomstate_ctor
	original_bitgen_ctor = numpy_random_pickle.__bit_generator_ctor

	def patched_randomstate_ctor(bit_generator_name: str = "MT19937", *args):
		return original_randomstate_ctor(bit_generator_name)

	def patched_bitgen_ctor(bit_generator_name: str = "MT19937", *args):
		return original_bitgen_ctor(bit_generator_name)

	numpy_random_pickle.__randomstate_ctor = patched_randomstate_ctor
	numpy_random_pickle.__bit_generator_ctor = patched_bitgen_ctor
	try:
		with state_path.open("rb") as file_handle:
			return dill.load(file_handle)
	finally:
		numpy_random_pickle.__randomstate_ctor = original_randomstate_ctor
		numpy_random_pickle.__bit_generator_ctor = original_bitgen_ctor


def get_xmm_posterior_samples(state_path: Path, sample_count: int) -> np.ndarray:
	"""Return the tail of the saved posterior samples used for the notebook profiles."""
	state = load_xmm_state(state_path)
	samples = state["particles"].get("x", flat=True)
	return samples[-sample_count:]


def percentile_triplet(profiles: List[np.ndarray], prefix: str) -> Dict[str, np.ndarray]:
	profile_p5, profile_p50, profile_p95 = np.percentile(profiles, [5, 50, 95], axis=0)
	return {
		f"{prefix}_p5": profile_p5,
		f"{prefix}_p50": profile_p50,
		f"{prefix}_p95": profile_p95,
	}


def compute_xmm_density_profiles(samples: np.ndarray, r_vals: np.ndarray) -> Dict[str, np.ndarray]:
	"""Build the XMM median and 90% interval density profiles for both subclusters."""
	n_profiles_main = []
	n_profiles_sub = []

	for sample in samples:
		(
			x_main_mcmc,
			y_main_mcmc,
			x_sub_mcmc,
			y_sub_mcmc,
			rs_main_mcmc,
			rc_sub_mcmc,
			beta_main_mcmc,
			beta_sub_mcmc,
			a_main_mcmc,
			a_sub_mcmc,
			a_bkg_mcmc,
		) = sample

		main_n_0 = np.sqrt(a_main_mcmc / F_S_MAIN)
		sub_n_0 = np.sqrt(
			a_sub_mcmc * scipy.special.gamma(3.0 * beta_sub_mcmc)
			/ (F_S_SUB * scipy.special.gamma(3.0 * beta_sub_mcmc - 0.5))
		)

		main_profile_density = Xray_gNFW(
			r=r_vals,
			n_0=main_n_0,
			r_s=rs_main_mcmc * XMM_PIXEL_SCALE_ARCSEC,
			alpha=ALPHA_MAIN,
			beta=beta_main_mcmc / 2.0,
			gamma=GAMMA_MAIN,
		)

		sub_profile_density = iso_beta(
			r=r_vals,
			c_0=sub_n_0,
			r_c=rc_sub_mcmc * XMM_PIXEL_SCALE_ARCSEC,
			beta=beta_sub_mcmc,
		)

		n_profiles_main.append(main_profile_density)
		n_profiles_sub.append(sub_profile_density)

	profile_summary = {}
	profile_summary.update(percentile_triplet(n_profiles_main, "n_main"))
	profile_summary.update(percentile_triplet(n_profiles_sub, "n_sub"))
	return profile_summary


def format_radius_axes(axis: plt.Axes, kpc_per_arcsec: float) -> None:
	"""Apply the shared arcsec and kpc radius axis formatting."""
	axis.set_xlim(10, 90)
	axis.set_xlabel('r (")')
	axis.set_xticks([10, 20, 30, 40, 50, 60, 70, 80, 90])
	axis.set_xticklabels(["10", "20", "30", "40", "50", "60", "70", "80", "90"])

	axis_top = axis.secondary_xaxis("top", functions=(lambda x: x * kpc_per_arcsec, lambda x: x / kpc_per_arcsec))
	axis_top.set_xlabel("r (kpc)")
	axis_top.set_xticks([100, 200, 300, 400, 500, 600, 700])
	axis_top.set_xticklabels(["100", "200", "300", "400", "500", "600", "700"])
# -------------------------------------------------------------------------------------------- #
def main() -> None:
	main_sz_params = pd.read_csv(location / "parameter_files/main_cluster_gNFW_sz.csv")
	main_cluster_sz_profile = SZ_gNFW(
		r=R_SZ,
		p_0=main_sz_params["P_0"].values[0],
		r_s=main_sz_params["r_s_arcsec"].values[0],
		alpha=main_sz_params["alpha"].values[0],
		beta=main_sz_params["beta"].values[0],
		gamma=main_sz_params["gamma"].values[0],
	)

	sub_sz_params = pd.read_csv(location / "parameter_files/subcluster_sph_isobeta_sz.csv")
	subcluster_sz_profile = iso_beta(
		r=R_SZ,
		c_0=sub_sz_params["P_0"].values[0],
		r_c=sub_sz_params["r_c_arcsec"].values[0],
		beta=sub_sz_params["beta"].values[0],
	)

	xmm_samples = get_xmm_posterior_samples(xmm_state_path, PROFILE_SAMPLE_COUNT)
	xmm_profiles = compute_xmm_density_profiles(xmm_samples, R_XMM)

	fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

	z = 1.189
	kpc_per_arcsec = Planck18.kpc_proper_per_arcmin(z).to(u.kpc / u.arcsec).value

	ax1.set_title("Number Density Profiles")
	ax1.loglog(R_XMM, xmm_profiles["n_main_p50"], c="navy", label="Main")
	ax1.fill_between(R_XMM, xmm_profiles["n_main_p5"], xmm_profiles["n_main_p95"], color="navy", alpha=0.25)
	ax1.loglog(R_XMM, xmm_profiles["n_sub_p50"], c="darkorange", label="West")
	ax1.fill_between(R_XMM, xmm_profiles["n_sub_p5"], xmm_profiles["n_sub_p95"], color="darkorange", alpha=0.25)
	ax1.set_ylabel("n (cm$^{-3}$)")
	ax1.legend()
	format_radius_axes(ax1, kpc_per_arcsec)

	ax2.set_title("Pressure Profiles")
	ax2.loglog(R_SZ, main_cluster_sz_profile, c="navy", label="Main")
	ax2.loglog(R_SZ, subcluster_sz_profile, c="darkorange", label="West")
	ax2.set_ylabel("P (keV cm$^{-3}$)")
	ax2.legend()
	format_radius_axes(ax2, kpc_per_arcsec)

	plt.tight_layout()
	plt.savefig(location / "plots/MOO_1142_density_pressure_profiles_90CL.png", dpi=300)
	# plt.show()


if __name__ == "__main__":
	main()