"""
Compare the SZ and X-ray pressure profiles of MOO 1142 and show the XMM 90% confidence intervals.

This combines the base two-panel comparison plot from compare_xray_SZ_profiles.py
with the posterior-profile construction used in radial_profile_and_uncert_n.ipynb.

Created in 2026-04.
"""
from __future__ import annotations

# Provenance in this file is recorded at the function or block level rather than per line.
# Source names used below:
# - compare_xray_SZ_pressure_profiles.py: base two-panel SZ/X-ray comparison figure.
# - radial_profile_and_uncert_n.ipynb: XMM posterior sampling and pressure-profile construction.
# - standalone 90% CL script: glue added here to load saved state and combine both workflows.

# Standard-library imports.
import builtins
import sys
from pathlib import Path

# Shared scientific Python stack used by the comparison script and notebook workflow.
from astropy import units as u
from astropy.cosmology import Planck18

import dill
import matplotlib.pyplot as plt
import numpy as np
import numpy.random._pickle as numpy_random_pickle
import pandas as pd
import scipy.special

try:
	from .profile_models import gNFW, iso_beta
except ImportError:
	from profile_models import gNFW, iso_beta

# -------------------------------------------------------------------------------------------- #
location = Path(__file__).resolve().parent

# Provenance: adapted from radial_profile_and_uncert_n.ipynb, saved-state posterior workflow.
# Intent: keep the XMM fit location and posterior-sampling constants together so the standalone
# script can reproduce the notebook pressure-profile construction outside the notebook session.
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
def load_xmm_state(state_path: Path) -> dict:
	# Provenance: standalone 90% CL script.
	# Intent: load the saved pocoMC state outside the original notebook session by patching the
	# legacy objects expected by the serialized RNG payload.
	"""Load the saved pocoMC state with compatibility shims for this environment."""
	sys.modules["__builtin__"] = builtins
	original_randomstate_ctor = numpy_random_pickle.__randomstate_ctor
	original_bitgen_ctor = numpy_random_pickle.__bit_generator_ctor

	def patched_randomstate_ctor(bit_generator_name: str = "MT19937", *args):
		# Ignore extra pickle args so older numpy random-state payloads still load.
		return original_randomstate_ctor(bit_generator_name)

	def patched_bitgen_ctor(bit_generator_name: str = "MT19937", *args):
		# Ignore extra pickle args for the bit-generator payload as well.
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
	# Provenance: adapted from radial_profile_and_uncert_n.ipynb posterior access.
	# Intent: reconstruct the same posterior tail used for the notebook pressure-profile summaries,
	# but from the serialized state instead of a live sampler object.
	"""Return the tail of the saved posterior samples used for the notebook profiles."""
	state = load_xmm_state(state_path)
	samples = state["particles"].get("x", flat=True)
	return samples[-sample_count:]

def compute_xmm_pressure_profiles(samples: np.ndarray, r_vals: np.ndarray) -> dict[str, np.ndarray]:
	# Provenance: adapted from radial_profile_and_uncert_n.ipynb, posterior-profile construction.
	# Intent: turn each posterior sample into a main-cluster and west-subcluster pressure profile,
	# then summarize the ensemble with the same 5/50/95 percentile curves used in the notebook.
	"""Build the XMM median and 90% interval pressure profiles for both subclusters."""
	p_profiles_main = []
	p_profiles_sub = []

	for sample in samples:
		# Match the posterior parameter ordering used by the saved notebook fit.
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

		# Keep the notebook normalization and shape conventions, with radii converted to arcsec.
		main_profile = gNFW(
			r=r_vals,
			p_0=np.sqrt(a_main_mcmc / F_S_MAIN) * KT_MAIN_KEV,
			r_s=rs_main_mcmc * XMM_PIXEL_SCALE_ARCSEC,
			alpha=ALPHA_MAIN,
			beta=beta_main_mcmc / 2.0,
			gamma=GAMMA_MAIN,
		)
		sub_profile = iso_beta(
			r=r_vals,
			p_0=np.sqrt(
				a_sub_mcmc * scipy.special.gamma(3.0 * beta_sub_mcmc)
				/ (F_S_SUB * scipy.special.gamma(3.0 * beta_sub_mcmc - 0.5))
			)
			* KT_SUB_KEV,
			r_c=rc_sub_mcmc * XMM_PIXEL_SCALE_ARCSEC,
			beta=beta_sub_mcmc,
		)

		p_profiles_main.append(main_profile)
		p_profiles_sub.append(sub_profile)

	p_p5_main, p_p50_main, p_p95_main = np.percentile(p_profiles_main, [5, 50, 95], axis=0)
	p_p5_sub, p_p50_sub, p_p95_sub = np.percentile(p_profiles_sub, [5, 50, 95], axis=0)
	return {
		"P_p5_main": p_p5_main,
		"P_p50_main": p_p50_main,
		"P_p95_main": p_p95_main,
		"P_p5_sub": p_p5_sub,
		"P_p50_sub": p_p50_sub,
		"P_p95_sub": p_p95_sub,
	}


def format_radius_axes(
	axis: plt.Axes,
	kpc_per_arcsec: float,
	*,
	hide_leftmost_label: bool = False,
	hide_rightmost_label: bool = False,
) -> None:
	# Provenance: adapted from compare_xray_SZ_pressure_profiles.py, shared radius-axis formatting.
	# Intent: keep both panels on the same arcsec and kpc axis scheme while allowing the seam labels
	# to be suppressed where the panels touch.
	"""Apply the shared arcsec and kpc radius axis formatting."""
	bottom_tick_labels = ["10", "20", "30", "40", "50", "60", "70", "80", "90"]
	if hide_leftmost_label:
		bottom_tick_labels[0] = ""
	if hide_rightmost_label:
		bottom_tick_labels[-1] = ""

	axis.set_xlim(10, 90)
	axis.set_xlabel('r (")')
	axis.set_xticks([10, 20, 30, 40, 50, 60, 70, 80, 90])
	axis.set_xticklabels(bottom_tick_labels)

	axis_top = axis.secondary_xaxis("top", functions=(lambda x: x * kpc_per_arcsec, lambda x: x / kpc_per_arcsec))
	axis_top.set_xlabel("r (kpc)")
	axis_top.set_xticks([100, 200, 300, 400, 500, 600, 700])
	axis_top.set_xticklabels(["100", "200", "300", "400", "500", "600", "700"])


def format_shared_pressure_axes(left_axis: plt.Axes, right_axis: plt.Axes) -> None:
	# Provenance: standalone 90% CL script, shared two-panel layout.
	# Intent: place one effective pressure scale between the panels while preserving readable ticks on
	# the outer left and right edges of the combined figure.
	"""Apply the shared log-pressure axis styling for the two-panel figure."""
	left_axis.set_ylim(1e-7, 1e0)
	left_axis.set_ylabel("P (keV cm$^{-3}$)")
	left_axis.set_zorder(2)
	left_axis.patch.set_visible(False)
	left_axis.tick_params(axis="y", which="both", left=True, labelleft=True, right=True, labelright=False)

	right_axis.yaxis.tick_right()
	right_axis.tick_params(axis="y", which="both", left=False, labelleft=False, right=True, labelright=True)
	right_axis.spines["left"].set_visible(False)
# -------------------------------------------------------------------------------------------- #
def main() -> None:
	# Provenance: adapted from compare_xray_SZ_pressure_profiles.py, main two-panel plotting setup.
	# Intent: keep the SZ best-fit curves and overall figure structure from the original comparison
	# script, then replace the XMM single best-fit curves with posterior median and 90% interval bands.
	main_sz_params = pd.read_csv(location / "parameter_files/main_cluster_gNFW_sz.csv")
	main_cluster_sz_profile = gNFW(
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
		p_0=sub_sz_params["P_0"].values[0],
		r_c=sub_sz_params["r_c_arcsec"].values[0],
		beta=sub_sz_params["beta"].values[0],
	)

	xmm_samples = get_xmm_posterior_samples(xmm_state_path, PROFILE_SAMPLE_COUNT)
	xmm_profiles = compute_xmm_pressure_profiles(xmm_samples, R_XMM)

	fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), sharey=True, gridspec_kw={"wspace": 0.0})

	z = 1.189
	kpc_per_arcsec = Planck18.kpc_proper_per_arcmin(z).to(u.kpc / u.arcsec).value

	ax1.set_title("Main Cluster Pressure Profile")
	ax1.loglog(R_SZ, main_cluster_sz_profile, c="blue", label="M2")
	ax1.loglog(R_XMM, xmm_profiles["P_p50_main"], c="red", label="XMM")
	ax1.fill_between(R_XMM, xmm_profiles["P_p5_main"], xmm_profiles["P_p95_main"], color="red", alpha=0.3)
	ax1.legend()
	format_radius_axes(ax1, kpc_per_arcsec, hide_rightmost_label=True)

	ax2.set_title("West Subcluster Pressure Profile")
	ax2.loglog(R_SZ, subcluster_sz_profile, c="blue", label="M2")
	ax2.loglog(R_XMM, xmm_profiles["P_p50_sub"], c="red", label="XMM")
	ax2.fill_between(R_XMM, xmm_profiles["P_p5_sub"], xmm_profiles["P_p95_sub"], color="red", alpha=0.3)
	ax2.legend()
	format_radius_axes(ax2, kpc_per_arcsec)
	format_shared_pressure_axes(ax1, ax2)

	fig.subplots_adjust(wspace=0.0)
	plt.savefig(location / "plots/pressure/MOO_1142_main+west_sz_xray_profiles_90CL_shared.png", dpi=300)
	#plt.show()

if __name__ == "__main__":
	main()