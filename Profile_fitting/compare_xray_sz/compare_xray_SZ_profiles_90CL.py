"""
Compare the SZ and X-ray pressure profiles of MOO 1142 and show the XMM 90% confidence intervals.

This combines the base two-panel comparison plot from compare_xray_SZ_profiles.py
with the posterior-profile construction used in radial_profile_and_uncert_n.ipynb.

Created in 2026-04.
"""
from __future__ import annotations  # A: added here for modern type hints; not taken from either source file.

# Provenance shorthand used below:
# C = adapted from compare_xray_SZ_profiles.py written by Emily Moravec.
# N = adapted from radial_profile_and_uncert_n.ipynb writted by Eleonora Barbavara.
# A = added here while turning the comparison into a standalone 90% CL script written with Co-pilot.

import builtins  # A: added here so dill can load the saved pocoMC state outside the original notebook session.
import sys  # A: added here for the __builtin__ compatibility shim used during state loading.
from pathlib import Path  # C line 16: same Path usage as the original comparison script.

from astropy import units as u  # C line 10 and N line 17: both source files use astropy units for kpc/arcsec conversion.
from astropy.cosmology import Planck18  # C line 7: same cosmology object as the original comparison script.

import dill  # A: added here to deserialize pmc_final.state directly.
import matplotlib.pyplot as plt  # C line 12 and N line 12: same plotting library used in both source files.
import numpy as np  # C line 14 and N line 11: same array library used in both source files.
import pandas as pd  # C line 15: same CSV-reading library as the original comparison script.
import scipy.special  # N line 18: the notebook imports scipy there, and this script uses scipy.special for the west-subcluster isobeta normalization.
import numpy.random._pickle as numpy_random_pickle  # A: added here to patch numpy RNG pickles when loading the saved state.
# -------------------------------------------------------------------------------------------- #
location = Path(__file__).resolve().parent 
xmm_fit_dir = location.parent / "XMM/Barbavara_fit_2026-04"  # A: added here to point from the comparison-script folder to the XMM fit folder.
xmm_state_path = xmm_fit_dir / "gnfw_circ+beta_circ_acfixed_NESTED/pmc_final.state"  # A: added here to load the saved posterior instead of reading XMM CSV best fits.

PROFILE_SAMPLE_COUNT = 1000  # N lines 264 and 283: notebook builds both pressure-profile ensembles from 1000 posterior samples.
R_SZ = np.arange(0.1, 100.0, 0.1)  # C line 33 with A tweak: same 0-100 arcsec plotting idea as the original script, but starts at 0.1 to avoid log(0).
R_XMM = np.arange(0.1, 100.0, 0.1)  # A informed by C line 33 and N line 258: use an arcsec grid and extend it to the same plotting range as the original script.

ALPHA_MAIN = 2.26  # N line 180: same fixed main-cluster alpha used in the XMM notebook fit setup.
GAMMA_MAIN = 0.465  # N line 181: same fixed main-cluster gamma used in the XMM notebook fit setup.
F_S_MAIN = 1.632e4  # N line 260: same main-cluster normalization factor used when turning A_main into density/pressure.
F_S_SUB = 2.597e5  # N line 279: same west-subcluster normalization factor used in the notebook.
KT_MAIN_KEV = 6.76  # N line 261: same main-cluster temperature used to convert density to pressure.
KT_SUB_KEV = 7.13  # N line 280: same west-subcluster temperature used to convert density to pressure.
XMM_PIXEL_SCALE_ARCSEC = 2.5  # A from config_MOO1142.yaml line 29: added here to convert posterior radii from pixels to arcsec.
# -------------------------------------------------------------------------------------------- #
def SZ_gNFW(r: np.ndarray, p_0: float, r_s: float, alpha: float, beta: float, gamma: float) -> np.ndarray:
	# C lines 20-22: same SZ gNFW expression as the original comparison script, with only naming/typing cleaned up here.
	denominator = ((r / r_s) ** gamma) * (1 + (r / r_s) ** alpha) ** ((beta - gamma) / alpha)
	return p_0 / denominator

def Xray_gNFW(r: np.ndarray, p_0: float, r_s: float, alpha: float, beta: float, gamma: float) -> np.ndarray:
	# C lines 24-26 and N line 252: same gNFW profile form used by both the original script and the notebook.
	expression = ((r / r_s) ** -gamma) * (1 + (r / r_s) ** alpha) ** ((gamma - beta) / alpha)
	return p_0 * expression

def iso_beta(r: np.ndarray, p_0: float, r_c: float, beta: float) -> np.ndarray:
	# C lines 28-30 and N line 255: same isobeta profile form used by both the original script and the notebook.
	expression = (1 + (r / r_c) ** 2) ** (-1.5 * beta)
	return p_0 * expression

def load_xmm_state(state_path: Path) -> dict:
	# A: this whole helper was added here because the notebook used sampler.posterior() live, while this script loads pmc_final.state directly.
	"""Load the saved pocoMC state with compatibility shims for this environment."""
	sys.modules["__builtin__"] = builtins  # A: reproduce the Python-2-style module name expected by the serialized state.
	original_randomstate_ctor = numpy_random_pickle.__randomstate_ctor  # A: save the original RNG unpickler before patching.
	original_bitgen_ctor = numpy_random_pickle.__bit_generator_ctor  # A: save the original bit-generator unpickler before patching.

	def patched_randomstate_ctor(bit_generator_name: str = "MT19937", *args):
		# A: ignore extra pickle args so older numpy random-state payloads still load.
		return original_randomstate_ctor(bit_generator_name)

	def patched_bitgen_ctor(bit_generator_name: str = "MT19937", *args):
		# A: ignore extra pickle args for the bit-generator payload as well.
		return original_bitgen_ctor(bit_generator_name)

	numpy_random_pickle.__randomstate_ctor = patched_randomstate_ctor  # A: temporarily install the compatibility shim.
	numpy_random_pickle.__bit_generator_ctor = patched_bitgen_ctor  # A: temporarily install the compatibility shim.
	try:
		# A: the standalone script loads the saved posterior state from disk because it does not have the notebook sampler object in memory.
		with state_path.open("rb") as file_handle:
			return dill.load(file_handle)
	finally:
		# A: restore the original numpy pickle handlers after the state is loaded.
		numpy_random_pickle.__randomstate_ctor = original_randomstate_ctor
		numpy_random_pickle.__bit_generator_ctor = original_bitgen_ctor

def get_xmm_posterior_samples(state_path: Path, sample_count: int) -> np.ndarray:
	# A with N context: the notebook used sampler.posterior() directly, so this helper recreates that posterior sample access from the saved state.
	"""Return the tail of the saved posterior samples used for the notebook profiles."""
	state = load_xmm_state(state_path)  # A: load the serialized pocoMC state instead of relying on a live notebook sampler.
	samples = state["particles"].get("x", flat=True)  # A: this is how the flattened posterior chain is exposed inside the saved state.
	return samples[-sample_count:]  # A chosen to mirror the notebook's 1000-sample pressure-profile construction as closely as possible.

def compute_xmm_pressure_profiles(samples: np.ndarray, r_vals: np.ndarray) -> dict[str, np.ndarray]:
	# N lines 260-289: this helper is the script form of the notebook's main/west pressure-profile loops.
	"""Build the XMM median and 90% interval pressure profiles for both subclusters."""
	p_profiles_main = []  # N lines 263 and 266-270: same list-then-percentile pattern used for the main-cluster pressure profiles.
	p_profiles_sub = []  # N lines 282 and 285-289: same list-then-percentile pattern used for the west-subcluster pressure profiles.

	for sample in samples:
		# N lines 265 and 284: keep the notebook's posterior parameter order when unpacking each sample.
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

		# N line 266 supplies the pressure normalization and beta/2 usage; A adds the pixel-to-arcsec correction for r_s.
		main_profile = Xray_gNFW(
			r=r_vals,
			p_0=np.sqrt(a_main_mcmc / F_S_MAIN) * KT_MAIN_KEV,
			r_s=rs_main_mcmc * XMM_PIXEL_SCALE_ARCSEC,
			alpha=ALPHA_MAIN,
			beta=beta_main_mcmc / 2.0,
			gamma=GAMMA_MAIN,
		)
		# N line 285 supplies the west-subcluster normalization; A adds the pixel-to-arcsec correction for r_c.
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

		p_profiles_main.append(main_profile)  # N lines 266-268 condensed into a standalone-script pressure-profile list.
		p_profiles_sub.append(sub_profile)  # N lines 285-287 condensed into a standalone-script pressure-profile list.

	# N lines 270 and 289: use the same [5, 50, 95] percentile summary as the notebook plots.
	p_p5_main, p_p50_main, p_p95_main = np.percentile(p_profiles_main, [5, 50, 95], axis=0)
	p_p5_sub, p_p50_sub, p_p95_sub = np.percentile(p_profiles_sub, [5, 50, 95], axis=0)
	return {
		# A: package the notebook percentile arrays into named outputs so the plotting code stays readable.
		"P_p5_main": p_p5_main,
		"P_p50_main": p_p50_main,
		"P_p95_main": p_p95_main,
		"P_p5_sub": p_p5_sub,
		"P_p50_sub": p_p50_sub,
		"P_p95_sub": p_p95_sub,
	}


def format_radius_axes(axis: plt.Axes, kpc_per_arcsec: float) -> None:
	# C lines 79-93 and 99-113: this helper factors out the repeated axis-formatting code from the original two-panel script.
	"""Apply the shared arcsec and kpc radius axis formatting."""
	axis.set_xlim(10, 90)  # C lines 79 and 99: same x-range as the original comparison figure.
	axis.set_xlabel('r (")')  # C lines 80 and 100: same arcsec x-axis label.
	axis.set_xticks([10, 20, 30, 40, 50, 60, 70, 80, 90])  # C lines 84 and 104: same major tick locations.
	axis.set_xticklabels(["10", "20", "30", "40", "50", "60", "70", "80", "90"])  # C lines 85 and 105: same tick labels.

	# C lines 88-91 and 108-111: same secondary kpc axis, just wrapped in a helper here.
	axis_top = axis.secondary_xaxis("top", functions=(lambda x: x * kpc_per_arcsec, lambda x: x / kpc_per_arcsec))
	axis_top.set_xlabel("r (kpc)")
	axis_top.set_xticks([100, 200, 300, 400, 500, 600, 700])
	axis_top.set_xticklabels(["100", "200", "300", "400", "500", "600", "700"])


def main() -> None:
	# C lines 37-58: keep the SZ side exactly as in the original comparison script.
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
		p_0=sub_sz_params["P_0"].values[0],
		r_c=sub_sz_params["r_c_arcsec"].values[0],
		beta=sub_sz_params["beta"].values[0],
	)

	# N lines 260-289 plus A state-loading helpers: replace the single CSV XMM curves with posterior median and 90% interval curves.
	xmm_samples = get_xmm_posterior_samples(xmm_state_path, PROFILE_SAMPLE_COUNT)
	xmm_profiles = compute_xmm_pressure_profiles(xmm_samples, R_XMM)

	fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))  # C line 69: same two-panel figure layout as the original script.

	z = 1.189  # C line 72: same adopted MOO 1142 redshift as the original comparison plot.
	kpc_per_arcsec = Planck18.kpc_proper_per_arcmin(z).to(u.kpc / u.arcsec).value  # C line 73: same kpc/arcsec conversion.

	# C line 76 and following: preserve the original main-panel title/layout, but swap the red CSV XMM curve for the notebook-style posterior median and band.
	ax1.set_title("Main Cluster Pressure Profile")
	ax1.loglog(R_SZ, main_cluster_sz_profile, c="blue", label="M2")
	ax1.loglog(R_XMM, xmm_profiles["P_p50_main"], c="red", label="XMM")
	ax1.fill_between(R_XMM, xmm_profiles["P_p5_main"], xmm_profiles["P_p95_main"], color="red", alpha=0.3)
	ax1.set_ylabel("P (keV cm$^{-3}$)")
	ax1.legend()
	format_radius_axes(ax1, kpc_per_arcsec)

	# C line 96 and following: preserve the original west-panel layout, but again replace the single CSV XMM curve with posterior median and band.
	ax2.set_title("West Subcluster Pressure Profile")
	ax2.loglog(R_SZ, subcluster_sz_profile, c="blue", label="M2")
	ax2.loglog(R_XMM, xmm_profiles["P_p50_sub"], c="red", label="XMM")
	ax2.fill_between(R_XMM, xmm_profiles["P_p5_sub"], xmm_profiles["P_p95_sub"], color="red", alpha=0.3)
	ax2.set_ylabel("P (keV cm$^{-3}$)")
	ax2.legend()
	format_radius_axes(ax2, kpc_per_arcsec)

	plt.tight_layout()  # C line 115: same final layout cleanup as the original script.
	plt.savefig(location / "plots/MOO_1142_main+west_sz_xray_profiles_90CL.png", dpi=300)  # C line 116 with A filename change to distinguish the 90% CL output.


if __name__ == "__main__":
	main()  # A: standard standalone-script entry point; the notebook did not need this wrapper.