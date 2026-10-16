"""
Compare the SZ and X-ray pressure profiles of MOO 1142 and represent the error with 90% confidence intervals.

This is a copy of compare_xray_SZ_pressure_profiles_errors.py that additionally computes and plots the
full (n_r, n_r) covariance of each pressure-profile ensemble, to check whether the 90% CL bands (which
only probe the diagonal of this covariance) are hiding strong off-diagonal correlations between radial bins.

Created in 2026-05, covariance diagnostics added 2026-09.
"""
from __future__ import annotations

# Standard-library imports.
import builtins
import sys
from pathlib import Path

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

# Globally relevant constants and parameters.
R_values = np.arange(0.1, 100.0, 0.1)
PROFILE_SAMPLE_COUNT = 1000

# XMM profile information
xmm_fit_dir = location.parent / "XMM/Barbavara_fit_2026-04"
xmm_state_path = xmm_fit_dir / "gnfw_circ+beta_circ_acfixed_NESTED/pmc_final.state"
ALPHA_MAIN = 2.26
GAMMA_MAIN = 0.465
F_S_MAIN = 1.632e4
F_S_SUB = 2.597e5
KT_MAIN_KEV = 6.76
KT_SUB_KEV = 7.13
XMM_PIXEL_SCALE_ARCSEC = 2.5

# M2 profile information
m2_fit_dir = location.parent / "M2/MOO1142_Xray_priors_subcluster"
m2_fit_cov_path = m2_fit_dir / "par_results_final_fit.dill"
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

def compute_profile_covariance(p_profiles: list[np.ndarray]) -> np.ndarray:
	# Provenance: colleague's suggested sanity check on 2026-09.
	# Intent: expose the full (n_r, n_r) covariance across radial bins, since the percentile bands
	# only summarize the diagonal and can hide strong bin-to-bin correlations.
	"""Return the (n_r, n_r) covariance matrix of a pressure-profile ensemble."""
	return np.cov(np.asarray(p_profiles).T)

def plot_profile_covariance(cov: np.ndarray, r_vals: np.ndarray, title: str, output_path: Path) -> None:
	"""Save an imshow diagnostic plot of a radial-bin covariance matrix."""
	fig, ax = plt.subplots(figsize=(6, 5))
	extent = (r_vals[0], r_vals[-1], r_vals[-1], r_vals[0])
	log_cov = np.where(np.abs(cov) > 0.0, np.log10(np.abs(cov)), np.nan)
	im = ax.imshow(log_cov, extent=extent, cmap="viridis")
	ax.set_title(title)
	ax.set_xlabel('r (")')
	ax.set_ylabel('r (")')
	fig.colorbar(im, ax=ax, label="log10(abs(covariance))")
	fig.savefig(output_path, dpi=300)
	plt.close(fig)

def compute_xmm_pressure_profiles(samples: np.ndarray, r_vals: np.ndarray) -> tuple[list[np.ndarray], list[np.ndarray]]:
	# Provenance: adapted from radial_profile_and_uncert_n.ipynb, posterior-profile construction.
	# Intent: turn each posterior sample into a main-cluster and west-subcluster pressure profile.
	"""Build the XMM main-cluster and west-subcluster pressure-profile ensembles."""
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

	return p_profiles_main, p_profiles_sub

def load_M2_covariance(cov_path: Path) -> dict[str, pd.DataFrame]:
	"""Load covariance of M2 fit."""
	with open(cov_path, "rb") as f:
		M2_fit_par_dict = dill.load(f)

	return M2_fit_par_dict

def report_M2_fit_parameters(
	M2_fit_par_dict: dict[str, pd.DataFrame],
	expected_values: np.ndarray,
	nonzero: tuple[np.ndarray, ...],
	cov: np.ndarray,
) -> None:
	"""Print constrained M2 fit parameters and their 1-sigma uncertainties."""
	mydiag = np.diag(cov)
	print("M2 fit parameters and 1-sigma uncertainties:")
	for pn, expected, unc in zip(M2_fit_par_dict["par_names"][nonzero], expected_values[nonzero], np.sqrt(mydiag)):
		print(pn, ':', np.round(expected, 4), '+/-', np.round(unc, 4))

def compute_M2_pressure_profiles(
	M2_fit_par_dict: dict[str, pd.DataFrame],
	r_vals: np.ndarray,
	*,
	report_fit_parameters: bool = True,
) -> tuple[list[np.ndarray], list[np.ndarray]]:
	expected_values = M2_fit_par_dict["parameters"]   # mean if PURELY Gaussian.

	# Check for zero rows in the covariance matrix, which would indicate unconstrained parameters/parameters that were not fit.
	rowsum = np.sum(M2_fit_par_dict["cov"], axis=0)
	#print(rowsum.shape)
	nonzero = np.where(rowsum != 0)

	# Extract the submatrix of the covariance corresponding to the nonzero rows/columns.
	cov_intermediate = M2_fit_par_dict["cov"][nonzero]
	cov = cov_intermediate[:,nonzero[0]]
	#print(cov.shape)

	# Optional: print out parameter values and errors
	if report_fit_parameters:
		report_M2_fit_parameters(M2_fit_par_dict, expected_values, nonzero, cov)

	# Generate samples from the multivariate normal distribution defined by the expected values and covariance.
	L = np.linalg.cholesky(cov)
	rng = np.random.default_rng()
	ngoodpar = cov.shape[0]
	rand_gen = rng.normal(size=(PROFILE_SAMPLE_COUNT, ngoodpar))  # Unscaled random values. (Or just scaled to unity)
	rotated = np.dot(L,rand_gen.T)  # Both scales the random values and rotates according to covariance (correlation). 
	# rotated is "Rotating about zero" so lets add the expected values back in to get the full distribution of parameter values
	generated = rotated.T + expected_values[nonzero]   # Add means back in, full distribution.

	# Optional: print out the generated parameter values for inspection.
	#print(M2_fit_par_dict['par_names'][nonzero])

	p_profiles_main = []
	p_profiles_sub = []

	for sample in generated:
		#print(sample.shape)
		pressure_main = gNFW(r_vals, sample[3], sample[4], sample[6], sample[7], sample[5])  # alpha, beta, gamma in gNFW call.
		pressure_sub = iso_beta(r_vals, sample[2], sample[0], sample[1])
		
		p_profiles_main.append(pressure_main)
		p_profiles_sub.append(pressure_sub)

	return p_profiles_main, p_profiles_sub

# -------------------------------------------------------------------------------------------- #
def main() -> None:
	cov_plot_dir = location / "plots/pressure/covariance"
	cov_plot_dir.mkdir(parents=True, exist_ok=True)

	# Read in the covariance of the M2 fit and summarize the pressure-profile uncertainty.
	m2_fit_cov = load_M2_covariance(m2_fit_cov_path)
	m2_p_profiles_main, m2_p_profiles_sub = compute_M2_pressure_profiles(m2_fit_cov, R_values)

	# Load the XMM posterior samples and compute the pressure profiles for each sample, then summarize with the median and 90% interval curves.
	xmm_samples = get_xmm_posterior_samples(xmm_state_path, PROFILE_SAMPLE_COUNT)
	xmm_p_profiles_main, xmm_p_profiles_sub = compute_xmm_pressure_profiles(xmm_samples, R_values)

	# Sanity-check the radial-bin covariance structure behind each 90% CL band.
	m2_cov_main = compute_profile_covariance(m2_p_profiles_main)
	m2_cov_sub = compute_profile_covariance(m2_p_profiles_sub)
	xmm_cov_main = compute_profile_covariance(xmm_p_profiles_main)
	xmm_cov_sub = compute_profile_covariance(xmm_p_profiles_sub)

	plot_profile_covariance(m2_cov_main, R_values, "M2 Main Cluster Profile Covariance", cov_plot_dir / "M2_main_cov.png")
	plot_profile_covariance(m2_cov_sub, R_values, "M2 West Subcluster Profile Covariance", cov_plot_dir / "M2_sub_cov.png")
	plot_profile_covariance(xmm_cov_main, R_values, "XMM Main Cluster Profile Covariance", cov_plot_dir / "XMM_main_cov.png")
	plot_profile_covariance(xmm_cov_sub, R_values, "XMM West Subcluster Profile Covariance", cov_plot_dir / "XMM_sub_cov.png")

if __name__ == "__main__":
	main()
