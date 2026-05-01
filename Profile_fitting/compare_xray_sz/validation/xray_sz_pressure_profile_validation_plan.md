# X-ray/SZ Pressure Profile Validation Plan

## Scope

This note summarizes what has already been validated for the combined X-ray and SZ pressure-profile workflow, what concerns remain in the current implementation, and what the next steps should be.

Relevant files:

- `profile_fitting/compare_xray_sz/compare_xray_SZ_pressure_profiles.py`
- `profile_fitting/compare_xray_sz/compare_xray_SZ_pressure_profiles_90CL.py`
- `profile_fitting/XMM/Barbavara_fit_2026-04/radial_profile_and_uncert_n.ipynb`
- `tests/test_compare_xray_SZ_pressure_profiles.py`

## What Has Been Accomplished

### 1. Verified that the analytic profile forms match the original script logic

The existing test file confirms that the shared `gNFW` and `iso_beta` helpers reproduce the expected analytic expressions used by the earlier comparison scripts.

Result:

- The common profile-model helpers are behaving consistently with the original pressure-profile script.

### 2. Verified that the legacy CSV-driven pressure profiles are still constructed correctly

The tests confirm that `build_pressure_profiles(...)` maps the loaded parameter values to the expected pressure-profile arrays for both the main cluster and the west subcluster.

Result:

- The original CSV-based line construction is still behaving as intended.

### 3. Verified that the posterior-derived XMM median parameters reproduce the legacy X-ray parameter files

New tests were added to compare the median of the saved posterior samples against the legacy X-ray CSV parameter files after applying the same conversions used in the notebook port.

Checks performed:

- Main-cluster pressure normalization `P_0`
- Main-cluster scale radius `r_s_arcsec`
- Main-cluster slope parameter `beta` after the notebook-style `beta / 2` conversion
- West-subcluster pressure normalization `P_0`
- West-subcluster core radius `r_c_arcsec`
- West-subcluster slope parameter `beta`

Result:

- The posterior-derived median parameters match the legacy X-ray CSV values at about the 1% level.
- This strongly suggests that the standalone script is converting posterior samples into physical profile parameters in a way that is compatible with the original X-ray profile files.

### 4. Verified that the old X-ray curves lie inside the new 90% posterior bands

New tests compare the legacy CSV-driven X-ray curves against the posterior-derived 5th to 95th percentile profile envelopes.

Result:

- The old X-ray curves lie entirely inside the new 90% posterior bands across the tested radius range.
- This indicates that the new banded plot is broadly compatible with the earlier single-line X-ray representation.

### 5. Verified that the plotted 90% bands actually behave like 90% percentile bands

New tests independently rebuild the pressure-profile ensembles from the saved posterior samples, then check what fraction of sampled profiles fall between the plotted lower and upper bands at each radius.

Result:

- The pointwise containment is about 90% as expected.
- This confirms that the 5th to 95th percentile band construction is numerically consistent with the intended central 90% posterior interval.

### 6. Confirmed that the saved state does not currently introduce a weighting mismatch

The original notebook called `sampler.posterior()` and exposed `weights`, so one concern was that the standalone script might be dropping important weighting information when reading the saved state directly.

Result:

- For the current `pmc_final.state`, the stored particle weights are effectively uniform.
- This means the current implementation is not obviously wrong on that point for this specific state.

### 7. Re-ran the focused validation tests successfully

The targeted test file was run after the new validation tests and documentation comments were added.

Result:

- `tests/test_compare_xray_SZ_pressure_profiles.py` passed fully.

## Current Red Flags and Open Concerns

### 1. The standalone script does not exactly reproduce the notebook's sample-selection logic

The notebook builds profile ensembles using `samples[-i]` inside a loop over 1000 samples. The standalone script instead takes the last 1000 samples as a contiguous slice.

Why this matters:

- These are not exactly the same sample sets.
- For the current saved state, the difference is small for most derived curves, but it is not identically the same logic.
- The largest measured discrepancy showed up in the west-subcluster lower percentile envelope.

Interpretation:

- This is the clearest concrete divergence from the original notebook logic.
- It is not necessarily a fatal bug, but it should be made explicit and either preserved intentionally or corrected intentionally.

### 2. The standalone script evaluates XMM profiles beyond the notebook's original radius range

The notebook constructs profiles on a radius grid from 1 to 60 arcsec, while the standalone script evaluates XMM profiles from 0.1 to 100 arcsec and then plots out to 90 arcsec.

Why this matters:

- The combined script is using the model beyond the notebook's original displayed radius domain.
- The XMM band shown outside the notebook range is therefore an extrapolation of the same fitted model, not a direct reproduction of the notebook output.

Interpretation:

- This may be acceptable for the comparison figure, but it is a real assumption that should be documented.

### 3. Several conversion constants are duplicated manually in the standalone script

The combined script hard-codes:

- `ALPHA_MAIN`
- `GAMMA_MAIN`
- `F_S_MAIN`
- `F_S_SUB`
- `KT_MAIN_KEV`
- `KT_SUB_KEV`
- `XMM_PIXEL_SCALE_ARCSEC`

Why this matters:

- These values are coupled to the fit setup and physical interpretation.
- If the fit is rerun with updated assumptions or different data products, the comparison script can quietly become stale.

Interpretation:

- This is a maintainability and reproducibility risk more than a proven current bug.

### 4. The term `90CL` is a naming convention, not a derived quantity from the code

The code constructs the band using the 5th and 95th percentiles of the posterior profile ensemble.

Why this matters:

- This is best interpreted statistically as a central 90% posterior interval.
- The filename shorthand `90CL` is understandable, but it should not be confused with a separately derived frequentist confidence construction.

Interpretation:

- This is mostly a terminology issue, not a numerical bug.

## Recommended Next Steps

### High Priority

1. Decide whether the goal is exact notebook reproduction or a cleaned-up standalone implementation.

Questions to answer:

- Should the script intentionally mimic the notebook's `samples[-i]` indexing behavior exactly?
- Or should it intentionally use the last 1000 samples as a contiguous slice and document that this is a cleaned-up choice?

Suggested action:

- Make the sample-selection rule explicit in code comments and tests.

2. Externalize the XMM conversion assumptions.

Suggested action:

- Read the pixel scale and other fit-dependent constants from a configuration or metadata source instead of hard-coding them directly in the comparison script.
- If some quantities truly are fixed analysis choices rather than fit metadata, document their provenance in a dedicated constants block or sidecar file.

### Medium Priority

3. Add one explicit regression test for notebook-equivalent sample selection.

Suggested action:

- Add a test that compares the current sample-selection method against the notebook-style selection and forces the project to choose one behavior intentionally.

4. Decide whether the plotted XMM radius range should match the notebook range.

Suggested action:

- Either restrict the XMM band to the notebook domain or document clearly that the larger domain is a deliberate model extrapolation for comparison with the SZ profile.

### Optional but Useful

5. Add an export or diagnostic script for direct inspection of plotted arrays.

Suggested action:

- Save the legacy X-ray curve, posterior median curve, and lower and upper percentile curves to CSV for direct comparison outside matplotlib.

6. Clarify the statistical language in outputs.

Suggested action:

- If the project wants strict terminology, refer to the plotted XMM band as the central 90% posterior interval or central 90% credible interval.
- If the project prefers the existing shorthand, keep `90CL` but note what it means computationally.

## Suggested Immediate Action

The most practical next move is:

1. Resolve the sample-selection question.
2. Move the fit-dependent constants out of the script body.
3. Decide whether the XMM radius range should be reproduction-based or comparison-oriented.

These three choices will remove the biggest current ambiguities without changing the overall structure of the comparison workflow.