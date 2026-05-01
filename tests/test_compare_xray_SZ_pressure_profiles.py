import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np
import pytest
import scipy.special
from numpy.testing import assert_allclose

from profile_fitting.compare_xray_sz.profile_models import gNFW, iso_beta
from profile_fitting.compare_xray_sz.compare_xray_SZ_pressure_profiles import (
    build_pressure_profiles,
    load_profile_parameters,
    make_pressure_profile_figure,
)
from profile_fitting.compare_xray_sz.compare_xray_SZ_pressure_profiles_90CL import (
    ALPHA_MAIN,
    F_S_MAIN,
    F_S_SUB,
    GAMMA_MAIN,
    KT_MAIN_KEV,
    KT_SUB_KEV,
    PROFILE_SAMPLE_COUNT,
    XMM_PIXEL_SCALE_ARCSEC,
    compute_xmm_pressure_profiles,
    get_xmm_posterior_samples,
    xmm_state_path,
)


def _require_xmm_state():
    if not xmm_state_path.exists():
        pytest.skip('Requires the saved XMM posterior state file')


def _build_xmm_pressure_ensembles(samples, r):
    main_profiles = []
    sub_profiles = []

    for sample in samples:
        (
            _x_main,
            _y_main,
            _x_sub,
            _y_sub,
            rs_main_mcmc,
            rc_sub_mcmc,
            beta_main_mcmc,
            beta_sub_mcmc,
            a_main_mcmc,
            a_sub_mcmc,
            _a_bkg_mcmc,
        ) = sample

        main_profiles.append(
            gNFW(
                r=r,
                p_0=np.sqrt(a_main_mcmc / F_S_MAIN) * KT_MAIN_KEV,
                r_s=rs_main_mcmc * XMM_PIXEL_SCALE_ARCSEC,
                alpha=ALPHA_MAIN,
                beta=beta_main_mcmc / 2.0,
                gamma=GAMMA_MAIN,
            )
        )
        sub_profiles.append(
            iso_beta(
                r=r,
                p_0=np.sqrt(
                    a_sub_mcmc * scipy.special.gamma(3.0 * beta_sub_mcmc)
                    / (F_S_SUB * scipy.special.gamma(3.0 * beta_sub_mcmc - 0.5))
                )
                * KT_SUB_KEV,
                r_c=rc_sub_mcmc * XMM_PIXEL_SCALE_ARCSEC,
                beta=beta_sub_mcmc,
            )
        )

    return np.asarray(main_profiles), np.asarray(sub_profiles)


def test_gnfw_matches_previous_profile_form():
    """Check that the shared gNFW helper reproduces the previous analytic profile form."""
    r = np.array([1.0, 5.0, 10.0, 20.0])
    params = dict(p_0=3.2, r_s=15.0, alpha=1.1, beta=5.4, gamma=0.3)
    expected = params['p_0'] * ((r / params['r_s']) ** -params['gamma']) * (
        1 + (r / params['r_s']) ** params['alpha']
    ) ** ((params['gamma'] - params['beta']) / params['alpha'])

    assert_allclose(gNFW(r, **params), expected)


def test_iso_beta_is_central_value_at_zero_radius_and_declines_outward():
    """Verify that the isothermal beta profile peaks at the center and decreases with radius."""
    r = np.array([0.0, 10.0, 30.0])
    profile = iso_beta(r, p_0=5.0, r_c=12.0, beta=0.8)

    assert profile[0] == 5.0
    assert np.all(np.diff(profile) < 0)


def test_build_pressure_profiles_returns_expected_arrays_from_lightweight_params():
    """Ensure profile construction maps lightweight parameter inputs to the expected model outputs."""
    r = np.array([1.0, 2.0, 5.0])
    params = {
        'main_sz': {'P_0': [1.0], 'r_s_arcsec': [10.0], 'alpha': [1.0], 'beta': [5.0], 'gamma': [0.5]},
        'main_xray': {'P_0': [2.0], 'r_s_arcsec': [12.0], 'alpha': [1.2], 'beta': [4.0], 'gamma': [0.4]},
        'sub_sz': {'P_0': [3.0], 'r_c_arcsec': [8.0], 'beta': [0.8]},
        'sub_xray': {'P_0': [4.0], 'r_c_arcsec': [9.0], 'beta': [0.9]},
    }

    profiles = build_pressure_profiles(r, params)

    assert_allclose(profiles['main_sz'], gNFW(r, 1.0, 10.0, 1.0, 5.0, 0.5))
    assert_allclose(profiles['main_xray'], gNFW(r, 2.0, 12.0, 1.2, 4.0, 0.4))
    assert_allclose(profiles['sub_sz'], iso_beta(r, 3.0, 8.0, 0.8))
    assert_allclose(profiles['sub_xray'], iso_beta(r, 4.0, 9.0, 0.9))


def test_make_pressure_profile_figure_creates_two_panels():
    """Confirm the plotting helper builds the expected two-panel pressure-profile figure."""
    r = np.array([10.0, 20.0, 40.0])
    profiles = {
        'main_sz': np.array([3.0, 2.0, 1.0]),
        'main_xray': np.array([2.5, 1.8, 0.9]),
        'sub_sz': np.array([1.5, 1.0, 0.5]),
        'sub_xray': np.array([1.4, 0.9, 0.4]),
    }

    fig = make_pressure_profile_figure(r, profiles)

    assert len(fig.axes) == 2
    assert fig.axes[0].get_title() == 'Main Cluster Pressure Profile'
    assert fig.axes[1].get_title() == 'Sub-cluster Pressure Profile'
    assert len(fig.axes[0].lines) == 2
    assert len(fig.axes[1].lines) == 2
    plt.close(fig)


def test_posterior_median_parameters_match_legacy_xray_parameter_files():
    """Verify the posterior-derived median XMM parameters agree with the legacy CSV values."""
    _require_xmm_state()
    params = load_profile_parameters()
    samples = get_xmm_posterior_samples(xmm_state_path, PROFILE_SAMPLE_COUNT)
    medians = np.median(samples, axis=0)

    # Check that the median main-cluster amplitude from the posterior converts back to
    # the same physical pressure normalization stored in the legacy X-ray parameter file.
    assert_allclose(
        np.sqrt(medians[8] / F_S_MAIN) * KT_MAIN_KEV,
        params['main_xray']['P_0'].values[0],
        rtol=0.01,
    )
    # Check that the median main-cluster scale radius from the posterior, after converting
    # from fit pixels to arcseconds, matches the legacy X-ray radius value.
    assert_allclose(
        medians[4] * XMM_PIXEL_SCALE_ARCSEC,
        params['main_xray']['r_s_arcsec'].values[0],
        rtol=0.01,
    )
    # Check that the main-cluster slope parameter is being translated the same way in both
    # code paths, including the notebook-to-script beta/2 convention.
    assert_allclose(
        medians[6] / 2.0,
        params['main_xray']['beta'].values[0],
        rtol=0.01,
    )
    # Check that the median west-subcluster amplitude from the posterior reproduces the
    # legacy X-ray pressure normalization after applying the same isobeta conversion.
    assert_allclose(
        np.sqrt(
            medians[9] * scipy.special.gamma(3.0 * medians[7])
            / (F_S_SUB * scipy.special.gamma(3.0 * medians[7] - 0.5))
        )
        * KT_SUB_KEV,
        params['sub_xray']['P_0'].values[0],
        rtol=0.01,
    )
    # Check that the median west-subcluster core radius from the posterior, converted from
    # pixels to arcseconds, matches the radius listed in the legacy X-ray file.
    assert_allclose(
        medians[5] * XMM_PIXEL_SCALE_ARCSEC,
        params['sub_xray']['r_c_arcsec'].values[0],
        rtol=0.01,
    )
    # Check that the west-subcluster beta slope from the posterior matches the legacy X-ray
    # beta value directly, since this branch does not use the beta/2 conversion.
    assert_allclose(
        medians[7],
        params['sub_xray']['beta'].values[0],
        rtol=0.01,
    )


def test_legacy_xray_curves_fall_inside_posterior_90_percent_bands():
    """Ensure the old CSV-driven XMM curves are consistent with the posterior-derived 90% bands."""
    _require_xmm_state()
    r = np.arange(0.1, 100.0, 0.1)
    params = load_profile_parameters()
    legacy_profiles = build_pressure_profiles(r, params)
    samples = get_xmm_posterior_samples(xmm_state_path, PROFILE_SAMPLE_COUNT)
    posterior_profiles = compute_xmm_pressure_profiles(samples, r)

    assert np.all(legacy_profiles['main_xray'] >= posterior_profiles['P_p5_main'])
    assert np.all(legacy_profiles['main_xray'] <= posterior_profiles['P_p95_main'])
    assert np.all(legacy_profiles['sub_xray'] >= posterior_profiles['P_p5_sub'])
    assert np.all(legacy_profiles['sub_xray'] <= posterior_profiles['P_p95_sub'])


def test_posterior_90_percent_bands_match_independent_profile_ensemble_coverage():
    """Check that the plotted 90% bands enclose 90% of independently rebuilt profiles pointwise."""
    _require_xmm_state()
    r = np.arange(0.1, 100.0, 0.1)
    samples = get_xmm_posterior_samples(xmm_state_path, PROFILE_SAMPLE_COUNT)
    main_profiles, sub_profiles = _build_xmm_pressure_ensembles(samples, r)
    posterior_profiles = compute_xmm_pressure_profiles(samples, r)

    assert np.all(posterior_profiles['P_p5_main'] <= posterior_profiles['P_p50_main'])
    assert np.all(posterior_profiles['P_p50_main'] <= posterior_profiles['P_p95_main'])
    assert np.all(posterior_profiles['P_p5_sub'] <= posterior_profiles['P_p50_sub'])
    assert np.all(posterior_profiles['P_p50_sub'] <= posterior_profiles['P_p95_sub'])

    main_coverage = (
        (main_profiles >= posterior_profiles['P_p5_main'])
        & (main_profiles <= posterior_profiles['P_p95_main'])
    ).mean(axis=0)
    sub_coverage = (
        (sub_profiles >= posterior_profiles['P_p5_sub'])
        & (sub_profiles <= posterior_profiles['P_p95_sub'])
    ).mean(axis=0)

    assert_allclose(main_coverage, 0.9, atol=1e-12)
    assert_allclose(sub_coverage, 0.9, atol=1e-3)