from pathlib import Path
from datetime import date
import shutil
import subprocess
import sys

import numpy as np
import pytest
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
from numpy.testing import assert_allclose, assert_array_equal

from members.specz.cross_match_ground_JWST.cross_match_utils import (
    closest_unique_matches,
    report_multiple_ground_candidates,
)


def coordinates(offsets):
    return SkyCoord(ra=np.array(offsets) * u.arcsec, dec=np.zeros(len(offsets)) * u.deg)


def test_uses_closest_available_alternative():
    jwst = coordinates([0.1, -0.3, 10.0])
    ground = coordinates([0.0, -1.2, 10.5])

    j, g, sep = closest_unique_matches(jwst, ground, 1.5)

    assert_array_equal(j, [0, 2, 1])
    assert_array_equal(g, [0, 2, 1])
    assert_allclose(sep, [0.1, 0.5, 0.9], atol=1e-8)


def test_closest_first_does_not_maximize_pair_count():
    jwst = coordinates([0.0, -0.7])
    ground = coordinates([-0.3, 0.5])

    j, g, sep = closest_unique_matches(jwst, ground, 0.6)

    assert_array_equal(j, [0])
    assert_array_equal(g, [0])
    assert_allclose(sep, [0.3], atol=1e-8)


def test_no_matches_and_radius_cut():
    jwst = coordinates([0.0, 10.0])
    ground = coordinates([0.2, 10.4])

    j, g, sep = closest_unique_matches(jwst, ground, 0.3)
    assert_array_equal(j, [0])
    assert_array_equal(g, [0])
    assert_allclose(sep, [0.2], atol=1e-8)

    j, g, sep = closest_unique_matches(jwst, ground, 0.1)
    assert len(j) == len(g) == len(sep) == 0


def test_equal_separations_use_row_indices():
    jwst = coordinates([0.0, 0.0])
    ground = coordinates([0.0, 0.0])

    j, g, sep = closest_unique_matches(jwst, ground, 1.5)

    assert_array_equal(j, [0, 1])
    assert_array_equal(g, [0, 1])
    assert_array_equal(sep, [0.0, 0.0])


def test_reports_alternative_match_and_ground_assignment(capsys):
    jwst = coordinates([0.1, -0.3])
    ground = coordinates([0.0, -1.2])
    j, g, _ = closest_unique_matches(jwst, ground, 1.5)

    report_multiple_ground_candidates(
        jwst, ground, 1.5, ['J0', 'J1'], ['G0', 'G1'], j, g
    )

    output = capsys.readouterr().out
    assert 'JWST J1 / ground G0 (sep = 0.300"): ground matched to JWST J0 (sep = 0.100")' in output
    assert 'JWST J1 / ground G1 (sep = 0.900"): accepted' in output
    assert 'JWST J1 remains unmatched' not in output


def test_reports_multiple_candidates_when_jwst_is_unmatched(capsys):
    jwst = coordinates([0.0, -0.3, 0.5])
    ground = coordinates([-0.3, 0.5])
    j, g, _ = closest_unique_matches(jwst, ground, 0.6)

    report_multiple_ground_candidates(
        jwst, ground, 0.6, ['J0', 'J1', 'J2'], ['G0', 'G1'], j, g
    )

    output = capsys.readouterr().out
    assert 'JWST J0 / ground G0 (sep = 0.300"): ground matched to JWST J1 (sep = 0.000")' in output
    assert 'JWST J0 / ground G1 (sep = 0.500"): ground matched to JWST J2 (sep = 0.000")' in output
    assert 'JWST J0 remains unmatched' in output


def test_single_candidate_does_not_produce_report(capsys):
    jwst = coordinates([0.0])
    ground = coordinates([0.2])
    j, g, _ = closest_unique_matches(jwst, ground, 1.5)

    report_multiple_ground_candidates(jwst, ground, 1.5, ['J0'], ['G0'], j, g)

    assert capsys.readouterr().out == ''


@pytest.mark.parametrize(
    'bcg_ra, expected_id',
    [(0.00011, 'J0'), (0.02001, 'G1'), (0.03, None)],
)
def test_combined_catalog_bcg_flag(tmp_path, bcg_ra, expected_id):
    source = Path(__file__).resolve().parents[1] / 'members/specz/cross_match_ground_JWST'
    specz = tmp_path / 'specz'
    work = specz / 'cross_match_ground_JWST'
    work.mkdir(parents=True)
    bcg_dir = tmp_path / 'BCG'
    bcg_dir.mkdir()
    for name in ('03_do_ground_JWST_cross_match.py', 'cross_match_utils.py'):
        shutil.copy2(str(source / name), str(work / name))

    def catalog(ids, ra):
        return Table({
            'ID': ids, 'RA': ra, 'Dec': [0.0, 0.0],
            'Spec-z': [1.0, 1.1], 'Q': ['3', '3'],
            'Type': ['galaxy', 'galaxy'], 'Source': ['survey', 'survey'],
        })

    catalog(['J0', 'J1'], [0.0001, 0.01]).write(
        str(specz / '2026-10-15.MOO_1142+1527.speczs.JWST.fits')
    )
    catalog(['G0', 'G1'], [0.00011, 0.02]).write(
        str(specz / '2026-11.MOO_1142+1527.speczs.allGroundBased.fits')
    )
    Table({'ra_spz': [bcg_ra], 'dec_spz': [0.0]}).write(
        str(bcg_dir / 'MOO_1142_bcg_info.txt'), format='ascii.csv'
    )

    result = subprocess.run(
        [sys.executable, str(work / '03_do_ground_JWST_cross_match.py')],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True,
        check=True,
    )
    outputs = list(specz.glob('*.MOO_1142+1527.speczs.all.fits'))
    assert len(outputs) == 1
    assert outputs[0].name == '{}.MOO_1142+1527.speczs.all.fits'.format(date.today().isoformat())
    combined = Table.read(str(outputs[0]))
    assert combined.colnames[-1] == 'BCG'
    assert 'sep_BCG' not in combined.colnames
    assert combined['BCG'].dtype.kind == 'b'
    assert list(combined['ID'][combined['BCG']]) == ([] if expected_id is None else [expected_id])
    assert set(combined['ID']) == {'J0', 'J1', 'G1'}
    if expected_id is None:
        assert 'WARNING: no source within 1.5" of the BCG' in result.stdout


def test_good_quality_catalog_excludes_only_numeric_one(tmp_path):
    source = Path(__file__).resolve().parents[1] / 'members/specz/cross_match_ground_JWST'
    specz = tmp_path / 'specz'
    work = specz / 'cross_match_ground_JWST'
    work.mkdir(parents=True)
    bcg_dir = tmp_path / 'BCG'
    bcg_dir.mkdir()
    for name in ('03_do_ground_JWST_cross_match.py', 'cross_match_utils.py'):
        shutil.copy2(str(source / name), str(work / name))

    def catalog(ids, ra, q):
        return Table({
            'ID': ids, 'RA': ra, 'Dec': np.zeros(len(ids)),
            'Spec-z': np.arange(len(ids), dtype=float), 'Q': q,
            'Type': ['galaxy'] * len(ids), 'Source': ['survey'] * len(ids),
        })

    catalog(['J0', 'J1', 'J2', 'J3'], [0.0, 0.01, 0.02, 0.03],
            ['A/B', '1', '1.00', '2']).write(
        str(specz / '2026-10-15.MOO_1142+1527.speczs.JWST.fits')
    )
    catalog(['G0', 'G1', 'G2'], [0.04, 0.05, 0.06],
            [1.0, 2.0, np.nan]).write(
        str(specz / '2026-11.MOO_1142+1527.speczs.allGroundBased.fits')
    )
    Table({'ra_spz': [0.01], 'dec_spz': [0.0]}).write(
        str(bcg_dir / 'MOO_1142_bcg_info.txt'), format='ascii.csv'
    )

    subprocess.run(
        [sys.executable, str(work / '03_do_ground_JWST_cross_match.py')],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True,
        check=True,
    )
    combined = Table.read(str(next(specz.glob('*.MOO_1142+1527.speczs.all.fits'))))
    good = Table.read(str(next(specz.glob('*.MOO_1142+1527.speczs.all.goodQ.fits'))))
    assert (specz / '{}.MOO_1142+1527.speczs.all.goodQ.fits'.format(date.today().isoformat())).exists()
    assert set(combined['ID']) == {'J0', 'J1', 'J2', 'J3', 'G0', 'G1', 'G2'}
    assert list(good['ID']) == [row['ID'] for row in combined if row['ID'] in {'J0', 'J3', 'G1', 'G2'}]
    assert good.colnames == combined.colnames
    assert not np.any(good['BCG'])
    assert list(combined['ID'][combined['BCG']]) == ['J1']
