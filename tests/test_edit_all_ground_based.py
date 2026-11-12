from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np
from astropy.table import Table
from numpy.testing import assert_allclose


SCRIPT = (Path(__file__).resolve().parents[1] / 'members' / 'specz' /
          'scripts' / 'edit_allGroundBased.py')
INPUT_NAME = '2026-10-19.MOO_1142+1527.speczs.internal.groundBasedAll.fits'
OUTPUT_NAME = '2026-11.MOO_1142+1527.speczs.internal.groundBasedAll.fits'


def run_script(tmp_path, rows):
    scripts = tmp_path / 'scripts'
    extra = tmp_path / 'extra'
    scripts.mkdir()
    extra.mkdir()
    script = scripts / SCRIPT.name
    shutil.copy(SCRIPT, script)
    table = Table(rows=rows, names=('Source', 'Q', 'Type'))
    table['RA'] = np.full(len(rows), 175.0)
    table['Dec'] = np.full(len(rows), 15.0)
    table['Spec-z'] = np.full(len(rows), 1.0)
    table.write(extra / INPUT_NAME)
    result = subprocess.run(
        [sys.executable, str(script)],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True,
    )
    return result, tmp_path / OUTPUT_NAME


def test_source_and_quality_dependent_errors(tmp_path):
    rows = [
        ('KG', 1, 'A'), ('KG', 2, 'E'), ('KG', 3, 'A+E'), ('KG', 3, 'star'),
        ('PN', 3, 'other'), ('PD', 3, 'star'),
        ('PN', 2, 'A'), ('PN', 2, 'E'), ('PN', 2, 'A+E'),
        ('PD', 2, 'A'), ('PD', 2, 'E'), ('PD', 2, 'A+E'),
        ('PN', 1, 'A'), ('PD', 2, 'star'), ('PD', 4, 'E'),
    ]
    result, output = run_script(tmp_path, rows)

    assert result.returncode == 0, result.stderr
    assert result.stdout.count(
        "you have a Palomar source that doesn't fit the standard boiler plate error estimation."
    ) == 3
    assert 'Row 12: Source=PN, Q=1, Type=A' in result.stdout
    table = Table.read(output)
    expected = [0.003, 0.0005, 0.0005, 0.0] + [0.002] * 8 + [np.nan] * 3
    assert_allclose(table['Spec-z_err_hi'], expected, equal_nan=True)
    assert_allclose(table['Spec-z_err_lo'], expected, equal_nan=True)


def test_unknown_non_palomar_type_still_fails(tmp_path):
    result, output = run_script(tmp_path, [('KG', 3, 'other')])

    assert result.returncode != 0
    assert "Unknown Type values: ['other']" in result.stderr
    assert not output.exists()
