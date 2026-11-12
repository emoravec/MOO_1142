"""
Cross-match JWST and ground sources within 1.5", accepting the closest available
pairs first. Write the matches and a combined table in which the JWST row is kept
for each duplicate. Flag the closest combined source within 1.5" of the BCG
and write an additional catalog excluding sources with Q equal to 1.
"""
import numpy as np
from datetime import date
from pathlib import Path

from astropy.coordinates import SkyCoord
from astropy.table import Table, hstack, vstack
from astropy import units as u
from cross_match_utils import closest_unique_matches, report_multiple_ground_candidates
# -------------------------------------------------------------------------------------------- #
RADIUS = 1.5  # arcsec
BCG_RADIUS = 1.5  # arcsec
EXPECTED_DUPS = 8

specz_dir = Path(__file__).resolve().parent.parent
jwst = Table.read(specz_dir / '2026-10-15.MOO_1142+1527.speczs.JWST.fits')
ground = Table.read(specz_dir / '2026-11.MOO_1142+1527.speczs.allGroundBased.fits')

jwst_coord = SkyCoord(jwst['RA'] * u.deg, jwst['Dec'] * u.deg, frame='icrs')
ground_coord = SkyCoord(ground['RA'] * u.deg, ground['Dec'] * u.deg, frame='icrs')

# Search all pairs within the radius, accepting the closest unused sources first.
# A source whose closest candidate is already used can match another available source.
jwst_i, ground_i, sep = closest_unique_matches(jwst_coord, ground_coord, RADIUS)
report_multiple_ground_candidates(jwst_coord, ground_coord, RADIUS,
                                  jwst['ID'], ground['ID'], jwst_i, ground_i)

n_dup = len(ground_i)
if n_dup != EXPECTED_DUPS:
    print(f'WARNING: found {n_dup} duplicates, expected {EXPECTED_DUPS}')

# Save each accepted pair side by side before converting types for the combined table.
duplicates = hstack([ground[ground_i], jwst[jwst_i]],
                    table_names=['ground', 'JWST'],
                    uniq_col_name='{col_name}_{table_name}')
duplicates['separation_arcsec'] = sep
duplicates.sort('separation_arcsec')
duplicates_path = Path(__file__).resolve().parent / 'JWST.Ground_duplicates.fits'
duplicates.write(duplicates_path, overwrite=True)
print(f'Wrote {n_dup} matches (r <= {RADIUS}") to {duplicates_path}')

# Drop the ground-based copy of each duplicate and keep the JWST row
keep = np.ones(len(ground), dtype=bool)
keep[ground_i] = False

# Make column dtypes consistent between the tables (Q is a string in JWST, int in ground)
for tbl in (jwst, ground):
    tbl['Q'] = tbl['Q'].astype(str)
    for col in ('Type', 'Source'):
        tbl[col] = tbl[col].astype(str)

# JWST rows plus the ground rows that are not duplicates
combined = vstack([jwst, ground[keep]], metadata_conflicts='silent')

# Flag the closest combined source to the BCG only if it is within the radius.
bcg = Table.read(specz_dir.parent / 'BCG' / 'MOO_1142_bcg_info.txt', format='ascii.csv')
bcg_coord = SkyCoord(bcg['ra_spz'][0] * u.deg, bcg['dec_spz'][0] * u.deg, frame='icrs')
combined_coord = SkyCoord(combined['RA'] * u.deg, combined['Dec'] * u.deg, frame='icrs')
combined['BCG'] = np.zeros(len(combined), dtype=bool)
if len(combined):
    bcg_sep = combined_coord.separation(bcg_coord).to(u.arcsec).value
    closest = np.argmin(bcg_sep)
    if bcg_sep[closest] <= BCG_RADIUS:
        combined['BCG'][closest] = True
        print(f'BCG: {combined["ID"][closest]} (sep = {bcg_sep[closest]:.3f}")')
    else:
        print(f'WARNING: no source within {BCG_RADIUS}" of the BCG '
              f'(closest separation = {bcg_sep[closest]:.3f}")')
else:
    print('WARNING: cannot identify the BCG in an empty combined catalog')

combined.sort('Spec-z')

today = date.today().isoformat()
out_name = f'{today}.MOO_1142+1527.speczs.all.fits'
combined.write(specz_dir / out_name, overwrite=True)
print(f'{n_dup} duplicates; wrote {len(combined)} sources to {specz_dir / out_name}')


def is_quality_one(value):
    # Q also contains categorical values such as A/B, which must be retained.
    try:
        return float(value) == 1.0
    except ValueError:
        return False


good_q = combined[[not is_quality_one(q) for q in combined['Q']]]
good_q_name = f'{today}.MOO_1142+1527.speczs.all.goodQ.fits'
good_q.write(specz_dir / good_q_name, overwrite=True)
print(f'Wrote {len(good_q)} sources (Q != 1.0) to {specz_dir / good_q_name}')
