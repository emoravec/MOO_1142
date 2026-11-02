"""
Cross-match the ground-based and JWST spectroscopic redshift tables (1.5" radius, one-to-one,
closest pairs first) and write the combined table. For duplicates the JWST row is kept.
"""
import numpy as np
from datetime import date
from pathlib import Path

from astropy.coordinates import SkyCoord
from astropy.table import Table, vstack
from astropy import units as u
# -------------------------------------------------------------------------------------------- #
RADIUS = 1.5  # arcsec
EXPECTED_DUPS = 8

specz_dir = Path(__file__).resolve().parent.parent
jwst = Table.read(specz_dir / '2026-10-15.MOO_1142+1527.speczs.JWST.fits')
ground = Table.read(specz_dir / '2026-10-19.MOO_1142+1527.speczs.internal.groundBasedAll.fits')

jwst_coord = SkyCoord(jwst['RA'] * u.deg, jwst['Dec'] * u.deg, frame='icrs')
ground_coord = SkyCoord(ground['RA'] * u.deg, ground['Dec'] * u.deg, frame='icrs')

# One-to-one matching:
# 1) search_around_sky returns EVERY (JWST, ground) pair within the radius. One source can
#    appear in several pairs (e.g. two ground sources both within 1.5" of one JWST source).
# 2) Sort the pairs by separation and walk through them from the closest to the farthest.
# 3) Accept a pair only if neither its JWST source nor its ground source is already in a
#    previously accepted pair; otherwise skip it. Because closer pairs are handled first,
#    each source ends up matched to its closest available partner and to only one source.
j_idx, g_idx, sep2d, _ = ground_coord.search_around_sky(jwst_coord, RADIUS * u.arcsec)
used_j, used_g = set(), set()
ground_dup = []  # indices of the ground rows that duplicate a JWST row
for k in np.argsort(sep2d.arcsec):
    j, g = j_idx[k], g_idx[k]
    if j in used_j or g in used_g:
        print(f'Skipping pair JWST {jwst["ID"][j]} / ground {ground["ID"][g]} '
              f'(sep = {sep2d.arcsec[k]:.2f}"): '
              f'{"JWST" if j in used_j else "ground"} source already matched')
        continue
    used_j.add(j)
    used_g.add(g)
    ground_dup.append(g)

n_dup = len(ground_dup)
if n_dup != EXPECTED_DUPS:
    print(f'WARNING: found {n_dup} duplicates, expected {EXPECTED_DUPS}')

# Drop the ground-based copy of each duplicate and keep the JWST row
keep = np.ones(len(ground), dtype=bool)
keep[ground_dup] = False

# Make column dtypes consistent between the tables (Q is a string in JWST, int in ground)
for tbl in (jwst, ground):
    tbl['Q'] = tbl['Q'].astype(str)
    for col in ('Type', 'Source'):
        tbl[col] = tbl[col].astype(str)

# JWST rows plus the ground rows that are not duplicates
combined = vstack([jwst, ground[keep]], metadata_conflicts='silent')

# Separation of each source from the BCG (arcsec), added as the last column
bcg = Table.read(specz_dir.parent / 'MOO_1142_bcg_info.txt', format='ascii.csv')
bcg_coord = SkyCoord(bcg['ra_spz'][0] * u.deg, bcg['dec_spz'][0] * u.deg, frame='icrs')
combined_coord = SkyCoord(combined['RA'] * u.deg, combined['Dec'] * u.deg, frame='icrs')
combined['sep_BCG'] = combined_coord.separation(bcg_coord).to(u.arcsec).value

combined.sort('Spec-z')

out_name = f'{date.today().isoformat()}.MOO_1142+1527.speczs.all.fits'
combined.write(specz_dir / out_name, overwrite=True)
print(f'{n_dup} duplicates; wrote {len(combined)} sources to {specz_dir / out_name}')
