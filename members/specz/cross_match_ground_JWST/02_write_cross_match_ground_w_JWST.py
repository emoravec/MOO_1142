"""
Write closest-first one-to-one ground/JWST cross-matches at 1", 2" and 2.5".
"""
from pathlib import Path

from astropy.coordinates import SkyCoord
from astropy.table import Table, hstack
from astropy import units as u
from cross_match_utils import closest_unique_matches, report_multiple_ground_candidates
# -------------------------------------------------------------------------------------------- #
specz_dir = Path(__file__).resolve().parent.parent
out_dir = Path(__file__).resolve().parent
jwst = Table.read(specz_dir / '2026-10-15.MOO_1142+1527.speczs.JWST.fits')
#ground = Table.read(specz_dir / '2026-10-19.MOO_1142+1527.speczs.internal.groundBasedAll.fits')
ground = Table.read(specz_dir / '2026-11.MOO_1142+1527.speczs.allGroundBased.fits')

jwst_coord = SkyCoord(jwst['RA'] * u.deg, jwst['Dec'] * u.deg, frame='icrs')
ground_coord = SkyCoord(ground['RA'] * u.deg, ground['Dec'] * u.deg, frame='icrs')

bcg = Table.read(specz_dir.parent / 'BCG' / 'MOO_1142_bcg_info.txt', format='ascii.csv')
bcg_coord = SkyCoord(bcg['ra_spz'][0] * u.deg, bcg['dec_spz'][0] * u.deg, frame='icrs')

for radius, name in [(1.0, 'ground_JWST_1as.fits'), (2.0, 'ground_JWST_2as.fits'), (2.5, 'ground_JWST_2.5as.fits')]:
    # Search all candidates within this radius, accepting closest unused pairs first.
    jwst_i, ground_i, sep = closest_unique_matches(jwst_coord, ground_coord, radius)
    report_multiple_ground_candidates(jwst_coord, ground_coord, radius,
                                      jwst['ID'], ground['ID'], jwst_i, ground_i)

    # Place the matched source rows side by side, distinguishing shared column names.
    tbl = hstack([ground[ground_i], jwst[jwst_i]],
                 table_names=['ground', 'JWST'],
                 uniq_col_name='{col_name}_{table_name}')
    tbl['separation_arcsec'] = sep
    tbl.rename_column('ID_ground', 'ID_ground')
    tbl.rename_column('ID_JWST', 'ID_JWST')
    tbl['dz'] = tbl['Spec-z_JWST'] - tbl['Spec-z_ground']
    # BCG separation measured from the ground-based position of each match
    tbl['sep_BCG'] = ground_coord[ground_i].separation(bcg_coord).to(u.arcsec).value
    tbl.sort('separation_arcsec')

    tbl.write(out_dir / name, overwrite=True)
    print(f'Wrote {len(tbl)} matches (r <= {radius}") to {out_dir / name}')
