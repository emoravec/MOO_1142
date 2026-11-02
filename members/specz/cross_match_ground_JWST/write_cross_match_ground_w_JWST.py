"""
Write the one-to-one ground-based/JWST cross-matches at 1" and 2" radii to FITS tables.
"""
import numpy as np
from pathlib import Path

from astropy.coordinates import SkyCoord
from astropy.table import Table, hstack
from astropy import units as u
# -------------------------------------------------------------------------------------------- #
specz_dir = Path(__file__).resolve().parent.parent
out_dir = Path(__file__).resolve().parent
jwst = Table.read(specz_dir / '2026-10-15.MOO_1142+1527.speczs.JWST.fits')
ground = Table.read(specz_dir / '2026-10-19.MOO_1142+1527.speczs.internal.groundBasedAll.fits')

jwst_coord = SkyCoord(jwst['RA'] * u.deg, jwst['Dec'] * u.deg, frame='icrs')
ground_coord = SkyCoord(ground['RA'] * u.deg, ground['Dec'] * u.deg, frame='icrs')

bcg = Table.read(specz_dir.parent / 'MOO_1142_bcg_info.txt', format='ascii.csv')
bcg_coord = SkyCoord(bcg['ra_spz'][0] * u.deg, bcg['dec_spz'][0] * u.deg, frame='icrs')

# Nearest ground-based neighbor for each JWST source
idx, sep2d, _ = jwst_coord.match_to_catalog_sky(ground_coord)
sep = sep2d.to(u.arcsec).value

# One-to-one: for each ground source keep only the closest JWST source
order = np.lexsort((sep, idx))
first = np.ones(len(order), dtype=bool)
first[1:] = idx[order][1:] != idx[order][:-1]
winner = np.zeros(len(idx), dtype=bool)
winner[order[first]] = True

for radius, name in [(1.0, 'ground_JWST_1as.fits'), (2.0, 'ground_JWST_2as.fits'), (2.5, 'ground_JWST_2.5as.fits')]:
    m = winner & (sep <= radius)
    jwst_i = np.where(m)[0]
    ground_i = idx[m]

    tbl = hstack([ground[ground_i], jwst[jwst_i]],
                 table_names=['ground', 'JWST'],
                 uniq_col_name='{col_name}_{table_name}')
    tbl['separation_arcsec'] = sep[m]
    tbl.rename_column('ID_ground', 'ID')
    tbl.remove_column('ID_JWST')
    tbl['dz'] = tbl['Spec-z_JWST'] - tbl['Spec-z_ground']
    # BCG separation measured from the ground-based position of each match
    tbl['sep_BCG'] = ground_coord[ground_i].separation(bcg_coord).to(u.arcsec).value
    tbl.sort('separation_arcsec')

    tbl.write(out_dir / name, overwrite=True)
    print(f'Wrote {len(tbl)} matches (r <= {radius}") to {out_dir / name}')
