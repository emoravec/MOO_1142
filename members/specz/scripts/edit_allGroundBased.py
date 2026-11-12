"""
Edit the ground-based spec-z FITS table: add IDs and source-dependent errors.
"""
from pathlib import Path

from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy import units as u
# -------------------------------------------------------------------------------------------- #
location = Path(__file__).resolve().parent.parent
in_path = location / 'extra' / '2026-10-19.MOO_1142+1527.speczs.internal.groundBasedAll.fits'
out_path = location / '2026-11.MOO_1142+1527.speczs.internal.groundBasedAll.fits'

tbl = Table.read(in_path)

for col in ('Cluster', 'Comment'):
    if col in tbl.colnames:
        tbl.remove_column(col)

error_by_type = {'A': 0.003, 'E': 0.0005, 'A+E': 0.0005, 'star': 0.0}
unknown_types = {row['Type'] for row in tbl
                 if row['Source'] not in ('PN', 'PD')
                 and row['Type'] not in error_by_type}
if unknown_types:
    raise ValueError(f"Unknown Type values: {sorted(unknown_types)}")

errors = []
for index, row in enumerate(tbl):
    if row['Source'] in ('PN', 'PD'):
        if row['Q'] == 3 or (row['Q'] == 2 and row['Type'] in ('E', 'A', 'A+E')):
            errors.append(0.002)
        else:
            print("you have a Palomar source that doesn't fit the standard boiler plate error estimation."
                  f" Row {index}: Source={row['Source']}, Q={row['Q']}, Type={row['Type']}")
            errors.append(float('nan'))
    else:
        errors.append(error_by_type[row['Type']])
specz_index = tbl.colnames.index('Spec-z')
tbl.add_column(errors, name='Spec-z_err_hi', index=specz_index + 1)
tbl.add_column(errors, name='Spec-z_err_lo', index=specz_index + 2)

coord = SkyCoord(ra=tbl['RA'], dec=tbl['Dec'], unit=(u.deg, u.deg), frame='icrs')
RA_hms = coord.ra.to_string(unit=u.hourangle, sep=':', precision=1, pad=True)
Dec_dms = coord.dec.to_string(unit=u.deg, sep=':', precision=0, alwayssign=True, pad=True)

# ID of the form JHHMMSS.SpDDMMSS
ID = ["J" + r[0:2] + r[3:5] + r[6:10] + "p" + d[1:3] + d[4:6] + d[7:9]
      for r, d in zip(RA_hms, Dec_dms)]

tbl.add_column(ID, name='ID', index=0)

tbl.write(out_path, overwrite=True)
print(f"Wrote {len(tbl)} rows to {out_path}")