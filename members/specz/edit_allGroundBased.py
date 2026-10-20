"""
Edit the ground-based spec-z FITS table: drop 'Cluster' and add an 'ID' column.
"""
from pathlib import Path

from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy import units as u
# -------------------------------------------------------------------------------------------- #
location = Path(__file__).resolve().parent
name = '2026-10-19.MOO_1142+1527.speczs.internal.groundBasedAll.fits'
in_path = location / 'extra' / name
out_path = location / name

tbl = Table.read(in_path)

for col in ('Cluster', 'Comment'):
    if col in tbl.colnames:
        tbl.remove_column(col)

coord = SkyCoord(ra=tbl['RA'], dec=tbl['Dec'], unit=(u.deg, u.deg), frame='icrs')
RA_hms = coord.ra.to_string(unit=u.hourangle, sep=':', precision=1, pad=True)
Dec_dms = coord.dec.to_string(unit=u.deg, sep=':', precision=0, alwayssign=True, pad=True)

# ID of the form JHHMMSS.SpDDMMSS
ID = ["J" + r[0:2] + r[3:5] + r[6:10] + "p" + d[1:3] + d[4:6] + d[7:9]
      for r, d in zip(RA_hms, Dec_dms)]

tbl.add_column(ID, name='ID', index=0)

tbl.write(out_path, overwrite=True)
print(f"Wrote {len(tbl)} rows to {out_path}")
