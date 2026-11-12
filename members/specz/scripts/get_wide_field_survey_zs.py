"""
Extract wide-field survey (SD and DS facility) spec-z data for MOO_1142+1527 and write to FITS.
"""
from pathlib import Path

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy import units as u

# Setup paths
script_location = Path(__file__).resolve().parent
input_path = Path('/Users/emoravec/Documents/Research/MaDCoWS/catalogs/specz/Specz_data_product_sorted_internalfirst.txt')
output_path = Path('/Users/emoravec/Documents/Research/MOO_1142/members/specz/2026-11.MOO_1142+1527.speczs.wide_field_surveys.fits')

# Read input data
tbl = Table.read(
    input_path,
    format='ascii.no_header',
    names=['Cluster_ID', 'RA', 'Dec', 'Spec-z', 'Q', 'Type', 'Facility', 'Date']
)

# Filter for MOO_1142+1527 and facilities SD or DS
mask = (tbl['Cluster_ID'] == 'MOO_1142+1527') & (tbl['Facility'].astype(str) != '-') & \
       ((tbl['Facility'] == 'SD') | (tbl['Facility'] == 'DS'))
tbl = tbl[mask]

if len(tbl) == 0:
    print("No rows match the filter criteria!")
    exit()

# Remove Cluster_ID column
tbl.remove_column('Cluster_ID')

# Convert Q column to float (convert '-' to np.nan)
Q_values = []
for q in tbl['Q']:
    if q == '-' or str(q).strip() == '-':
        Q_values.append(np.nan)
    else:
        try:
            Q_values.append(float(q))
        except (ValueError, TypeError):
            Q_values.append(np.nan)
tbl['Q'] = np.array(Q_values, dtype=np.float64)

# Create ID column from RA/Dec coordinates
coord = SkyCoord(ra=tbl['RA'], dec=tbl['Dec'], unit=(u.deg, u.deg), frame='icrs')
RA_hms = coord.ra.to_string(unit=u.hourangle, sep=':', precision=1, pad=True)
Dec_dms = coord.dec.to_string(unit=u.deg, sep=':', precision=0, alwayssign=True, pad=True)

# ID of the form JHHMMSS.SpDDMMSS
ID = ["J" + r[0:2] + r[3:5] + r[6:10] + "p" + d[1:3] + d[4:6] + d[7:9]
      for r, d in zip(RA_hms, Dec_dms)]

tbl.add_column(ID, name='ID', index=0)

# Add error columns with value 0, right after Spec-z column
specz_index = tbl.colnames.index('Spec-z')
errors = np.zeros(len(tbl))
tbl.add_column(errors, name='Spec-z_err_hi', index=specz_index + 1)
tbl.add_column(errors, name='Spec-z_err_lo', index=specz_index + 2)

# Write to FITS
tbl.write(output_path, overwrite=True)
print(f"Wrote {len(tbl)} rows to {output_path}")
