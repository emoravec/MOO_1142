"""
Convert the JWST spec-z csv into a FITS table.
"""
import numpy as np
import pandas as pd
from pathlib import Path

from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy import units as u
# -------------------------------------------------------------------------------------------- #
location = Path(__file__).resolve().parent.parent
csv_path = location / '2026-10-15.MOO_1142+1527.speczs.JWST.csv'
fits_path = location / '2026-10-15.MOO_1142+1527.speczs.JWST.fits'

df = pd.read_csv(csv_path)
RA_hms = df['RA_hms'].astype(str).str.strip()
Dec_dms = df['Dec_dms'].astype(str).str.strip()

# ID of the form JHHMMSS.SpDDMMSS
ID = ("J" + RA_hms.str[0:2] + RA_hms.str[3:5] + RA_hms.str[6:10]
      + "p" + Dec_dms.str[1:3] + Dec_dms.str[4:6] + Dec_dms.str[7:9])

coord = SkyCoord(
    ra=RA_hms.tolist(),
    dec=Dec_dms.tolist(),
    unit=(u.hourangle, u.deg),
    frame="icrs",
)

# Use z_grism when present, otherwise z_prism
use_grism = df['z_grism'].notna().to_numpy()
specz = np.where(use_grism, df['z_grism'], df['z_prism'])
specz_err_hi = np.where(use_grism, df['z_grism_err_hi'], df['z_prism_err_hi'])
specz_err_lo = np.where(use_grism, df['z_grism_err_lo'], df['z_prism_err_lo'])
source = np.where(use_grism, 'JG', 'JP')

n = len(df)
tbl = Table()
tbl['ID'] = ID.to_numpy().astype(str)
tbl['RA'] = coord.ra.deg
tbl['Dec'] = coord.dec.deg
tbl['Spec-z'] = specz.astype(float)
tbl['Spec-z_err_hi'] = specz_err_hi.astype(float)
tbl['Spec-z_err_lo'] = specz_err_lo.astype(float)
tbl['Source'] = source
tbl['Q'] = np.full(n, 'A/B')
tbl['Type'] = np.full(n, 'NA')
tbl['Date'] = np.full(n, '2025-06-07')

tbl.write(fits_path, overwrite=True)
print(f"Wrote {n} rows to {fits_path}")
