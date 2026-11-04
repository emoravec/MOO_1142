"""
Extract the MOO_1142+1527 BCG row from the MC1 BCG catalog.
"""
from pathlib import Path

from astropy.table import Table

catalog_path = Path('/Users/emoravec/Documents/Research/MaDCoWS/catalogs/photz/BCG/BCGs.fits')
output_path = Path(__file__).resolve().parent / 'MOO_1142_MC1_BCG.fits'
cluster_name = 'MOO_1142+1527'

catalog = Table.read(catalog_path)
matches = catalog[catalog['Cluster_name'] == cluster_name]

if len(matches) != 1:
    raise ValueError(
        f'Expected exactly one row for {cluster_name!r}, found {len(matches)}'
    )

matches.write(output_path, overwrite=True)
print(f'Wrote {len(matches)} row to {output_path}')
