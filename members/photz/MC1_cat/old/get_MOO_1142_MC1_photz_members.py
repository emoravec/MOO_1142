"""
Get MOO 1142 MC1 photometric redshift members.

Created in 2025-10
"""
import numpy as np
from pathlib import Path
from astropy.table import Table

output_dir = Path(__file__).resolve().parent

# Read in the photometric redshift catalog
photz_catalog = Table.read('/Users/emoravec/Documents/Research/MaDCoWS/catalogs/photz/MaDCoWS_all_photz_R2_2025-08-22.fits')

print(f"Total sources in photz catalog: {len(photz_catalog)}")

# Filter for MOO_1142+1527 cluster with INT_PZ >= 0.3
moo_1142_members = photz_catalog[(photz_catalog['Cluster_name'] == 'MOO_1142+1527') & 
                                 (photz_catalog['INT_PZ'] >= 0.3)]

print(f"MOO_1142+1527 sources with INT_PZ >= 0.3: {len(moo_1142_members)}")
print(f"INT_PZ range: {moo_1142_members['INT_PZ'].min():.3f} - {moo_1142_members['INT_PZ'].max():.3f}")

# Save the filtered members to the current directory
output_path = output_dir / 'MOO_1142_R2_photz_members_IntPzgt03.fits'
moo_1142_members.write(output_path, overwrite=True)
print(f"Saved {len(moo_1142_members)} members to: {output_path}")
