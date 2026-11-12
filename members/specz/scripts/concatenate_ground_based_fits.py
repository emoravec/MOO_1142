"""
Concatenate two FITS tables: internal ground-based and wide-field surveys.
"""
from pathlib import Path
import numpy as np

from astropy.table import Table, vstack

# Setup paths
input_path1 = Path('/Users/emoravec/Documents/Research/MOO_1142/members/specz/2026-11.MOO_1142+1527.speczs.internal.groundBasedAll.fits')
input_path2 = Path('/Users/emoravec/Documents/Research/MOO_1142/members/specz/2026-11.MOO_1142+1527.speczs.wide_field_surveys.fits')
output_path = Path('/Users/emoravec/Documents/Research/MOO_1142/members/specz/2026-11.MOO_1142+1527.speczs.allGroundBased.fits')

# Read the two FITS tables
tbl1 = Table.read(input_path1)
tbl2 = Table.read(input_path2)

print(f"Table 1 ({input_path1.name}): {len(tbl1)} rows, columns: {tbl1.colnames}")
print(f"Table 2 ({input_path2.name}): {len(tbl2)} rows, columns: {tbl2.colnames}")

# Standardize column types and naming
# Convert Q column to float to accommodate NaN values from wide_field_surveys
tbl1['Q'] = tbl1['Q'].astype(np.float64)
tbl2['Q'] = tbl2['Q'].astype(np.float64)

# Rename 'Facility' to 'Source' in table 2 for consistency
if 'Facility' in tbl2.colnames:
    tbl2.rename_column('Facility', 'Source')

# Add missing columns as needed
if 'Source' not in tbl1.colnames:
    tbl1['Source'] = [''] * len(tbl1)
if 'Source' not in tbl2.colnames:
    tbl2['Source'] = [''] * len(tbl2)

# Ensure column order is the same - keep columns that are in both tables
common_cols = [col for col in tbl1.colnames if col in tbl2.colnames]
tbl1 = tbl1[common_cols]
tbl2 = tbl2[common_cols]

# Concatenate the tables
combined = vstack([tbl1, tbl2])

print(f"Combined table: {len(combined)} rows, columns: {combined.colnames}")

# Write to output FITS file
combined.write(output_path, overwrite=True)
print(f"Wrote combined table to {output_path}")
