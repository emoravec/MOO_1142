"""
Cross-match ground-based and wide-field survey spectroscopic redshift catalogs.

Reads two FITS catalogs, performs a spatial cross-match using SkyCoord,
and outputs matched sources with separation information.
"""

from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.table import Table, hstack
import astropy.units as u
import numpy as np


def cross_match_catalogs(ground_file, wide_file, output_file, 
                         coord_cols=('RA', 'Dec'), sep_unit='arcsec'):
    """
    Cross-match two spectroscopic redshift catalogs spatially.
    
    Parameters
    ----------
    ground_file : str
        Path to ground-based catalog FITS file
    wide_file : str
        Path to wide-field surveys catalog FITS file
    output_file : str
        Path for output FITS file with matches
    coord_cols : tuple
        Names of (RA, Dec) columns in the FITS files
    sep_unit : str
        Unit for separation ('arcsec', 'arcmin', or 'deg')
    """
    
    # Read the catalogs
    print(f"Reading ground-based catalog: {ground_file}")
    with fits.open(ground_file) as hdul:
        ground_data = hdul[1].data
        ground_table = Table(ground_data)
    
    print(f"Reading wide-field surveys catalog: {wide_file}")
    with fits.open(wide_file) as hdul:
        wide_data = hdul[1].data
        wide_table = Table(wide_data)
    
    print(f"Ground-based sources: {len(ground_table)}")
    print(f"Wide-field sources: {len(wide_table)}")
    
    # Create SkyCoord objects
    ra_col, dec_col = coord_cols
    ground_coords = SkyCoord(ra=ground_table[ra_col], 
                             dec=ground_table[dec_col], 
                             unit=u.deg)
    wide_coords = SkyCoord(ra=wide_table[ra_col], 
                           dec=wide_table[dec_col], 
                           unit=u.deg)
    
    # Cross-match: find closest ground-based match for each wide-field source
    idx, sep2d, _ = wide_coords.match_to_catalog_sky(ground_coords)
    
    # Convert separation to requested unit
    if sep_unit == 'arcsec':
        separation = sep2d.to(u.arcsec)
    elif sep_unit == 'arcmin':
        separation = sep2d.to(u.arcmin)
    elif sep_unit == 'deg':
        separation = sep2d.to(u.deg)
    else:
        raise ValueError(f"Unknown separation unit: {sep_unit}")
    
    # Create output table with all columns from both catalogs
    output_table = hstack([wide_table, ground_table[idx]], 
                          table_names=['wide', 'ground'])
    
    # Add separation column
    sep_col_name = f'separation_{sep_unit}'
    output_table[sep_col_name] = separation.value
    
    # Write to FITS
    print(f"\nWriting {len(output_table)} matched sources to: {output_file}")
    print(f"Separation column: {sep_col_name}")
    output_table.write(output_file, format='fits', overwrite=True)
    
    # Print statistics
    print(f"\nCross-match Statistics:")
    print(f"  Minimum separation: {separation.min():.4f}")
    print(f"  Maximum separation: {separation.max():.4f}")
    print(f"  Mean separation: {separation.mean():.4f}")
    print(f"  Median separation: {np.median(separation.value):.4f}")


if __name__ == '__main__':
    ground_file = '/Users/emoravec/Documents/Research/MOO_1142/members/specz/2026-11.MOO_1142+1527.speczs.internal.groundBasedAll.fits'
    wide_file = '/Users/emoravec/Documents/Research/MOO_1142/members/specz/2026-11.MOO_1142+1527.speczs.wide_field_surveys.fits'
    output_file = '/Users/emoravec/Documents/Research/MOO_1142/members/specz/extra/ground_based_matches.fits'
    
    cross_match_catalogs(ground_file, wide_file, output_file)
