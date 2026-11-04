"""
Finds MOO 1142 BCG radio and spitzer information from radio source catalog.

Created in 2023-011
"""
import numpy as np

from astropy.coordinates import SkyCoord
from astropy.table import hstack,Table,vstack
from astropy import units as u
# -------------------------------------------------------------------------------------------- #
location = '/Users/emoravec/Documents/Research/merging_clusters/'
sources = Table.read(location + 'catalogs/vla21b/vla_spz_LS/MOO_1142+1527_vla_spz_LS_matches.fits')
bcg = SkyCoord('11:42:47.434 15:27:11.29',unit=(u.hourangle,u.deg))
member_cat = SkyCoord(sources['Radio_RA'], sources['Radio_Dec'], unit='degree')
idx, d2d, d3d = bcg.match_to_catalog_sky(member_cat)

max_sep = 1.0 * u.arcsec
sep_constraint = d2d <= max_sep
if sep_constraint == True:
    bcg_match = sources[idx]

output_t = Table(bcg_match)
bcg_match.write(location + 'analysis/MOO_1142/members/MOO_1142_bcg_info.txt',format='ascii.csv')
