"""
Cross match the JWST spec-z catalog with the ground-based spec-z catalog over a
range of matching radii (2.5" down to 0.1" in 0.1" steps). Left panel: number of
matches vs. radius. Right panel: standard deviation of the matched-pair
separations vs. radius.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from pathlib import Path

from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy import units as u
from cross_match_utils import closest_unique_matches
# -------------------------------------------------------------------------------------------- #
specz_dir = Path(__file__).resolve().parent.parent
jwst_path = specz_dir / '2026-10-15.MOO_1142+1527.speczs.JWST.fits'
# ground_path = specz_dir / '2026-10-19.MOO_1142+1527.speczs.internal.groundBasedAll.fits'
ground_path = specz_dir / '2026-11.MOO_1142+1527.speczs.allGroundBased.fits'
pdf_path = Path(__file__).resolve().parent / 'cross_match_ground_w_JWST.pdf'

jwst = Table.read(jwst_path)
ground = Table.read(ground_path)

jwst_coord = SkyCoord(jwst['RA'] * u.deg, jwst['Dec'] * u.deg, frame='icrs')
ground_coord = SkyCoord(ground['RA'] * u.deg, ground['Dec'] * u.deg, frame='icrs')

radii = np.round(np.arange(0.1, 2.5 + 0.05, 0.1), 1)
n_match = np.zeros(len(radii), dtype=int)
std_sep = np.full(len(radii), np.nan)
for i, r in enumerate(radii):
    # Search all pairs within this radius and accept closest unused sources first.
    _, _, sep = closest_unique_matches(jwst_coord, ground_coord, r)
    n_match[i] = len(sep)
    # The sample standard deviation requires at least two accepted pairs.
    if n_match[i] > 1:
        std_sep[i] = np.std(sep, ddof=1)
    print(f'r = {r:.1f}": {n_match[i]} matches, std = {std_sep[i]:.3f}"')

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))
ax1.plot(radii, n_match, 'o-')
ax1.set_xlabel('Match radius [arcsec]')
ax1.set_ylabel('Number of matches')
ax1.yaxis.set_major_locator(MultipleLocator(1))
ax1.xaxis.set_major_locator(MultipleLocator(0.5))
ax1.xaxis.set_minor_locator(MultipleLocator(0.1))

ax2.plot(radii, std_sep, 'o-')
ax2.set_xlabel('Match radius [arcsec]')
ax2.set_ylabel('Std. [arcsec]')
ax2.yaxis.set_major_locator(MultipleLocator(0.1))
ax2.xaxis.set_major_locator(MultipleLocator(0.5))
ax2.xaxis.set_minor_locator(MultipleLocator(0.1))
ax2.grid(True)

fig.tight_layout()
fig.savefig(pdf_path)
print(f'Saved {pdf_path}')
