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
# -------------------------------------------------------------------------------------------- #
specz_dir = Path(__file__).resolve().parent.parent
jwst_path = specz_dir / '2026-10-15.MOO_1142+1527.speczs.JWST.fits'
ground_path = specz_dir / '2026-10-19.MOO_1142+1527.speczs.internal.groundBasedAll.fits'
pdf_path = Path(__file__).resolve().parent / 'cross_match_ground_w_JWST.pdf'

jwst = Table.read(jwst_path)
ground = Table.read(ground_path)

jwst_coord = SkyCoord(jwst['RA'] * u.deg, jwst['Dec'] * u.deg, frame='icrs')
ground_coord = SkyCoord(ground['RA'] * u.deg, ground['Dec'] * u.deg, frame='icrs')

# Nearest ground-based neighbor for each JWST source
idx, sep2d, _ = jwst_coord.match_to_catalog_sky(ground_coord)
sep = sep2d.to(u.arcsec).value

# One-to-one: for each ground source keep only the closest JWST source. This is
# radius independent, since the closest JWST source is within r whenever any is.

# Sort JWST sources by matched ground index, ties broken by smallest separation,
# so each group sharing a ground source starts with its closest JWST source.
order = np.lexsort((sep, idx))
# Flag the first entry of each group (where the ground index changes) as the winner.
first = np.ones(len(order), dtype=bool)
first[1:] = idx[order][1:] != idx[order][:-1]
# Map winners back to the original JWST ordering.
unique = np.zeros(len(idx), dtype=bool)
unique[order[first]] = True
# Losers get an infinite separation so they never pass the sep <= r cut.
sep = np.where(unique, sep, np.inf)

radii = np.round(np.arange(0.1, 2.5 + 0.05, 0.1), 1)
n_match = np.zeros(len(radii), dtype=int)
std_sep = np.full(len(radii), np.nan)
for i, r in enumerate(radii):
    m = sep <= r
    n_match[i] = m.sum()
    if n_match[i] > 1:
        std_sep[i] = np.std(sep[m], ddof=1)
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
