"""One-to-one matching of JWST and ground-based coordinates."""
import numpy as np
from astropy import units as u


def closest_unique_matches(jwst_coord, ground_coord, radius):
    """Return JWST indices, ground indices and separations in arcsec.

    Consider all pairs within the radius, accepting closest pairs first when
    neither source has already been used. Equal separations are resolved by
    JWST row index, then ground row index.
    """
    # search_around_sky returns indices into its argument first, then into self.
    # These arrays contain every candidate pair within the requested radius.
    ground_i, jwst_i, sep2d, _ = jwst_coord.search_around_sky(
        ground_coord, radius * u.arcsec
    )
    sep = sep2d.to_value(u.arcsec)
    # Process shortest separations first; row indices deterministically break ties.
    order = np.lexsort((ground_i, jwst_i, sep))
    used_jwst, used_ground = set(), set()
    accepted = []
    for pair in order:
        j, g = jwst_i[pair], ground_i[pair]
        # A rejected source can still match a later, more distant available source.
        if j not in used_jwst and g not in used_ground:
            used_jwst.add(j)
            used_ground.add(g)
            accepted.append(pair)

    accepted = np.asarray(accepted, dtype=int)
    return jwst_i[accepted], ground_i[accepted], sep[accepted]


def report_multiple_ground_candidates(jwst_coord, ground_coord, radius,
                                      jwst_ids, ground_ids, matched_jwst, matched_ground):
    """Print all candidates for JWST sources with multiple ground neighbors."""
    ground_i, jwst_i, sep2d, _ = jwst_coord.search_around_sky(
        ground_coord, radius * u.arcsec
    )
    sep = sep2d.to_value(u.arcsec)
    ground_matches = dict(zip(matched_ground, matched_jwst))
    jwst_matches = dict(zip(matched_jwst, matched_ground))
    sources, counts = np.unique(jwst_i, return_counts=True)

    for j in sources[counts > 1]:
        print(f'JWST {jwst_ids[j]} has multiple ground candidates within {radius}":')
        # Show candidates in separation order and identify the final assignment.
        candidates = np.flatnonzero(jwst_i == j)
        candidates = candidates[np.lexsort((ground_i[candidates], sep[candidates]))]
        for pair in candidates:
            g = ground_i[pair]
            if jwst_matches.get(j) == g:
                status = 'accepted'
            elif g in ground_matches:
                other_j = ground_matches[g]
                other_sep = jwst_coord[other_j].separation(ground_coord[g]).arcsec
                status = (f'ground matched to JWST {jwst_ids[other_j]} '
                          f'(sep = {other_sep:.3f}")')
            else:
                status = 'not selected'
            print(f'  JWST {jwst_ids[j]} / ground {ground_ids[g]} '
                  f'(sep = {sep[pair]:.3f}"): {status}')
        if j not in jwst_matches:
            print(f'  JWST {jwst_ids[j]} remains unmatched')
