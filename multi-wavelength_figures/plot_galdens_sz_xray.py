"""
Create a 2x2 figure (shared axes) of MOO J1142+1527 multi-wavelength maps.

Top left: INT_PZ-weighted galaxy surface density with photz sources and the BCG.
Top right: XMM X-ray image with MUSTANG-2 SZ contours.
Bottom left: galaxy density with INT_PZ >= 0.3 sources, X-ray (black) and SZ (red) contours.
Bottom right: same without sources, with the X-ray main and west subcluster centers marked.
"""
from pathlib import Path
import warnings

import aplpy
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, PowerNorm
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy.utils.exceptions import AstropyDeprecationWarning
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from reproject import reproject_interp
from scipy.ndimage import gaussian_filter

try:
    from palettable.cartocolors.sequential import Mint_7 as pal4
    cmap = pal4.mpl_colormap.copy()
except ImportError:
    cmap = plt.get_cmap('YlGn').copy()

plt.rc('font', **{'family': 'serif', 'size': 11})
cmap.set_under('white')
# -------------------------------------------------------------------------------------------- #
warnings.filterwarnings(
    'ignore',
    category=AstropyDeprecationWarning,
    module='astropy.visualization.wcsaxes.coordinate_helpers',
)
warnings.filterwarnings(
    'ignore',
    message='You passed a edgecolor/edgecolors.*for an unfilled marker',
    module='aplpy.*',
)
# -------------------------------------------------------------------------------------------- #
BASE = Path('/Users/emoravec/Documents/Research/MOO_1142')
PHOTZ_PATH = BASE / 'members/photz/MC1_cat/catalogs/MOO_1142+1527.fits'
DENSITY_PATH = BASE / 'members/photz/MC1_cat/galdens/MOO_1142+1527_final_full_IntPzdens.fits'
BCG_PATH = BASE / 'members/BCG/MOO_1142_MC1_BCG.fits'
OUTPUT_PATH = BASE / 'multi-wavelength_figures/pdfs/MOO_1142_galdens_sz_xray.pdf'

TOP_RADIUS = 2.5 * u.arcmin
BOTTOM_RADIUS = 1 * u.arcmin
SCALE_BAR = 30 * u.arcsec
DENSITY_STEP = 0.5
DENSITY_FLOOR = 0.5
INT_PZ_OUTLINE = 0.3
OVERLAY_MIN_SIZE = 10
XMM_PATH = BASE / 'xray/images/XMM/smoothed/XMM_comb-net-center_smoothed_10arcsec.fits'
M2_PATH = BASE / (
    'M2/images/minkasi/2025-11_JackOS/cutouts/smoothed/'
    'MOO_1142_signal_4p5_by_eye_ctr_smoothed_10arcsec.fits'
)
XMM_VMIN = 1.7e-6
XMM_VMAX = 1.2e-4
SZ_LEVELS = np.round(np.arange(-5.5e-4, -5.1e-5, 5e-5), 12)
SZ_LEVELS_NEG = SZ_LEVELS
SZ_LEVELS_POS = np.round(np.arange(5e-5, 2.51e-4, 5e-5), 12)
XRAY_CONTOUR_MIN = 1e-5
XRAY_CONTOUR_COUNT = 15
XRAY_MAIN = SkyCoord('11:42:47.4384 +15:27:12.96', unit=(u.hourangle, u.deg))
XRAY_WEST = SkyCoord('11:42:45.216 +15:27:20.16', unit=(u.hourangle, u.deg))
MC1_cluster_center = SkyCoord(175.68728, 15.45493, unit='deg')
by_eye_center = SkyCoord('11:42:46.45 15:27:12.38', unit=(u.hourangle, u.deg))
ds9_a = LinearSegmentedColormap.from_list(
    'ds9_a',
    [
        (0.00, (0.0, 0.0, 0.0)),
        (0.25, (0.0, 1.0, 0.0)),
        (0.50, (1.0, 0.0, 1.0)),
        (0.64, (1.0, 0.0, 0.5)),
        (0.77, (1.0, 0.0, 0.0)),
        (1.00, (1.0, 1.0, 0.0)),
    ],
    N=256,
)
# -------------------------------------------------------------------------------------------- #
def build_density_levels(density_data):
    density_max = np.nanmax(density_data)
    if not np.isfinite(density_max):
        raise ValueError('Density map contains no finite values.')

    density_upper = max(
        DENSITY_FLOOR + DENSITY_STEP,
        np.ceil(density_max / DENSITY_STEP) * DENSITY_STEP,
    )
    level_count = int(round((density_upper - DENSITY_FLOOR) / DENSITY_STEP)) + 1
    return DENSITY_FLOOR + DENSITY_STEP * np.arange(level_count)


def get_sources(members):
    """Valid photz sources with marker sizes scaled by INT_PZ."""
    galaxy_ra = np.asarray(np.ma.filled(members['GAL_RA'], np.nan), dtype=float)
    galaxy_dec = np.asarray(np.ma.filled(members['GAL_DEC'], np.nan), dtype=float)
    int_pz = np.asarray(np.ma.filled(members['INT_PZ'], np.nan), dtype=float)
    valid = np.isfinite(galaxy_ra) & np.isfinite(galaxy_dec) & np.isfinite(int_pz)
    valid &= int_pz >= 0
    galaxy_ra, galaxy_dec, int_pz = galaxy_ra[valid], galaxy_dec[valid], int_pz[valid]
    max_int_pz = int_pz.max() if int_pz.size else 0.0
    marker_sizes = 100 * int_pz / max_int_pz if max_int_pz > 0 else np.zeros(len(int_pz))
    return galaxy_ra, galaxy_dec, int_pz, marker_sizes


def plot_galdens(hdu, members, bcg, figure, subplot):
    density_data, density_header = hdu[0].data, hdu[0].header
    density_levels = build_density_levels(density_data)
    density_hdu = fits.PrimaryHDU(data=density_data, header=density_header)
    density_img = aplpy.FITSFigure(data=density_hdu, figure=figure, subplot=subplot)

    bcg_ra = float(bcg['GAL_RA'][0])
    bcg_dec = float(bcg['GAL_DEC'][0])
    density_img.recenter(
        x=by_eye_center.ra.deg,
        y=by_eye_center.dec.deg,
        radius=TOP_RADIUS.to_value(u.degree),
    )
    density_img.show_contour(
        data=density_hdu,
        cmap=cmap,
        filled=True,
        levels=density_levels,
    )
    density_contour = density_img._layers['contour_set_1']
    density_contour.cmap.set_under('white')
    density_img.ax.tick_params(axis='both', which='both', direction='in', color='k')

    galaxy_ra, galaxy_dec, int_pz, marker_sizes = get_sources(members)
    if len(galaxy_ra):
        density_img.show_markers(
            galaxy_ra,
            galaxy_dec,
            coords_frame='world',
            s=marker_sizes,
            alpha=0.3,
            edgecolor='none',
            facecolor='black',
        )
        # Same sizes, so the outline stays on the circle perimeter
        high = int_pz >= INT_PZ_OUTLINE
        if high.any():
            density_img.show_markers(
                galaxy_ra[high],
                galaxy_dec[high],
                coords_frame='world',
                s=marker_sizes[high],
                edgecolor='black',
                facecolor='none',
                linewidths=0.8,
            )

    density_img.show_markers(
        bcg_ra,
        bcg_dec,
        coords_frame='world',
        s=120,
        marker='*',
        edgecolor='yellow',
        facecolor='none',
        linewidths=1.5,
    )
    density_img.axis_labels.hide()
    return density_img, density_contour


def reproject_xray(wcs, shape):
    with fits.open(XMM_PATH) as hdu:
        xray, _ = reproject_interp((hdu[0].data, hdu[0].header), wcs, shape_out=shape)
    return xray


def load_sz():
    with fits.open(M2_PATH) as hdu:
        sz_data = np.nan_to_num(np.squeeze(hdu[0].data))
        sz_wcs = WCS(hdu[0].header, naxis=2)
    return gaussian_filter(sz_data, 1), sz_wcs


def plot_overlay(ax, xray, density_data, density_levels, sz_data, sz_wcs, bcg, members=None):
    """Galaxy density, X-ray (black) and SZ (red) contours; sources only if members is given."""
    ax.contourf(density_data, levels=density_levels, cmap=cmap, zorder=1)
    if members is not None:
        galaxy_ra, galaxy_dec, int_pz, _ = get_sources(members)
        high = int_pz >= INT_PZ_OUTLINE
        # Sizes span the INT_PZ >= cut sources only: cut -> OVERLAY_MIN_SIZE, max -> 100
        span = max(int_pz.max() - INT_PZ_OUTLINE, 1e-12)
        sizes = OVERLAY_MIN_SIZE + (100 - OVERLAY_MIN_SIZE) * (int_pz[high] - INT_PZ_OUTLINE) / span
        ax.scatter(
            galaxy_ra[high],
            galaxy_dec[high],
            transform=ax.get_transform('world'),
            s=sizes,
            facecolor='black',
            edgecolor='none',
            alpha=0.4,
            zorder=1.5,
        )
    # X-ray levels run linearly up to the maximum within the displayed field
    (x0, x1), (y0, y1) = ax.get_xlim(), ax.get_ylim()
    field = xray[int(max(y0, 0)):int(y1) + 1, int(max(x0, 0)):int(x1) + 1]
    xray_levels = np.linspace(XRAY_CONTOUR_MIN, np.nanmax(field), XRAY_CONTOUR_COUNT)
    ax.contour(xray, levels=xray_levels, colors='black', linestyles='solid', linewidths=1.0, zorder=2)
    for levels, linestyle in ((SZ_LEVELS_NEG, 'dashed'), (SZ_LEVELS_POS, 'solid')):
        ax.contour(
            sz_data,
            levels=levels,
            colors='red',
            linestyles=linestyle,
            linewidths=1.2,
            transform=ax.get_transform(sz_wcs),
            zorder=3,
        )
    ax.scatter(
        float(bcg['GAL_RA'][0]),
        float(bcg['GAL_DEC'][0]),
        transform=ax.get_transform('world'),
        s=120,
        marker='*',
        edgecolor='yellow',
        facecolor='none',
        linewidths=1.5,
        zorder=4,
    )


def mark_xray_centers(ax):
    ax.scatter(
        MC1_cluster_center.ra.deg,
        MC1_cluster_center.dec.deg,
        transform=ax.get_transform('world'),
        s=100,
        marker='+',
        color='darkorange',
        linewidths=2.5,
        zorder=5,
    )
    for center, color in ((XRAY_MAIN, 'magenta'), (XRAY_WEST, 'mediumslateblue')):
        ax.scatter(
            center.ra.deg,
            center.dec.deg,
            transform=ax.get_transform('world'),
            s=100,
            marker='x',
            color=color,
            linewidths=2.5,
            zorder=5,
        )


def plot_xray_sz(ax, xray, sz_data, sz_wcs):
    """Draw the X-ray image and SZ contours on an axes that shares panel 1's WCS."""
    ax.imshow(
        xray,
        origin='lower',
        cmap=ds9_a,
        norm=PowerNorm(gamma=0.5, vmin=XMM_VMIN, vmax=XMM_VMAX, clip=True),
    )

    ax.contour(
        sz_data,
        levels=SZ_LEVELS,
        colors='white',
        linestyles='dashed',
        linewidths=1.2,
        transform=ax.get_transform(sz_wcs),
    )


def add_scale_bar(ax, pixel_scale, color):
    """Draw a SCALE_BAR-long bar in the upper right corner of the current view."""
    (x0, x1), (y0, y1) = ax.get_xlim(), ax.get_ylim()
    length = SCALE_BAR.to_value(u.deg) / pixel_scale
    x_end = x1 - 0.05 * (x1 - x0)
    y_bar = y1 - 0.08 * (y1 - y0)
    ax.plot([x_end - length, x_end], [y_bar, y_bar], color=color, lw=2, solid_capstyle='butt', zorder=10)
    ax.text(
        x_end - length / 2,
        y_bar - 0.02 * (y1 - y0),
        f'{SCALE_BAR.to_value(u.arcsec):g}"',
        color=color,
        ha='center',
        va='top',
        zorder=10,
    )


def make_figure():
    members = Table.read(PHOTZ_PATH)
    bcg = Table.read(BCG_PATH)

    figure = plt.figure(figsize=(10, 9))
    figure.subplots_adjust(left=0.09, right=0.93, bottom=0.07, top=0.96, wspace=0.05, hspace=0.05)

    with fits.open(DENSITY_PATH) as hdu:
        density_img, _ = plot_galdens(
            hdu=hdu, members=members, bcg=bcg, figure=figure, subplot=(2, 2, 1)
        )
        wcs = density_img._wcs
        density_data = hdu[0].data
    ref_ax = density_img.ax

    # Each row shares axes; the bottom row is zoomed in relative to the top row
    xray_ax = figure.add_subplot(2, 2, 2, projection=wcs, sharex=ref_ax, sharey=ref_ax)
    overlay_ax = figure.add_subplot(2, 2, 3, projection=wcs)
    empty_ax = figure.add_subplot(2, 2, 4, projection=wcs, sharex=overlay_ax, sharey=overlay_ax)
    axes = [ref_ax, xray_ax, overlay_ax, empty_ax]
    xray = reproject_xray(wcs, density_data.shape)
    sz_data, sz_wcs = load_sz()
    plot_xray_sz(xray_ax, xray, sz_data, sz_wcs)
    for ax in axes:
        ax.set_box_aspect(1)
        ax.tick_params(axis='both', which='both', direction='in', color='k')

    zoom = (BOTTOM_RADIUS / TOP_RADIUS).decompose().value
    (x0, x1), (y0, y1) = ref_ax.get_xlim(), ref_ax.get_ylim()
    (x_mid, y_mid), (x_half, y_half) = ((x0 + x1) / 2, (y0 + y1) / 2), ((x1 - x0) / 2, (y1 - y0) / 2)
    for ax in (xray_ax, ref_ax):
        ax.set_xlim(x0, x1)
        ax.set_ylim(y0, y1)
    overlay_ax.set_xlim(x_mid - zoom * x_half, x_mid + zoom * x_half)
    overlay_ax.set_ylim(y_mid - zoom * y_half, y_mid + zoom * y_half)
    plot_overlay(
        overlay_ax, xray, density_data, build_density_levels(density_data), sz_data, sz_wcs, bcg, members
    )
    plot_overlay(
        empty_ax, xray, density_data, build_density_levels(density_data), sz_data, sz_wcs, bcg
    )
    mark_xray_centers(empty_ax)
    for ax in (overlay_ax, empty_ax):
        ax.set_xlim(x_mid - zoom * x_half, x_mid + zoom * x_half)
        ax.set_ylim(y_mid - zoom * y_half, y_mid + zoom * y_half)
    xray_ax.tick_params(axis='both', which='both', color='white')

    pixel_scale = proj_plane_pixel_scales(wcs)[0]
    for ax, color in zip(axes, ('black', 'white', 'black', 'black')):
        add_scale_bar(ax, pixel_scale, color)
    other_axes = [empty_ax]
    # Hide inner tick labels: x only on the bottom row, y only on the left column
    for ax in (ref_ax, xray_ax):
        ax.coords[0].set_ticklabel_visible(False)
    for ax in (xray_ax, other_axes[0]):
        ax.coords[1].set_ticklabel_visible(False)
    for ax in [overlay_ax, *other_axes]:
        ax.coords[0].set_axislabel('')
        ax.coords[1].set_axislabel('')

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(OUTPUT_PATH, format='pdf', dpi=300, bbox_inches='tight')
    plt.close(figure)
    print(f'Saved figure to: {OUTPUT_PATH}')


if __name__ == '__main__':
    make_figure()
