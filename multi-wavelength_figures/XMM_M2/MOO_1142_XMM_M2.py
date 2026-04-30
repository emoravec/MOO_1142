"""
Plot XMM image and M2 image which were both previously smoothed to the M2 beam for MOO 1142 in a two-panel APLpy figure.

Created in 2026-04
"""
from __future__ import annotations

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import ScalarFormatter
from pathlib import Path

import aplpy

from astropy import units as u
from astropy.coordinates import SkyCoord
# -------------------------------------------------------------------------------------------- #
ds9_a = LinearSegmentedColormap.from_list(
    "ds9_a",
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

def make_contour_levels(start: float, stop: float, step: float) -> list[float]:
	"""Generate inclusive contour levels in increasing order."""
	if step <= 0:
		raise ValueError("step must be positive")
	if stop < start:
		raise ValueError("stop must be greater than or equal to start")

	levels: list[float] = []
	current = start
	epsilon = step / 10
	while current <= stop + epsilon:
		levels.append(round(current, 12))
		current += step
	return levels
# -------------------------------------------------------------------------------------------- #
location = Path(__file__).resolve().parents[2]
xmm_image = location / "xray/images/XMM/smoothed/XMM_comb-net-center_smoothed_10arcsec.fits"
m2_image = location / "M2/images/minkasi/2025-11_JackOS/cutouts/smoothed/MOO_1142_signal_4p5_by_eye_ctr_smoothed_10arcsec.fits"
#m2_image = location / "M2/images/minkasi/2025-11_JackOS/cutouts/smoothed/MOO_1142_PS_sub_4p5_by_eye_ctr_smoothed_10arcsec.fits"
output_dir = location / "multi-wavelength_figures/XMM_M2/pdfs"

# Image centers
# from /Users/emoravec/Documents/Research/MOO_1142/profile_fitting/XMM/XMM_fit_coords_EBarbavara.txt
# which is Eleonora's fit to the main and west subcluster
main_cluster = SkyCoord(175.69766238,15.45364413,unit='deg')
west_cluster = SkyCoord(175.68837814,15.45566318,unit='deg')
by_eye_center = SkyCoord('11:42:46.45 15:27:12.38', unit=(u.hourangle,u.deg))
# Size of image
RA_WIDTH = 4.0 * u.arcmin
DEC_HEIGHT = 3.0 * u.arcmin
# XMM adaptively smoothed image display parameters
XMM_VMIN = 2E-6
XMM_VMAX = 8E-5
# M2 image display parameters
# M2_VMIN = -0.0006 # native image with no smoothing: -0.0008
# M2_VMAX = 6E-5 # native image with no smoothing: 0.0006
M2_VMIN = -3e-4 #-6e-4
M2_VMAX = 3e-4 #1e-4
M2_CONTOUR_LEVELS_NEG = make_contour_levels(start=-5.5e-4, stop=-5e-5, step=5e-5)
M2_CONTOUR_LEVELS_POS = make_contour_levels(start=5e-5, stop=2.5e-4, step=5e-5)

# BCG
bcg = SkyCoord(175.69784, 15.45318, unit='deg')

# -------------------------------------------------------------------------------------------- #
def style_panel(panel: aplpy.FITSFigure, show_y_axis: bool = True) -> None:
	"""Apply common styling to an APLpy panel."""
	panel.ticks.set_color("black")
	panel.tick_labels.set_font(size=11)
	panel.axis_labels.set_font(size=12)
	panel.frame.set_color("black")
	if not show_y_axis:
		panel.axis_labels.hide_y()
		panel.tick_labels.hide_y()

def show_points_of_interest(panel: aplpy.FITSFigure) -> None:
	"""Overlay the BCG and subcluster centers on a panel."""
	panel.show_markers(bcg.ra, bcg.dec, coords_frame='world', s=100, marker='*', edgecolor='k', linewidth=1)
	#panel.show_markers(main_cluster.ra, main_cluster.dec, coords_frame='world', s=100, marker='x', facecolor='k', linewidth=2)
	#panel.show_markers(west_cluster.ra, west_cluster.dec, coords_frame='world', s=100, marker='+', facecolor='k', linewidth=2)

def format_colorbar_scientific(panel: aplpy.FITSFigure, power: int) -> None:
	"""Force an APLpy colorbar to use a specific scientific-notation power."""
	formatter = ScalarFormatter(useMathText=False)
	formatter.set_scientific(True)
	formatter.set_powerlimits((power, power))
	panel.colorbar._colorbar.formatter = formatter
	panel.colorbar._colorbar.update_ticks()

def make_zero_white_cmap(vmin: float, vmax: float, base_cmap: str = "RdBu_r") -> LinearSegmentedColormap:
	"""Shift a diverging colormap so zero stays white for asymmetric limits."""
	if not (vmin < 0 < vmax):
		raise ValueError("Expected asymmetric limits with zero inside the display range.")

	zero_position = -vmin / (vmax - vmin)
	base = plt.get_cmap(base_cmap)
	negative_steps = max(2, int(round(256 * zero_position)))
	positive_steps = max(2, 256 - negative_steps)

	negative_segment = [
		(zero_position * index / (negative_steps - 1), base(0.5 * index / (negative_steps - 1)))
		for index in range(negative_steps)
	]
	positive_segment = [
		(
			zero_position + (1.0 - zero_position) * index / (positive_steps - 1),
			base(0.5 + 0.5 * index / (positive_steps - 1)),
		)
		for index in range(1, positive_steps)
	]

	return LinearSegmentedColormap.from_list(
		f"{base_cmap}_zero_white",
		negative_segment + positive_segment,
		N=256,
	)

def main():
	# Figure
	fig = plt.figure(figsize=(16, 7))

	# XMM
	left_panel = aplpy.FITSFigure(str(xmm_image), figure=fig, subplot=[0.07, 0.12, 0.41, 0.8])
	left_panel.show_colorscale(
		stretch="sqrt",
		#cmap = inferno
		#cmap=plt.cm.nipy_spectral,
		cmap=ds9_a,
		vmin=XMM_VMIN,
		vmax=XMM_VMAX,
	)
	left_panel.recenter(
		by_eye_center.ra,
		by_eye_center.dec,
		width=RA_WIDTH.to_value(u.deg),
		height=DEC_HEIGHT.to_value(u.deg),
	)
	left_panel.add_colorbar()
	left_panel.colorbar.set_axis_label_text("Intensity (cts/pix/s)")
	format_colorbar_scientific(left_panel, -5)
	#left_panel.set_title("XMM (0.4-7.2 keV)")
	left_panel.axis_labels.set_xtext("Right Ascension")
	left_panel.axis_labels.set_ytext("Declination")
	style_panel(left_panel)

	# MUSTANG-2
	right_panel = aplpy.FITSFigure(
		str(m2_image),
		figure=fig,
		subplot=[0.54, 0.12, 0.41, 0.8],
		convention="calabretta",
	)
	right_panel.show_colorscale(
		stretch="linear",
		cmap="RdBu_r",
		#cmap=make_zero_white_cmap(M2_VMIN, M2_VMAX),
		vmin=M2_VMIN,
		vmax=M2_VMAX,
		smooth=None,
	)
	right_panel.recenter(
		by_eye_center.ra,
		by_eye_center.dec,
		width=RA_WIDTH.to_value(u.deg),
		height=DEC_HEIGHT.to_value(u.deg),
	)
	right_panel.show_contour(
		str(m2_image),
		levels=M2_CONTOUR_LEVELS_NEG,
		colors="yellow",
		linestyles="dashed",
		linewidths=1.2,
		smooth=1,
		zorder=1,
		convention="calabretta",
	)
	# right_panel.show_contour(
	# 	str(m2_image),
	# 	levels=M2_CONTOUR_LEVELS_POS,
	# 	colors="yellow",
	# 	linestyles="solid",
	# 	linewidths=1.2,
	# 	smooth=1,
	# 	zorder=1,
	# )
	right_panel.add_colorbar()
	right_panel.colorbar.set_axis_label_text("Compton y")
	format_colorbar_scientific(right_panel, -4)
	#right_panel.set_title("MUSTANG-2 (90 GHz)")
	right_panel.axis_labels.set_xtext("Right Ascension")
	style_panel(right_panel, show_y_axis=False)

	# Show points of interest
	show_points_of_interest(left_panel)
	show_points_of_interest(right_panel)

	# Save/show plot
	fig.savefig(output_dir / "MOO_1142_XMM_M2_PS_sub_contours_BCG.pdf", dpi=300, bbox_inches="tight")
	#plt.show()

if __name__ == "__main__":
	main()