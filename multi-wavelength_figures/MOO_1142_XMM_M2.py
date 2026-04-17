"""
Plot XMM image and M2 image which were both previously smoothed to the M2 beam for MOO 1142 in a two-panel APLpy figure.

Created in 2026-04
"""
from __future__ import annotations

import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from pathlib import Path

import aplpy

from astropy import units as u
from astropy.coordinates import SkyCoord
# -------------------------------------------------------------------------------------------- #
location = Path(__file__).resolve().parent.parent
xmm_image = location / "xray/images/XMM/smoothed/XMM_comb-net-center_smoothed_10arcsec.fits"
#m2_image = location / "M2/images/minkasi/2025-11_JackOS/cutouts/smoothed/MOO_1142_signal_4p5_by_eye_ctr_smoothed_10arcsec.fits"
m2_image = location / "M2/images/minkasi/2025-11_JackOS/cutouts/smoothed/MOO_1142_PS_sub_4p5_by_eye_ctr_smoothed_10arcsec.fits"

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
XMM_VMIN = 0
XMM_VMAX = 0.00011
# M2 image display parameters
M2_VMIN = -0.0006 # native image with no smoothing: -0.0008
M2_VMAX = 6E-5 # native image with no smoothing: 0.0006

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
	panel.show_markers(bcg.ra, bcg.dec, coords_frame='world', s=100, marker='*', edgecolor='yellow', linewidth=1)
	panel.show_markers(main_cluster.ra, main_cluster.dec, coords_frame='world', s=100, marker='x', facecolor='k', linewidth=2)
	panel.show_markers(west_cluster.ra, west_cluster.dec, coords_frame='world', s=100, marker='+', facecolor='k', linewidth=2)

def format_colorbar_scientific(panel: aplpy.FITSFigure, power: int) -> None:
	"""Force an APLpy colorbar to use a specific scientific-notation power."""
	formatter = ScalarFormatter(useMathText=False)
	formatter.set_scientific(True)
	formatter.set_powerlimits((power, power))
	panel.colorbar._colorbar.formatter = formatter
	panel.colorbar._colorbar.update_ticks()

def main():
	# Figure
	fig = plt.figure(figsize=(12, 5.4))

	# XMM
	left_panel = aplpy.FITSFigure(str(xmm_image), figure=fig, subplot=[0.07, 0.12, 0.41, 0.8])
	left_panel.show_colorscale(
		stretch="sqrt",
		#cmap = inferno
		cmap=plt.cm.nipy_spectral,
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
	left_panel.set_title("XMM (0.4-7.2 keV)")
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
		#cmap=plt.cm.nipy_spectral,
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
	right_panel.add_colorbar()
	right_panel.colorbar.set_axis_label_text("Compton y")
	format_colorbar_scientific(right_panel, -4)
	right_panel.set_title("MUSTANG-2 (90 GHz)")
	right_panel.axis_labels.set_xtext("Right Ascension")
	style_panel(right_panel, show_y_axis=False)

	# Show points of interest
	show_points_of_interest(left_panel)
	show_points_of_interest(right_panel)

	# Save/show plot
	#fig.savefig(location / "multi-wavelength_figures/XMM_M2/MOO_1142_XMM_comb-net_M2_PS_sub_RB_BCG_subclusters.pdf", dpi=300, bbox_inches="tight")
	plt.show()

if __name__ == "__main__":
	main()