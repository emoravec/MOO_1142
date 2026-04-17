"""
Plot XMM adaptively smoothed image and M2 data (using aplpy smooth) for MOO 1142 in a two-panel APLpy figure.

Created in 2026-04
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

import aplpy

from astropy import units as u
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
# -------------------------------------------------------------------------------------------- #
location = Path(__file__).resolve().parents[2]
xmm_image = location / "xray/images/XMM/XMM_adapt-400-7200.fits"
#xmm_image = location / "xray/images/XMM/XMM_comb-net-center.fits"
m2_image = location / "M2/images/minkasi/2025-11_JackOS/MOO_1142_signal.fits"

# Image centers
# from /Users/emoravec/Documents/Research/MOO_1142/profile_fitting/XMM/XMM_fit_coords_EBarbavara.txt
# which is Eleonora's fit to the main and west subcluster
main = SkyCoord(175.69766238,15.45364413,unit='deg')
west = SkyCoord(175.68837814,15.45566318,unit='deg')
by_eye_center = SkyCoord('11:42:46.45 15:27:12.38', unit=(u.hourangle,u.deg))
# Size of image
RA_WIDTH = 4.0 * u.arcmin
DEC_HEIGHT = 3.0 * u.arcmin
# XMM adaptively smoothed image display parameters
XMM_VMIN = 1.0
XMM_VMAX = 320.0
XMM_VMID = 0.0
# M2 image display parameters
M2_FWHM = 9.0 * u.arcsec
M2_VMIN = -0.0008 # native image with no smoothing: -0.0008
M2_VMAX = 0.000023 # native image with no smoothing: 0.0006
M2_CALCULATE_SMOOTHING = False
# Manual APLpy smoothing value
MANAUAL_APLPY_SMOOTHING = 3
# -------------------------------------------------------------------------------------------- #
def get_gaussian_sigma_pixels(image_path: Path, smooth_fwhm: u.Quantity) -> float:
	"""Convert a Gaussian FWHM on-sky smoothing scale to sigma in image pixels."""
	header = fits.getheader(image_path)
	wcs = WCS(header).celestial
	## Calculate the pixel scale in arcsec/pixel using the WCS information from the image header.
	# proj_plane_pixel_scales(wcs) returns the pixel scale along each WCS axis, usually in degrees per pixel. There are typically two axes, so we take the mean and convert to arcsec per pixel. 
	pixel_scale = np.mean(proj_plane_pixel_scales(wcs)) * u.deg
	print(f"Image WCS pixel scale: {pixel_scale.to_value(u.arcsec):.3f}\"/pixel")
	fwhm_pixels = (smooth_fwhm / pixel_scale).decompose().value
	print(f"Image pixel scale: {pixel_scale.to_value(u.arcsec):.3f}\"/pixel, "f"requested FWHM: {smooth_fwhm.to_value(u.arcsec):.1f}\" = {fwhm_pixels:.3f} pixels")

	## Calculate the required sigma in pixels for the Gaussian kernel.
	# APLpy's Gaussian kernel takes an input of sigma in pixels, but if you are working with a FWHM, you need to convert requested FWHM to a sigma in pixels.
	# For example, for the 2025-11 M2 map the pixel scale is 2 arcsec per pixel, which means that a 9 arcsec Gaussian FWHM corresponds to 4.5 pixels FWHM. 
	# And it is know that the relationship between FWHM and sigma for a Gaussian is FWHM = 2 * sqrt(2 * ln(2)) * sigma or sigma = FWHM / (2 * sqrt(2 * ln(2))). ~= FWHM / 2.3548.
	# Which means that in this case the required sigma in pixels is approximately = FWHM / (2 * sqrt(2 * ln(2))) ~= 4.5 / 2.3548 ~= 1.91 pixels.
	sigma_pix = fwhm_pixels / (2.0 * np.sqrt(2.0 * np.log(2.0)))
	print(f"Calculated Gaussian sigma for APLpy kernel: {sigma_pix:.3f} pixels")
	return float(sigma_pix)

def get_aplpy_smooth_pixels(image_path: Path, smooth_fwhm: u.Quantity) -> int:
	"""Return the nearest integer sigma in pixels for APLpy's smooth= argument."""
	sigma_pixels = get_gaussian_sigma_pixels(image_path, smooth_fwhm)
	# APLpy 2.1/2.2 uses the smooth value to build an integer-sized kernel internally,
	# and Astropy requires that kernel size to be odd. Since APLpy sets x_size=y_size=5*smooth,
	# the usable smooth values are therefore odd integers, so we approximate by the nearest odd sigma.
	nearest_integer = int(np.round(sigma_pixels))
	if nearest_integer % 2 == 0:
		lower_odd = max(1, nearest_integer - 1)
		upper_odd = nearest_integer + 1
		if abs(sigma_pixels - lower_odd) <= abs(sigma_pixels - upper_odd):
			nearest_integer = lower_odd
		else:
			nearest_integer = upper_odd
	return nearest_integer

def get_xmm_log_limits(data: np.ndarray) -> tuple[float, float, float]:
	"""Return log-safe display limits for the XMM image."""
	positive = data[np.isfinite(data) & (data > 0)]
	vmin = float(np.min(positive))
	vmax = float(np.nanmax(data))
	return vmin, vmax, 0.0

def get_linear_limits(data: np.ndarray) -> tuple[float, float]:
	"""Return simple linear display limits from the finite data range."""
	finite = data[np.isfinite(data)]
	return float(np.nanmin(finite)), float(np.nanmax(finite))

def style_panel(panel: aplpy.FITSFigure, show_y_axis: bool = True) -> None:
	"""Apply common styling to an APLpy panel."""
	panel.ticks.set_color("black")
	panel.tick_labels.set_font(size=11)
	panel.axis_labels.set_font(size=12)
	panel.frame.set_color("black")
	if not show_y_axis:
		panel.axis_labels.hide_y()
		panel.tick_labels.hide_y()

def main():
	xmm_data = np.squeeze(fits.getdata(xmm_image)).astype(float)
	m2_data = np.squeeze(fits.getdata(m2_image)).astype(float)
	if M2_CALCULATE_SMOOTHING:
		m2_sigma_pixels = get_gaussian_sigma_pixels(m2_image, M2_FWHM)
		m2_aplpy_smooth = get_aplpy_smooth_pixels(m2_image, M2_FWHM)
		print(
			f"Applying M2 Gaussian smoothing: {M2_FWHM.to_value(u.arcsec):.1f}\" FWHM "
			f"= {m2_sigma_pixels:.3f} pixel sigma; using APLpy smooth={m2_aplpy_smooth} "
			f"as the nearest-pixel approximation"
		)
	else:
		m2_aplpy_smooth = MANAUAL_APLPY_SMOOTHING
		print(f"M2 Gaussian smoothing not applied; using default smooth={MANAUAL_APLPY_SMOOTHING}")

	# Figure
	fig = plt.figure(figsize=(12, 5.4))

	# XMM
	left_panel = aplpy.FITSFigure(str(xmm_image), figure=fig, subplot=[0.07, 0.12, 0.41, 0.8])
	left_panel.show_colorscale(
		stretch="log",
		#cmap = inferno
		cmap=plt.cm.nipy_spectral,
		vmin=XMM_VMIN,
		vmax=XMM_VMAX,
		vmid=XMM_VMID,
	)
	left_panel.recenter(
		by_eye_center.ra,
		by_eye_center.dec,
		width=RA_WIDTH.to_value(u.deg),
		height=DEC_HEIGHT.to_value(u.deg),
	)
	left_panel.add_colorbar()
	left_panel.colorbar.set_axis_label_text("Intensity (cts/pix/s)")
	left_panel.set_title("XMM 0.4-7.2 keV")
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
		#stretch="sqrt",
		#cmap="RdBu_r",
		cmap=plt.cm.nipy_spectral,
		vmin=M2_VMIN,
		vmax=M2_VMAX,
		smooth=m2_aplpy_smooth,
		kernel="gauss",
	)
	right_panel.recenter(
		by_eye_center.ra,
		by_eye_center.dec,
		width=RA_WIDTH.to_value(u.deg),
		height=DEC_HEIGHT.to_value(u.deg),
	)
	right_panel.add_colorbar()
	right_panel.colorbar.set_axis_label_text("Compton y")
	right_panel.set_title("MUSTANG-2, 90 GHz")
	right_panel.axis_labels.set_xtext("Right Ascension")
	style_panel(right_panel, show_y_axis=False)

	# Save/show plot
	fig.savefig(location / "multi-wavelength_figures/XMM_M2/MOO_1142_XMM_comb-net_M2_smooth=3.pdf", dpi=300, bbox_inches="tight")
	#plt.show()

if __name__ == "__main__":
	main()