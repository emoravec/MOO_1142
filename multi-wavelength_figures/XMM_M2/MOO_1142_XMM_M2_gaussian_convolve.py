"""
Gaussian-smooth the XMM and MUSTANG-2 images for MOO 1142.

The smoothing kernel sigma is derived from a requested on-sky FWHM using the
image pixel scale from the FITS header CDELT1 keyword:

sigma_pix = fwhm_pixels / (2.0 * np.sqrt(2.0 * np.log(2.0)))

Created in 2026-04.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from astropy import units as u
from astropy.convolution import Gaussian2DKernel, convolve
from astropy.io import fits
# -------------------------------------------------------------------------------------------- #
location = Path(__file__).resolve().parent.parent
xmm_image = location / "xray/images/XMM/XMM_comb-net-center.fits"
xmm_output_dir = location / "xray/images/XMM/smoothed"
m2_image = location / "M2/images/minkasi/2025-11_JackOS/cutouts/MOO_1142_signal_4p5_by_eye_ctr.fits"
m2_output_dir = location / "M2/images/minkasi/2025-11_JackOS/cutouts/smoothed"

SMOOTH_FWHM = 10.0 * u.arcsec
# -------------------------------------------------------------------------------------------- #
def get_pixel_scale_arcsec(header: fits.Header) -> u.Quantity:
	"""Return the absolute pixel scale from CDELT1 in arcsec/pixel."""
	if "CDELT1" not in header:
		raise KeyError("CDELT1 keyword not found in FITS header")

	pixel_scale_arcsec = abs(header["CDELT1"]) * u.deg.to(u.arcsec) * u.arcsec
	return pixel_scale_arcsec

def get_sigma_pixels(image_path: Path, smooth_fwhm: u.Quantity) -> float:
	"""Convert an on-sky FWHM to a Gaussian sigma in image pixels."""
	header = fits.getheader(image_path)
	pixel_scale_arcsec = get_pixel_scale_arcsec(header)
	fwhm_pixels = smooth_fwhm.to_value(u.arcsec) / pixel_scale_arcsec.to_value(u.arcsec)
	sigma_pix = fwhm_pixels / (2.0 * np.sqrt(2.0 * np.log(2.0)))
	rounded_sigma_pix = round(float(sigma_pix), 2)
	print(40 * "-")

	print(
		f"{image_path.name}: pixel scale = {pixel_scale_arcsec.to_value(u.arcsec):.3f} arcsec/pixel, "
		f"FWHM = {smooth_fwhm.to_value(u.arcsec):.2f} arcsec = {fwhm_pixels:.3f} pixels, "
		f"sigma = {rounded_sigma_pix:.2f} pixels"
	)
	return rounded_sigma_pix

def build_gaussian_kernel(sigma_pix: float) -> Gaussian2DKernel:
	"""Build an explicitly centered odd-sized Gaussian kernel."""
	# Convert sigma to an integer kernel radius in pixels
	# To do this, use a 4-sigma radius so the truncated Gaussian tail is negligible
	# Make the full kernel size odd so the kernel remains centered on a pixel.
	# Gaussian2DKernel defaults to a size of floor(8 * stddev + 1), which is effectively about a 4 * sigma radius on each side.
	kernel_radius = max(1, int(np.ceil(4.0 * sigma_pix)))
	# turn kernel_radius into the full odd-sized kernel width needed by Gaussian2DKernel.
	kernel_size = 2 * kernel_radius + 1
	return Gaussian2DKernel(
		sigma_pix,
		x_size=kernel_size,
		y_size=kernel_size,
		mode="center",
	)

def smooth_fits_image(image_path: Path, smooth_fwhm: u.Quantity, destination: Path) -> Path:
	"""Convolve a FITS image with a Gaussian kernel and write the smoothed FITS file."""
	with fits.open(image_path) as hdul:
		header = hdul[0].header.copy()
		data = np.asarray(hdul[0].data, dtype=float)

	if data.ndim != 2:
		raise ValueError(f"Expected a 2D image in {image_path}, found shape {data.shape}")

	sigma_pix = get_sigma_pixels(image_path, smooth_fwhm)
	kernel = build_gaussian_kernel(sigma_pix)
	smoothed = convolve(
		data,
		kernel,
		boundary="extend",
		normalize_kernel=True,
		nan_treatment="interpolate",
		preserve_nan=True,
	)

	header["HISTORY"] = (
		f"Gaussian smoothed with FWHM={smooth_fwhm.to_value(u.arcsec):.2f} arcsec "
		f"(sigma={sigma_pix:.2f} pix)"
	)
	fits.PrimaryHDU(data=smoothed.astype(np.float32), header=header).writeto(destination, overwrite=True)
	print(f"Wrote {destination}")
	return destination

def build_output_path(image_path: Path, smooth_fwhm: u.Quantity, output_dir: Path) -> Path:
	"""Create a descriptive output filename in the requested output directory."""
	fwhm_label = int(smooth_fwhm.to_value(u.arcsec))
	return output_dir / f"{image_path.stem}_smoothed_{fwhm_label}arcsec.fits"

def main() -> None:
	# image_configs = (
	# 	(xmm_image, xmm_output_dir),
	# 	(m2_image, m2_output_dir),
	# )

	image_configs = (
		(location / "M2/images/minkasi/2025-11_JackOS/cutouts/MOO_1142_PS_sub_4p5_by_eye_ctr.fits", location / "M2/images/minkasi/2025-11_JackOS/cutouts/smoothed"),
		(location / "M2/images/minkasi/2025-11_JackOS/MOO_1142_PS_sub.fits", location / "M2/images/minkasi/2025-11_JackOS/smoothed"),
		(location / "M2/images/minkasi/2025-11_JackOS/MOO_1142_signal.fits", location / "M2/images/minkasi/2025-11_JackOS/smoothed")
	)

	for image_path, output_dir in image_configs:
		output_dir.mkdir(parents=True, exist_ok=True)
		destination = build_output_path(image_path, SMOOTH_FWHM, output_dir)
		smooth_fits_image(image_path, SMOOTH_FWHM, destination)

if __name__ == "__main__":
	main()