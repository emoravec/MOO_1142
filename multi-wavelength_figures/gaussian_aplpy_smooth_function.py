# Function from MOO_1142_multi-wavelength_plots.py that calculates the exact Gaussian smoothing of the M2 map and rounds it to the nearest odd number for aplpy smoothing
def smooth_m2_map_direct(image_path: Path, smooth_fwhm: u.Quantity):
	"""Exact alternative: smooth the M2 map directly instead of using APLpy smooth=."""
	from astropy.convolution import Gaussian2DKernel, convolve

	with fits.open(image_path) as hdul:
		header = hdul[0].header.copy()
		data = np.squeeze(np.asarray(hdul[0].data, dtype=float))

	sigma_pixels = get_gaussian_sigma_pixels(image_path, smooth_fwhm)
	kernel_size = int(np.ceil(sigma_pixels * 5))
	if kernel_size % 2 == 0:
		kernel_size += 1
	kernel = Gaussian2DKernel(sigma_pixels, x_size=kernel_size, y_size=kernel_size)
	smoothed = convolve(data, kernel, boundary="extend")
	return fits.PrimaryHDU(data=smoothed.astype(np.float32), header=header)