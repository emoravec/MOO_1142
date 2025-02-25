XMM image log

XMM_adapt-400-7200.fits - adaptively smoothed image - don't use for analysis, only for plotting - from Craig
XMM_comb-net-center_incorrect.fits - same as XMM_comb-net-image but contains nans instead of 0 so creates issues with GGM filtering
XMM_comb-net-center.fits - the same as XMM_comb-net-image.fits but cutout and centered on cluster - from Craig March 5/6th, 2023 - use this one for (ggm) filtering!
XMM_comb-net-image.fits - the full XMM image - but has nans instead of 0s so causes large regions of white with GGM
XMM_comb-obj-im-400-7200.fits - from Tony, from Craig - I believe this is the exact same as XMM_comb-net-center.fits

XMM-comb-net-cen-M2.fits -  XMM image convolved to ~M2 resolution in fits and as a pdf. Made and sent by Craig on 2025-02-24.

XMM-2025-02: XMM images from Craig. I downloaded these from slack on 2025-02-25. I renamed his folder `MOO_XMM' to ``XMM-2025-02`. 
	- `comb-net-ctr.fits` is a zoomed in version of the `XMM_comb-net-image.fits`. `comb-net-ctr.fits` does have the nan regions but they are less apparent though somehow? `comb-net-ctr.fits` is similar to `XMM_comb-net-image.fits` but they are trimmed differently.
	- From DS9 looks like `adapt-ctr.fits` is the adaptively smoothed, centered image (`comb-net-ctr.fits`). 