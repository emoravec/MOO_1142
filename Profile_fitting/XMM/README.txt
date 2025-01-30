Xpeak_circle_nosub.dat - The XMM surface brightness profile for MOO 1142 made by Craig Sarazin in early 2025. See email from Craig on 2025-01-05 entitled "First XMM X-ray profile".

Columns for Xpeak_circle_nosub.dat
rm, ri, ro, SB, rate, rate_err, SB_err

where
rm = mean radius of annulus (“)
ri = inner radius (“)
ro = outer radius (“)
SB = total XMM surface brightness (counts/s/arcsec**2)
[rate = total count rate in counts/s]
[rate_err = uncertainty in rate]
SB_err = uncertainty in SB

I think you will only need columns 1-4 and 7.   I need the rates for some of my software.

ri and and ro - I wouldn’t say “error”.  ri and ro are the inner and outer radii of the annulus used to extract the surface brightness, and are not “uncertain”.   But, I put a |--|  bracketin the plots to show the radii covered by the annulus. rm is just the average of ri and ro, so in that sense they are symmetric.