2023-03
	* 2025-03-15.MOO_1142+1527.spec_members.txt -> spectroscopic members from Adam Stanford sent on March 15th, 2023 "Here is the current list which includes the lower quality redshifts" -> originally was called m1142p1527.members.15mar23.lst
	* MOVED THESE TO THE "OLD" FOLDER

2023-08
	* 2023-08-28.MOO_1142+1527.redshifts.txt (originally called m1142p1527.redshifts.28aug23.txt)  -> spectroscopic members from Adam Stanford sent on August 28th, 2023 -> Attached is a summary of all of our redshifts on m1142.  I think there are ~20 spec-z members already before including ones from the two masks observed in Dec 2022 and Jan 2023.  My initial reductions of those masks could be improved so not sure about the final yields but I think there are about 10 new spec-z members in those two masks.  We did not get anything from the masks designed for the March 2023 nights which were lost mainly because of weather.
	* 2023-08-28.MOO_1142+1527.speczs.txt - Emily Moravec removed anything from 2023-08-28.MOO_1142+1527.redshifts.txt was not marked A-C to make 
		* 2025-07-29: Somehow 2023-08-28.MOO_1142+1527.speczs.txt doesn't have all A-C from 2023-08-28.MOO_1142+1527.redshifts.txt? Not sure if I applied some other cuts somehow?
		* 2025-07-29: Emily Moravec took 2023-08-28.MOO_1142+1527.redshifts.txt and applied a quality filter to only have sources with Q == 'A-D' to 2023-08-28.MOO_1142+1527.redshifts.goodQ.txt 
	* MOVED THESE TO THE "OLD" FOLDER

* 2025-07
	* 2025-07-30.MOO_1142+1527.internal.redshifts.lst was previously named M1142+1527.internal.redshifts-30jul25.txt
		* sent by Adam Stanford on 2025-07-30
		* most updated version of the speczs for MOO 1142 at the time. He said, "This was the compilation I made back in January - there is an extensive header outlining various issues that remain."
	* ran speczs/convert_txt_to_fits.py on 2025-07-30.MOO_1142+1527.internal.redshifts.lst to create output .fits and .csv
		* 2025-07-30.MOO_1142+1527.speczs.fits is just the fits version of 2025-07-30.MOO_1142+1527.internal.redshifts.lst
		* 2025-07-30.MOO_1142+1527.speczs.goodQ.csv and .fits are only the quality 1, 2, and 3 redshifts from 2025-07-30.MOO_1142+1527.internal.redshifts.lst
	* 2025-07-30: Adam Stanford said that, "I’d suggest using 0.003 for the uncertainties."

* 2025-08
	* 2025-08-14.MOO_1142+1527.speczs.goodQ.pyrostat.members.14aug25.txt. Using convert_txt_to_fits.py, I wrote out the speczs with good quality redshifts (quality 3-1) from 2025-07-30.MOO_1142+1527.internal.redshifts.lst to send to Anthony. He then ran pyrostat on this file. Explanation of pyrostat: "the code basically does an iterative sigma-clipping to compute the velocity dispersion and keep galaxies that are within 3 sigma of the mean velocity. " After running pyrostat, he sent this list back to me and said "The biweight redshift is 1.1934 with 36 members."

* 2025-10
	* Ran /Users/emoravec/Documents/Research/MOO_1142/members/specz_velocity.py on 2025-08-14.MOO_1142+1527.speczs.goodQ.pyrostat.members.14aug25.txt to calculate the velocity of the sources.

* 2026-06
	* There had been progress on the MC1 catalog paper so I now have a list of the specz redshifts that we are including in the MC1 paper. Additionally there were some errors that Mark found in the catalog so he worked with Adam to correct those to produce /Users/emoravec/Documents/Research/MaDCoWS/catalogs/specz/MC1_specz_programs/All.public.redshifts.29jun26.lst (has coordinate erros fixed but not the coordinates to the exact same number of sig figs). I took All.public.redshifts.29jun26.lst and selected only the sources that were for MOO 1142 and that produced 2026-06-29.MOO_1142+1527.speczs.from.all.public.redshifts.fits. I gave 2026-06-29.MOO_1142+1527.speczs.from.all.public.redshifts.fits to Delaney White to combine with the new speczs that she calculated.

* 2026-10
	* Keck ground-based spectroscopy
		* extra/2026-10-07.MOO_1142+1527.speczs.DW_combined.fits = 2026-06-29.MOO_1142+1527.speczs.from.all.public.redshifts.fits + the new redshifts that Delaney calculated.
			* Made by Delaney White in october 2026. This has all the MOO 1142 redshifts from Adam and new redshifts from Delaney. The new redshifts are the ones with "source" == "Keck LRIS red". The rest are from the specz table going in the MC1 paper which I think was /Users/emoravec/Documents/Research/MaDCoWS/catalogs/specz/MC1_specz_programs/All.public.redshifts.29jun26.lst
			* She said, "I have not removed any potential duplicates because I wasn't sure what threshold you all were interested in using. When I did a search for duplicates using a search radius of 1.5" it suggested 12 duplicates. Most seem to have redshifts in agreement, but those that don't generally have lower Q of 1 or 2 for both catalogs, so they may need to be removed when moving forward."
			* moved this to the extra folder because not particularly useful as it is
		* 2026-10-07.MOO_1142+1527.specz_dups.fits - Delaney cross matched all redshifts she had within 1.5" and got 12 matches = duplicates
		* 2026-10-08.MOO_1142+1527.speczs.DW_all.8oct26.fits
			* Made by Delaney White in october 2026. She reduced three new MOO 1142 masks and produced 43 redshifts. This is the full list of those redshifts with some extra columns I asked for (wanted it to be in the standard format that we plan to in the MC1 paper). 
	* JWST
		* 2026-10-08.MOO_1142+1527.speczs.JWST.fits - list of JWST speczs from the MOO 1142 JWST program led by Yen-Ting Lin. Anthony asked Wren Heiligenstein to send them to me from the paper.
		* 2026-10-15.MOO_1142+1527.speczs.JWST.grism.lst and 2026-10-15.MOO_1142+1527.speczs.JWST.prism.lst - lists of individual spectra from JWST grism and prism respectively. This has the full quality information for the sources that ended up in 2026-10-08.MOO_1142+1527.speczs.JWST.fits.

2026-10-14: renamed files in this folder so that they were ordered based on date (changed format of file from MOOXXXX.ddMMMYY.fits to YYYY-MM-DD.MOO_1142+1527.fits)

