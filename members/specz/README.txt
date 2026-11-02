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
		* 2026-10-08.MOO_1142+1527.speczs.DW_all.8oct26.fits
			* Made by Delaney White in october 2026. She reduced three new MOO 1142 masks and produced 43 redshifts. This is the full list of those redshifts with some extra columns I asked for (wanted it to be in the standard format that we plan to in the MC1 paper). 
		* 2026-10-19.MOO_1142+1527.speczs.internal.groundBasedAll.fits
			* Made by Delaney White. Table of 124 sources that is all the ground based spectroscopy - Adams internal list + all new Keck spectroscopy that Delaney reduced. This has the duplicates resolved. 
			* There were 13 "duplicates" - 13 original internal sources that had a match with 12 of Delaneys new catalog of 43 sources. Two "original" sources were matched with one of Delaney's new sources. 
			* The general procedure of resolving duplicates is as follows: We identified duplicates across catalogs by using a search radius of 1.5”. With this matching criteria, 13 sources from the original catalog matched to 12 sources in the new catalog. This suggests that two sources in the previous catalog were duplicates already as they reported similar redshift values. To determine which case to keep, we compared redshift values and the quality flag, “Q”, assigned to each source. For each duplicate case, we kept the redshift value that corresponded to the estimate that had a higher Q value. If cases had the same Q value, the latest estimate for the source redshift was used. We also estimated duplicates to have a significant redshift difference when $\delta z > 0.01 (1+z)$. Out of the 13 duplicate sources, 11 had nonsignificant redshift differences. For the two sources with significantly different redshifts, one case showed a higher Q value with the latest estimate, while the other had Q=1 for both estimates and would ultimately be removed from the catalog.
			* I made this from extra/2026-10-19.MOO_1142+1527.speczs.internal.groundBasedAll.fits and edited the columns to be ready to combine with JWST (e.g., to have an ID column)
		* 2026-10-22.MOO_1142+1527.speczs.internal.groundBasedDups.fits
			* Made by Delaney White. A table that reports the duplicates that were found between the original internal catalog and her new redshifts. She said that Ireport z_final and Q_final columns that capture the final redshift and Q taken between the duplicates. I also included a column here "DUP_TYPE" that describes if the source would be considered a cluster member (CM) or non-cluster member (NCM) using a velocity dispersion of 10,000 km/s and if the redshift difference was significant (SIG_DIFF) or insignificant (NO_SIG_DIFF) given the >0.01(1+z) criteria Adam shared.
			* Delaney also reports the separation between the original specz and the new specz.
			* Delaney kept the final duplicate table limited to the 12 matches. For the case that had two matches from the original catalog, the columns that report values for the original catalog list both values for these two matches and the z_final and Q_final column list the value that was taken given the process I described previously for resolving duplicates.  
			
	* JWST
		* 2026-10-08.MOO_1142+1527.speczs.JWST.csv - list of JWST speczs from the MOO 1142 JWST program led by Yen-Ting Lin. Anthony asked Wren Heiligenstein to send them to me from the paper.
		* 2026-10-08.MOO_1142+1527.speczs.JWST.fits - is a standardized version of 2026-10-08.MOO_1142+1527.speczs.JWST.csv to prepare for combination with the ground based spectroscopy.
		* 2026-10-15.MOO_1142+1527.speczs.JWST.grism.lst and 2026-10-15.MOO_1142+1527.speczs.JWST.prism.lst - lists of individual spectra from JWST grism and prism respectively. This has the full quality information for the sources that ended up in 2026-10-08.MOO_1142+1527.speczs.JWST.fits.

* 2026-11
	* cross matching ground with JWST - see cross_match_ground_JWST
		* ground_JWST_1as.fits, ground_JWST_2as.fits, ground_JWST_2p5as.fits - ground and JWST duplicates within 1, 2, and 2.5" (7, 8 and 9 matches).
	* 2026-11.02.MOO_1142+1527.speczs.all.fits - all ground based + JWST redshifts with duplicates resolved

2026-10-14: renamed files in this folder so that they were ordered based on date (changed format of file from MOOXXXX.ddMMMYY.fits to YYYY-MM-DD.MOO_1142+1527.fits)

Facility codes: ~two letter code that denotes the first letter of the telescope and then the instrument that were used to take the spectrum. 
GN = Gemini North with GMOS
GO = Gran Telescopio Canarias OSIRIS, 
JG =  JWST grism
JP = JWST prism
KD = Keck DEIMOS, 
KG means that the spectrum was taken during multiple Keck and Gemini observing runs with either LRIS or DEIMOS at Keck, and GMOS at Gemini, 
KL = Keck LRIS, 
KLr = Keck LRIS red - Given the high redshift of this cluster, only data from the LRIS red side camera (0.3 - 1.0 microns) were used and reduced 
KM = Keck MOSFIRE, 
PD = Palomar DBSP, 
PN = Palomar NGPS
