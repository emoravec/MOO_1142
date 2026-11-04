Photometric redshifts

HST
	* make_photz_fits.py - makes a fits file from Mark's photz text files (specrank.MOO1142.wide_nohst.rank.txt and specrank.MOO1142.hst_core.rank.txt)

	* Files from Mark Brodwin on March 30th, 2023
		- The catalog I’m referring to is one I created years ago for Perlmutter to spectroscopically target cluster galaxies for his SN program.  I’m attaching two catalogs.  The first is just in the HST region in the cluster core, but has the best photo-z’s.  The second is a wider catalog using ground-based data that lets you go further out in radius.  The ranking was based on integrated redshift probability, but also weighted somewhat toward brighter galaxies to make the observations easier.  It should be fine for your selection.  Just take object from the top of the lists.
			- specrank.MOO1142.hst_core.rank.txt (originally called specrank.MOO1142.rank)
			- specrank.MOO1142.wide_nohst.rank.txt (originally called specrank.MOO1142_wide_nohst.rank)

	* The following files were made from these specrank files from Mark
		- MOO_1142.photz_hst_core_IntPz0.3* -> hst core with IntPz > 0.3 
		- MOO_1142.photz_wide_2am* -> hst core within 2am of cluster center
		- MOO_1142.photz_wide_IntPz* -> hst core with the IntPz >= to that in file name
		- MOO_1142.photz_wide_IntPz0.6_ALMA.fits -> used for ALMA CO proposal submitted in Spring 2023

		* Files that are from Adam on August 29th, 2023 that show which phot-sz he has been using
			- I’ve been using a few phot-z compilations that Mark made, which do appear to be somewhat different from your set - mainly the attached where the filenames are more or less self-explanatory.
				- m1142_HST_photz_12jan23.lst
				- m1142.wideall.photz.20jan23.lst

MC1_cat - files from the photz analysis in the MC1 catalog paper
	* catalogs
		- MOO_1142+1527.fits: the final photz catalog from the MC1 catalog work (/Users/emoravec/Documents/Research/MaDCoWS/catalogs/photz/from_mark/final/MOO_1142+1527.fits)
	* galdens
		- MOO_1142+1527_final_full_IntPzdens.fits: the final integrated probability surface density from the MC1 catalog work (/Users/emoravec/Documents/Research/MaDCoWS/catalogs/photz/density_centroids/galdens_fits/final_full_map/MOO_1142+1346_IntPzdens.fits). This is using all IntPzs so includes all the information. Used a sigma=30" and a pixel size of 2" (what we decided on for the MC1 catalog work).
	* old 
		* catalogs - Files of the form *photz*.fits without the _galdens are the catalogs of the members. 
		* galdens - files of the form *_galdens*.fits -> smoothed member distributions using the galdens code adapted from Luca.
			- MOO_1142+1527.galdens_smoothed_Luca.fits -> not sure - what Luca used in Dicker2020?

