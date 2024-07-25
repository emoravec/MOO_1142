Explanation of files for members of MOO 1142

**Photometric redshifts** -> photz folder
make_photz_fits.py - makes a fits file from Mark's photz text files (specrank.MOO1142.wide_nohst.rank.txt and specrank.MOO1142.hst_core.rank.txt)

Files that are from Adam on August 29th, 2023 that show which phot-sz he has been using
- I’ve been using a few phot-z compilations that Mark made, which do appear to be somewhat different from your set - mainly the attached where the filenames are more or less self-explanatory.
	- m1142_HST_photz_12jan23.lst
	- m1142.wideall.photz.20jan23.lst

Files from Mark Brodwin on March 30th, 2023
- The catalog I’m referring to is one I created years ago for Perlmutter to spectroscopically target cluster galaxies for his SN program.  I’m attaching two catalogs.  The first is just in the HST region in the cluster core, but has the best photo-z’s.  The second is a wider catalog using ground-based data that lets you go further out in radius.  The ranking was based on integrated redshift probability, but also weighted somewhat toward brighter galaxies to make the observations easier.  It should be fine for your selection.  Just take object from the top of the lists.
	- specrank.MOO1142.hst_core.rank.txt (originally called specrank.MOO1142.rank)
	- specrank.MOO1142.wide_nohst.rank.txt (originally called specrank.MOO1142_wide_nohst.rank)

The following files were made from these specrank files from Mark
	- MOO_1142.photz_hst_core_IntPz0.3* -> hst core with IntPz > 0.3 
	- MOO_1142.photz_wide_2am* -> hst core within 2am of cluster center
	- MOO_1142.photz_wide_IntPz* -> hst core with the IntPz >= to that in file name
	- MOO_1142.photz_wide_IntPz0.6_ALMA.fits -> used for ALMA CO proposal submitted in Spring 2023

Files of the form *_galdens*.fits -> smoothed member distributions using the galdens code adapted from Luca.
Files of the form *photz*.fits without the _galdens are the catalogs of the members. 

MOO_1142+1527.galdens_smoothed_Luca.fits -> not sure - what Luca used in Dicker2020?
MOO_1142+1527.MC_filtered_members.reg -> members that are identified from the Gonzalez+19 search (color-color analysis) which are in MOO_1142+1527.spitzer_ps.cat_filtered.out

**Spectroscopic redshifts** -> specz folder

	* MOO_1142+1527.spec_members.15mar23.txt -> spectroscopic members from Adam Stanford sent on March 15th, 2023 "Here is the current list which includes the lower quality redshifts" -> originally was called m1142p1527.members.15mar23.lst

	* m1142p1527.redshifts.28aug23.txt  -> spectroscopic members from Adam Stanford sent on August 28th, 2023 -> Attached is a summary of all of our redshifts on m1142.  I think there are ~20 spec-z members already before including ones from the two masks observed in Dec 2022 and Jan 2023.  My initial reductions of those masks could be improved so not sure about the final yields but I think there are about 10 new spec-z members in those two masks.  We did not get anything from the masks designed for the March 2023 nights which were lost mainly because of weather.
		* Emily Moravec removed anything that was not marked A-C to make MOO_1142+1527.spec_members.28aug23.txt

	* Use /Users/emoravec/Documents/Research/merging_clusters/analysis/vla21b/cluster_members/make_ds9_region_members.py to convert to fits file.

**color-color members** -> color-color folder
Files in this folder were made from color-color analysis

	* MOO_1142+1527.box_cluster_members.fits -> those within the box in the color-color diagram that I defined in the vla21b color-color analysis that chose members -> I think that these form the MOO_1142+1527.Moravec_members.reg file
	MOO_1142+1527.box_members.galdens.fits -> the above but run through Luca's galdens algorithm

	* MOO_1142+1527.spitzer_ps.cat_filtered.out -> I *think* this is the catalog of members that Bandon made which was a recreatinon of the color color cuts used to make the original MaDCoWS cluster member catalogs
		* I *think* this goes with MOO_1142+1527.MC_filtered_members.reg where I created a ds9 region file of MOO_1142+1527.spitzer_ps.cat_filtered.out?

	* MOO_1142+1527.Moravec_members.reg -> unsure what this is without looking deeper into it


