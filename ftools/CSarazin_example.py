#!/usr/bin/python
###################################################
#
# Extracts images for photon energy slices for Stefano
#
# Requires that heainit have been run
# Requires that sasinit have been run
#
# Inputs (Standard input):
#
#	direct = directory for output files (must already exist
#
#	prefix = start of file names for input and output (e.g., mos2S002)
#
#	slices = names of the file with the energy slices.  Each file has
#
#		elow = lower limit energy in eV
#		ehigh = upper limit energy in eV
#		chanl = lower limit channel in spectrum
#		chanh = upper limit channel in spectrum
#		rowl = lower limit row in spectral files
#		rowh = upper limit row in spectral files
#
# Output:
#
#	file bins.dat
#
#		rowi, rowf, bini, binf, ei, ef
#
#			rowi = starting row of bin
#			rowf = ending row of bin
#			bini = starting spectral bin
#			binf = ending spectral bin
#			ei = starting photon energy in eV
#			ef = ending photon energy in eV
#
#####################################################

import math
import subprocess
from subprocess import Popen, PIPE
import sys
import readline
from array import *
import os.path
# from astropy.io import fits

# read in output directory
print " "
print "Dirctory for final output (must already exit"
print " "

# read in directory
direct = str(raw_input("Directory for output: "))
# Make directory for output if needed
dirok = os.path.isdir(direct)
if not dirok:
     subprocess.call(["mkdir", direct])

# read in prefix for file names
print " "
print "Prefix for file names (e.g., mos2S002)"
print " "

# read in the prefix
prefix = str(raw_input("Prefix for file names: "))

print " "
print "File with energy slices for images"
print " "

# read in the file with energy slices
slices = prefix + "_bins.dat"

# Some needed parameters
# x & y size of output images
xsize = str(450)
ysize = str(450)

# x & y minimums and maximums values for image (physical coords)
xmin = str(21830)
xmax = str(30829)
ymin = str(23277)
ymax = str(32276)

# x & y bin size in arcsec
xbin = str(1)
ybin = str(1)

# name and extension for events file
events = "table=" + prefix + "-filt.fits:EVENTS"

# expressions for image binning
xbins = "ximagebinsize=" + xbin
ybins = "yimagebinsize=" + xbin

# expressions for image size
xsizes="ximagesize=" + xsize
ysizes="yimagesize=" + ysize
xmins = "ximagemin=" + xmin
xmaxs = "ximagemax=" + xmax
ymins = "yimagemin=" + ymin
ymaxs = "yimagemax=" + ymax

# name and extension for events file
events = "table=" + prefix + "-filt.fits:EVENTS"

# name of cleaned events file
clean = prefix + "-filt.fits"
cleans = "eventset=" + clean

# name for the mask file
mask = prefix + "_mask.fits[1]"

# location of AGN (RA, Dec, in degrees)
agnx = str(216.6341583)
agny = str(35.13763722)
agnxs = "x=" + agnx
agnys = "y=" + agny

# psf size
psfxs = "xsize=200"
psfys = "ysize=200"
# psf offsets
psfxo = str(131)
psfyo = str(121)
psfxos = "xoffset=" + psfxo
psfyos = "yoffset=" + psfyo

# file names for background spectra
pbgspec = prefix + "_pbg.dat"
xbgspec = prefix + "_xbg.dat"

# file name for pbg template
pbgtmpl = prefix + "_pbg_tmpl.fits"

# open slices file
fbin = open(slices,'r')

# start the big loop
for line in iter(fbin):
    r1, r2, r3, r4, r5, r6 = line.split()
    elow = str(r1)
    ehigh = str(r2)
    chanl = int(r3)
    chanh = int(r4)
    rowl = str(r5)
    rowh = str(r6)

# prefix for slice image names
    impre = prefix + "-" + elow + "-" + ehigh + "-"

# string for energy selection
#    esel = "expression='(PI in [" + elow + ":" + ehigh + "])'"
    esel = "expression=(PI in [" + elow + ":" + ehigh + "])"

# remove any old unmasked image files

    subprocess.call(["/bin/rm", "-f", "nomask.fits"])

# run evselect to create the data image with 1 arcsec pixels

    subprocess.call(["evselect", events, "withimageset=yes",
	"filtertype=expression", "xcolumn=X", "ycolumn=Y", xbins, ybins,
	"squarepixels=yes", xsizes, ysizes, "imagebinning=imageSize", xmins,
	xmaxs, "withxranges=yes", ymins, ymaxs, "withyranges=yes",
	"imagedatatype=Int32", "withimagedatatype=yes", "raimagecenter=0",
	"decimagecenter=0", "withcelestialcenter=no", esel,
	"imageset=nomask.fits"])

# strings for pi selection
    pil = "pimin=" + elow
    pih = "pimax=" + ehigh

# remove any old unmasked exp files

    subprocess.call(["/bin/rm", "-f", "exp_nomask.fits"])

# make the exposure map

    subprocess.call(["eexpmap", "imageset=nomask.fits",
	"attitudeset=atthk.fits", cleans, "expimageset=exp_nomask.fits", pil,pih])

# mask the image and exposure map

# name for data image
    image = impre + "image.fits"
    subprocess.call(["/bin/rm", "-f", image])

# make data image
    subprocess.call(["farith", "nomask.fits", mask, image, "MUL"])

# namne for exposure image
    exp = impre + "exp.fits"
    subprocess.call(["/bin/rm", "-f", exp])

# make exp image 
    subprocess.call(["farith", "exp_nomask.fits", mask, exp, "MUL"])

# mean energy for psf
    emean = str(int(0.5*(float(elow)+float(ehigh))))
    emeans = "energy=" + emean

# name for psf
    psf = impre + "psf.fits"
# name for masked psf
    psfm = impre + "psfm.fits"
# remove any old psfs
    subprocess.call(["/bin/rm", "-f", psf])
    subprocess.call(["/bin/rm", "-f", psfm])
    subprocess.call(["/bin/rm", "-f", "psf.fits"])
    subprocess.call(["/bin/rm", "-f", "psfm.fits"])


# make the psf
    subprocess.call(["psfgen", "image=nomask.fits", emeans, "level=ELLBETA",
	"coordtype=eqpos", agnxs, agnys, psfxs, psfys, "output=psf.fits"])
# put the psf in the blank image
    subprocess.call(["fimgmerge", "blank.fits", "psf.fits", psf, psfxos, psfyos])
# mask the psf
    subprocess.call(["farith", psf, mask, psfm, "MUL"])

# BACKGROUND

# remove any old spectrum files
    subprocess.call(["/bin/rm", "-f", "pbg.fits"])
    subprocess.call(["/bin/rm", "-f", "xbg.fits"])

# Particle background

    fpbg = open(pbgspec,'r')
    pbgct = float(0.0)
    for line2 in iter(fpbg):
	r1, r2, r3, r4, r5 = line2.split()
	row = str(r1)
	chan = int(r2)
	count = float(r3)
	elow = int(1000.0*float(r4))
	ehigh = int(1000.0*float(r5))
	if chan >= chanl and chan <= chanh:
	    pbgct = pbgct + count
    fpbg.close()
    subprocess.call(["fcarith", pbgtmpl, str(pbgct), "pbg.fits", "MUL"])

# X-ray BG

# total exposure in map

    subprocess.call(["fimgstat", exp, "INDEF", "INDEF"])
    imsum = float(subprocess.check_output(["pget", "fimgstat", "sum"]))

# total counts in X-ray background spectrum

    fxbg = open(xbgspec,'r')
    xbgct = float(0.0)
    for line3 in iter(fxbg):
	r1, r2, r3, r4, r5 = line3.split()
	row = str(r1)
	chan = int(r2)
	count = float(r3)
	elow = int(float(r4))
	ehigh = int(float(r5))
	if chan >= chanl and chan <= chanh:
	    xbgct = xbgct + count
    fxbg.close()
    ratio = xbgct/imsum
    subprocess.call(["fcarith", exp, str(ratio), "xbg.fits", "MUL"])

# total background

# name for background image
    back = impre + "back.fits"
    subprocess.call(["/bin/rm", "-f", back])

    subprocess.call(["farith", "xbg.fits", "pbg.fits", back, "ADD"])

# Move output to subdirectory

    subprocess.call(["mv", image, direct])
    subprocess.call(["mv", exp, direct])
    subprocess.call(["mv", psf, direct])
    subprocess.call(["mv", psfm, direct])
    subprocess.call(["mv", back, direct])

# close bin file
fbin.close()

print "Completed"

sys.exit()


