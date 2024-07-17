import numpy as np
#from matplotlib import pyplot as plt
import minkasi.minkasi_all as minkasi

import time
import glob
import os
#reload(minkasi)
#plt.ion()

# Record the start time
start_time = time.time()

#set file root for output maps
#outroot = "/Users/emoravec/Documents/Research/merging_clusters/analysis/MOO_1142/M2/images/minkasi/test/PREFIX" 
#outroot = "/Users/emoravec/Documents/Research/merging_clusters/analysis/MOO_1142/M2/images/minkasi/good/MOO_1142_good" 
#outroot = "/Users/emoravec/Documents/Research/merging_clusters/analysis/MOO_1142/M2/images/minkasi/remove_no/MOO_1142_maybe_good" 
outroot = "/Users/emoravec/Documents/Research/merging_clusters/analysis/MOO_1142/M2/images/minkasi/all/MOO_1142_all" 

#find tod files we want to map
idir = "/Users/emoravec/Documents/Research/merging_clusters/analysis/MOO_1142/M2/TODs/data/" 
tod_names=glob.glob(idir+'/Sig*.fits') # all
#tod_names=glob.glob(idir+'/Signal*23B_005_01*.fits')
if len(tod_names)==0:
    print('We found no TOD files.  Double check your path?')
    assert(1==0)

# #Use presets by source to automatically get and cut TODs
# #that were manually flagged for removal
# bad_tod, _ = minkasi.get_bad_tods("moo1142")
# tod_names = minkasi.cut_blacklist(tod_names, bad_tod)
# print('The Number of TODs being USED is:',len(tod_names))

#if running MPI, you would want to split up files between processes
#one easy way is to say to this:
tod_names=tod_names[minkasi.myrank::minkasi.nproc]
#NB - minkasi checks to see if MPI is around, if not
#it sets rank to 0 an nproc to 1, so this would still
#run in a non-MPI environment

todvec=minkasi.TodVec()

#loop over each file, and read it.
for fname in tod_names:
    t1=time.time()
    dat=minkasi.read_tod_from_fits(fname)
    t2=time.time()
    minkasi.truncate_tod(dat) #truncate_tod chops samples from the end to make
                              #the length happy for ffts
    minkasi.downsample_tod(dat)   #sometimes we have faster sampled data than we need.
                                  #this fixes that.  You don't need to, though.
    minkasi.truncate_tod(dat)    #since our length changed, make sure we have a happy length

    #figure out a guess at common mode #and (assumed) linear detector drifts/offset
    #drifts/offsets are removed, which is important for mode finding.  CM is *not* removed.
    dd=minkasi.fit_cm_plus_poly(dat['dat_calib'])  
    
    dat['dat_calib']=dd
    t3=time.time()
    tod=minkasi.Tod(dat)
    todvec.add_tod(tod)
    print('took ',t2-t1,' ',t3-t2,' seconds to read and downsample file ',fname)

#make a template map with desired pixel size an limits that cover the data
#todvec.lims() is MPI-aware and will return global limits, not just
#the ones from private TODs
lims=todvec.lims()
pixsize=2.0/3600*np.pi/180
map=minkasi.SkyMap(lims,pixsize)

#once we have a map, we can figure out the pixellization of the data.  Save that
#so we don't have to recompute.  Also calculate a noise model.  The one here
#(and currently the only supported one) is to rotate the data into SVD space, then
#smooth the power spectrum of each mode.  Other models would not be hard
#to implement.  The smoothing could do with a bit of tuning as well.
for tod in todvec.tods:
    ipix=map.get_pix(tod)
    tod.info['ipix']=ipix
    tod.set_noise(minkasi.NoiseSmoothedSVD) 

#get the hit count map.  We use this as a preconditioner
#which helps small-scale convergence quite a bit.
hits=minkasi.make_hits(todvec,map)

#setup the mapset.  In general this can have many things
#in addition to map(s) of the sky, but for now we'll just 
#use a single skymap.
mapset=minkasi.Mapset()
mapset.add_map(map)

#make A^T N^1 d.  TODs need to understand what to do with maps
#but maps don't necessarily need to understand what to do with TODs, 
#hence putting make_rhs in the vector of TODs. 
#Again, make_rhs is MPI-aware, so this should do the right thing
#if you run with many processes.

rhs=mapset.copy()
todvec.make_rhs(rhs)

#this is our starting guess.  Default to starting at 0,
#but you could start with a better guess if you have one.
x0=rhs.copy()
x0.clear()

#preconditioner is 1/ hit count map.  helps a lot for 
#convergence.
precon=mapset.copy()
tmp=hits.map.copy()
ii=tmp>0
tmp[ii]=1.0/tmp[ii]
precon.maps[0].map[:]=tmp[:]

#run PCG!

save_iters=[1,5,15,20,25,50]

#run PCG!
mapset_out=minkasi.run_pcg(rhs,x0,todvec,precon,maxiter=50, save_iters=save_iters, outroot = outroot)

if minkasi.myrank==0:
    mapset_out.maps[0].write(outroot+'_initial.fits') #and write out the map as a FITS file
else:
    print('not writing map on process ',minkasi.myrank)

#This first pass of fitting is a reasonable approximation. However, the SVD noise estimation we use invariable includes some signal in the noise model.
#In other words, our noise model assumes our TODs are 100% noise, performs an SVD, and uses the first N modes of that SVD as the noise estime. However,
#some portion of those modes are actually signal. To counteract this, we take the initial best-fit map (which is our best estimate of the signal), 
#subtract it from a copy of the TODs, and reestimate the noise. Subtracting the best-fit map should remove the largest signal modes, making our estimate of the
#noise cleaner. We can then remake the map and itterate on this process. Note the best-fit map is not subtracted from the TODs themselves for the purposes
#of mapmaking. That would remove signal from the final map. It is only removed from a copy of the TODs (technically a copy of the data in the TODs) used
#in noise estimation.

npass = 5
for niter in range(npass):
    maxiter = 26 + 25 * (niter + 1)
    # first, re-do the noise with the current best-guess map
    for tod in todvec.tods:
        mat = 0 * tod.info["dat_calib"]
        for mm in mapset_out.maps:
            mm.map2tod(tod, mat)
        tod.set_noise(
            minkasi.NoiseSmoothedSVD, dat=tod.info["dat_calib"] - mat, #note here we are feeding tod.info["dat_calib] - mat to dat, as opposed to
                                                                       #say copying tod.info["dat_calib"], subtracting mat, and passing it. 
                                                                       #This is the safest way to avoid changing the values stored in tod.info
        )
 
    priorset = None

    rhs = mapset.copy() 
    todvec.make_rhs(rhs)

    precon = mapset.copy() #Hist map preconditioner
    precon.maps[0].map[:] = hits.map[:]
    mapset_out = minkasi.run_pcg_wprior(
        rhs,
        mapset,
        todvec,
        priorset,
        precon,
        maxiter=maxiter,
        outroot= str(outroot) + "_niter_" + repr(niter + 1),
        save_iters=save_iters,
    )
    if minkasi.myrank == 0:
        mapset_out.maps[0].write(
            str(outroot) + "_niter_" + repr(niter + 1) + ".fits"
        )

minkasi.barrier()

if minkasi.myrank==0:
    print("Maps written to :", outroot)

# Record the end time
end_time = time.time()

# Calculate the elapsed time
elapsed_time = end_time - start_time

# Convert the elapsed time to hours and minutes
hours = int(elapsed_time // 3600)
minutes = int((elapsed_time % 3600) // 60)
seconds = int(elapsed_time % 60)

# Print the elapsed time in hours and minutes
print(f"Elapsed time: {hours} hours, {minutes} minutes, and {seconds} seconds")

