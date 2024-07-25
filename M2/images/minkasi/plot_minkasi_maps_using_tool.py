"""
Makes plots of minkasi maps using witch plotting tool written by Jack.

Created in 2024-07
"""
import numpy as np
import time
from witch.plotting import plot_cluster
print('Imports complete.')
# ------------------ #
moo1142_ra, moo1142_dec = np.rad2deg([3.06642436, 0.26972541])

def plot_minkasi_maps(N,M,tod_class,save_path):
    for n in N:
        for m in M:
            n_m_file_name = 'MOO_1142_'+tod_class+'_niter_'+str(n)+'_'+str(m)
            print(n_m_file_name)
            cluster_fits_path = tod_class_path + n_m_file_name +'.fits'
            img = plot_cluster("MOOJ1142", 
                               cluster_fits_path, 
                               ra=moo1142_ra, dec=moo1142_dec, 
                               pix_size = 3.0, 
                               radius = 3.0, 
                               smooth = 3, 
                               ncontours=21, 
                               plot_r = False, 
                               units = "uK_cmb")
            img.save(save_path + 'pdfs/witch_tool/' + n_m_file_name +'.pdf', format='pdf', dpi=300)

# ------------------ #
location = '/Users/emoravec/Documents/Research/merging_clusters/analysis/MOO_1142/M2/images/minkasi/'
tod_class = 'good'
tod_class_path = location + tod_class + '/'
N_array = [1,2,3,4,5]
M_array = [1,5,15,25,50]
# ------------------ #
# Record the start time
start_time = time.time()

# do plotting
plot_minkasi_maps(N_array,M_array,tod_class,tod_class_path)

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
