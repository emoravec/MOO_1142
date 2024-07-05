import aplpy
import matplotlib
import matplotlib.pyplot as pyplot
import numpy as np
import time
# -------------------------------------------------------------------------------------------- #
# colormap parameters blue-orange colormap
import matplotlib.cm as cm

bottom = cm.get_cmap('Oranges', 128)
top = cm.get_cmap('Blues_r', 128)
newcolors = np.vstack((top(np.linspace(0, 1, 128)),bottom(np.linspace(0, 1, 128))))

from matplotlib.colors import ListedColormap
cm.register_cmap('OrangeBlue', cmap = ListedColormap(newcolors))
cmap = 'OrangeBlue'

# vmin=-3
# vmax=+3
# ------------------ #
def plot_minkasi_maps(N,M,tod_class,save_path):
    for n in N:
        for m in M:
            n_m_file_name = 'MOO_1142_'+tod_class+'_niter_'+str(n)+'_'+str(m)
            print(n_m_file_name)
            cluster = tod_class_path + n_m_file_name +'.fits'
            # ------------------ #
            fig = pyplot.figure(figsize=(5, 5))
            img = aplpy.FITSFigure(cluster, hdu=0, figure=fig, downsample=1, smooth=False, convention='calabretta')
            img.set_theme('publication')
            img.set_title(n_m_file_name)

            img.show_colorscale(cmap=cmap,stretch='linear',smooth=1)#,vmin=vmin,vmax=vmax - pixels in the minkasi maps are 2-3"/pixels so a smoothing of smooth=2-4 times that roughly gives the M2 beam of 9" - 3 seems like a lot and 1 seems like not enough?

            img.recenter(moo1142_ra, moo1142_dec, radius=2.0/60.0)
            img.ax.tick_params(axis='both',which='both',direction='in')
            img.add_scalebar(0.5/60.0, '30\"', color='black')

            img.add_beam(major=9.0/3600.0, minor=9.0/3600.0, angle=0)
            img.beam.set_color('white')
            img.beam.set_edgecolor('green')
            img.beam.set_facecolor('white')
            img.beam.set_corner('bottom left')

            img.add_colorbar('right')
            img.colorbar.set_width(0.12)
            img.colorbar.set_axis_label_text('K$_{rj}$') # minkasi maps are usually brightness temperature, K_rj (Tony M July 3rd, 2024 Slack)

            # clevels=[-250e-6, -200e-6, -150e-6, -100e-6, -50e-6, 0, 50e-6, 100e-6, 150e-6, 200e-6, 250e-6]
            # img.show_contour(cluster,colors="gray",levels=clevels, returnlevels=True,convention='calabretta',smooth=3)

            img.save(save_path + 'pdfs/' + n_m_file_name +'.pdf', format='pdf', dpi=300)
# ------------------ #
location = '/Users/emoravec/Documents/Research/merging_clusters/analysis/MOO_1142/M2/images/minkasi/'
moo1142_ra, moo1142_dec = np.rad2deg([3.06642436, 0.26972541])
tod_class = 'all'
tod_class_path = location + tod_class + '/'
N_array = [1,2,3,4,5]
M_array = [1,5,15,25,50]
# ------------------ #
# Record the start time
start_time = time.time()

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
