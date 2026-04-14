import imageio.v2 as imageio
import glob

location ='/Users/emoravec/Documents/Research/merging_clusters/analysis/MOO_1142/ggm'
files = sorted(glob.glob(location + '/M2/pngs/moo1142_M2.ggm.*.png', recursive = True))
gif_images = []
for f in files:
    gif_images.append(imageio.imread(f))
imageio.mimsave(location + '/M2/pngs/moo1142_M2.ggm.gif',gif_images,duration=1)