import imageio.v2 as imageio
import glob
from pathlib import Path

location = Path(__file__).resolve().parents[2] / 'ggm'
files = sorted(glob.glob(str(location / 'M2/pngs/moo1142_M2.ggm.*.png'), recursive = True))
gif_images = []
for f in files:
    gif_images.append(imageio.imread(f))
imageio.mimsave(location / 'M2/pngs/moo1142_M2.ggm.gif',gif_images,duration=1)