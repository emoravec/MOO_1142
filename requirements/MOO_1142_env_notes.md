April 2026

2026-04-17: when I created a new MOO_1142 conda env based on what I knew I needed from the merging_clusters env
I did the following first

conda create -n MOO_1142 python=3.12 scipy pandas matplotlib pip ipdb ipython ipykernel -y
conda activate MOO_1142
conda install -c conda-forge astropy=5 dill corner
pip install aplpy==2.2.0 pocomc

Then started to run scripts to see if plotting worked. There was some sort of error with importing aplpy. The robot told me I could either downgrade matplotlib or upgrade astropy to astropy>=7. I chose to upgrade astropy. I also realized I needed palettable for plotting the density distributions. Created the env files.

I will note that it looks like one of these packages requires pyregion so might have issues once I try to use that (if I do). See note in Obsidian titled "Installing aplpy in a new environment"

2026-04-16 - when I split /Users/emoravec/Documents/Research/merging_clusters/MOO_1142 into its own director and git repo at /Users/emoravec/Documents/Research/MOO_1142/