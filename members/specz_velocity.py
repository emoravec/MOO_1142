"""
Plots the velocity distributions of cluster members.

Created in 2025-10
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path

from astropy import constants as const
from astropy.table import Table
# -------------------------------------------------------------------------------------------- #
location = Path('/Users/emoravec/Documents/Research/merging_clusters/analysis/MOO_1142/members')
def calculate_velocity(z):
    """Calculate the velocity corresponding to a given redshift."""
    c = const.c.to('km/s').value  # Speed of light in km/s
    velocity = c * (((1+z)**2 - 1) / ((1 + z)**2+1))
    return velocity

def calculate_delta_velocity(z_galaxy):
    """Calculate the velocity of a galaxy relative to the cluster redshift."""
    v_MOO_1142 = calculate_velocity(1.189)
    v_galaxy = calculate_velocity(z_galaxy)
    delta_velocity = v_MOO_1142 - v_galaxy
    return np.rint(delta_velocity).astype(int)
# -------------------------------------------------------------------------------------------- #
# open pyrostat members file
with open(location / 'specz/MOO_1142+1527.speczs.goodQ.pyrostat.members.14aug25.txt', 'r') as f:
    first_line = f.readline().strip()
    column_names = first_line.replace('#', '').strip().split()
# Read the data, skipping the comment line
specz_members = pd.read_csv(location / 'specz/MOO_1142+1527.speczs.goodQ.pyrostat.members.14aug25.txt', 
                           delim_whitespace=True, 
                           comment='#',
                           names=column_names)
# calculate delta velocities
delta_v = calculate_delta_velocity(specz_members['Spec-z'].values)
specz_members['Delta_v'] = delta_v
# # save updated members file with delta velocities
# specz_members.to_csv(location / 'specz/MOO_1142+1527.speczs.goodQ.pyrostat.members.with_delta_v.25oct25.csv', index=False)

# Create histogram of Delta_v
plt.figure(figsize=(10, 6))
#bins = np.arange(-3000, 3500, 500)  # -3000 to +3000 in increments of 500
bins = np.arange(-3000, 3200, 200)  # -3000 to +3000 in increments of 200
plt.hist(specz_members['Delta_v'], bins=bins, alpha=0.7, color='skyblue', edgecolor='black')
plt.xlabel('$\Delta_v$ (km/s)')
plt.ylabel('N')
plt.title('MOO 1142+1527 Spec-z Members')
plt.grid(True, alpha=0.3)
plt.tight_layout()
# Save the plot
plt.savefig(location / 'plots/velocity_distribution/MOO_1142_specz_members_delta_v_hist_bins=200.pdf', format='pdf', dpi=300, bbox_inches='tight')
