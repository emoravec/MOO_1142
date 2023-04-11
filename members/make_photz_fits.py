"""

Created in 2023-04
"""
import numpy as np
import pandas as pd
from pathlib import Path

import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy import units as u
# -------------------------------------------------------------------------------------------- #
location = Path('/Users/emoravec/Documents/Research/merging_clusters/analysis/MOO_1142/members')

def handle_row(row):
    # ID        RA            DEC           z      Mag    Int_Pz    Rad(')   Rank
    # Unpack to separate variables
    id, ra, dec, z, mag, int_pz, rad, rank = row
    # Unpack front and back, pulling dec out specifically
    front, dec, back = row[:3], row[3], row[4:]
    # modify dec in-place
    row[3] = "something else"
    return [id, ra, dec, z, mag, int_pz, rad, rank]
#     rows = [handle_row(line) for line in table_lines]

def read_photz(cat_path):
    '''Reads in the phot-z catalog and makes fits file'''
    with open(cat_path) as file:
        text = file.read()

    # Split the file text into a list of strings (by newline)
    lines = text.splitlines()
    # The second line is the headers
    header_str = lines[0]
    # Everything else is the table
    table_lines = lines[1:]
    # Convert each line to a list of floats by splitting it on whitespace.
    # If the string is 'nan', use np.NaN, otherwise convert fo float
    rows = [[float(s) for s in line.split()] for line in table_lines]
    cat = Table(names=header_str[1:].split(), rows=rows)
    # Convert the id column from float to integer
    cat["Rank"] = cat["Rank"].astype(int)

    return cat

# -------------------------------------------------------------------------------------------- #
# EXECUTE #
wide_cat = read_photz(location / 'specrank.MOO1142.wide_nohst.rank.txt')
print('Len Wide:', len(wide_cat))
wide_tbl = Table(wide_cat)
wide = wide_tbl[wide_tbl['Int_Pz']>0.5]
wide.sort(keys='Int_Pz',reverse=True)
print('Len Wide Int_Pz cut:', len(wide))
#wide = wide[wide["Rad(')"]<2.0]
#print('Len Wide 2am cut:', len(wide))
wide['RA'].name = 'ra'
wide['DEC'].name = 'dec'
wide.write(location/'MOO_1142.photz_wide_IntPz0.5.fits',overwrite=True)

# hst_core_cat = read_photz(location / 'specrank.MOO1142.hst_core.rank.txt')
# hst_core_tbl = Table(hst_core_cat)
# print('Len HST core:', len(hst_core_tbl))
# hst_core= hst_core_tbl[hst_core_tbl['Int_Pz']>0.3]
# hst_core.sort(keys='Int_Pz',reverse=True)
# print('Len HST core Int_Pz cut:', len(hst_core))
# hst_core['RA'].name = 'ra'
# hst_core['DEC'].name = 'dec'
# hst_core.write(location/'MOO_1142.photz_hst_core.fits',overwrite=True)