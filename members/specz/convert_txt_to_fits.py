"""
A function to read in text files of speczs and convert to FITS file and pandas DataFrame.
"""
import numpy as np
import pandas as pd
from pathlib import Path

import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy.table import Table,vstack
from astropy import units as u
# -------------------------------------------------------------------------------------------- #
location = Path('/Users/emoravec/Documents/Research/merging_clusters/analysis/MOO_1142/members/specz/')

def read_specz_text_file(table_path,header_start,header_end,column_names_line,table_start,table_end,column_names_for_numeric_value):
    with open(table_path) as file:
        text = file.read()
    # Split the file text into a list of strings (by newline)
    lines = text.splitlines()
    # Define header and column names
    header = lines[header_start:header_end]
    column_names_line = lines[column_names_line][1:] # assuming a hash is at the beginning of the line that has the colun names in it
    column_names = column_names_line.split()
    # table data
    if table_end == ':':
        table_lines = lines[table_start:]
    else:
        table_lines = lines[table_start:table_end]
    rows = [line.split() for line in table_lines]

    # ## Debug reading in table
    # max_cols = len(column_names)
    # print(f"Expected number of columns: {max_cols}")
    # print(f"Column names: {column_names}")

    # for i, row in enumerate(rows):
    #     if len(row) != max_cols:
    #         print(f"Row {i + table_start} has {len(row)} columns: {row}")

    # create table
    tbl = Table(names=column_names, rows=rows)
    # convert dtype string to numeric/float for specified columns
    for col_name in column_names_for_numeric_value:
        tbl[col_name] = tbl[col_name].astype(float)

    return tbl
# -------------------------------------------------------------------------------------------- #
MOO_1142_2025_07_speczs = read_specz_text_file(
    table_path=location / 'M1142+1527.internal.redshifts-30jul25.lst',
    header_start=0,
    header_end=10,
    column_names_line=11,
    table_start=12,
    table_end=':',
    column_names_for_numeric_value=['RA','Dec','Spec-z','Q']
)
# Add Spec-z_err column between Spec-z and Q
specz_index = list(MOO_1142_2025_07_speczs.colnames).index('Spec-z')
MOO_1142_2025_07_speczs.add_column(0.003, name='Spec-z_err', index=specz_index + 1)

# Print the distribution of Q values
print("Distribution of Q values:")
unique_q, counts = np.unique(MOO_1142_2025_07_speczs['Q'], return_counts=True)
for q_val, count in zip(unique_q, counts):
    print(f"Q = {q_val}: {count} sources")
print(f"\nTotal sources: {len(MOO_1142_2025_07_speczs)}")

# Write the table to a FITS file
MOO_1142_2025_07_speczs.write(location / 'MOO_1142+1527.speczs.30jul25.fits', overwrite=True)

# Get just the good quality sources 3/2 and add 1s which are guesses
MOO_1142_2025_07_speczs_good = MOO_1142_2025_07_speczs[np.isin(MOO_1142_2025_07_speczs['Q'], [1, 2, 3])]
MOO_1142_2025_07_speczs_good.write(location / 'MOO_1142+1527.speczs.30jul25.goodQ.fits', overwrite=True)
MOO_1142_2025_07_speczs_good.write(location / 'MOO_1142+1527.speczs.30jul25.goodQ.csv', overwrite=True)