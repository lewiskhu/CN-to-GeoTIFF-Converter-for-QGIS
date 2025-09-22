#!/usr/bin/python

import pandas as pd
import glob
import os
import re

# Directory containing .nc files
data_dir = r"D:\Lewis\Global phenology maping\Global_SIF_OCO2_MODIS_1863\data"
files = glob.glob(os.path.join(data_dir, '*.nc'))
print("Found files:", files)

# Output CSV file
output_file = r"D:\Lewis\Global phenology maping\output_phen\sif_to_tiff\SIF_OCO2_ANN_upload_metadata.csv"

# Columns for CSV
output_cols = ['id_no', 'system:time_start', 'system:time_end']

# Function to get start/end date for each file
def get_file_gee_timestamp(f, start=True):
    date_patt = r'(?<=ann_)\d{6}[ab]'
    date_str = re.search(date_patt, f).group()
    date = '%s-%s-%s' % (date_str[:4],
                         date_str[4:6],
                         str((1 + (15 * (date_str[-1] == 'b')))).zfill(2))
    if not start:
        if date_str[-1] == 'a':
            date = date[:8] + '16'
        else:
            if int(date[5:7]) < 12:
                date = date[:5] + str(int(date[5:7]) + 1).zfill(2) + '-01'
            else:
                date = str(int(date[:4]) + 1) + '-01-01'
    date = date + 'T00:00:00'
    return date

# Build metadata rows
rows = []
for f in files:
    row = {
        'id_no': os.path.splitext(os.path.basename(f))[0],
        'system:time_start': get_file_gee_timestamp(f),
        'system:time_end': get_file_gee_timestamp(f, start=False)
    }
    rows.append(row)

# Convert to DataFrame and save CSV
output_df = pd.DataFrame(rows, columns=output_cols)
print(output_df)
output_df.to_csv(output_file, index=False)
print("Metadata CSV saved to:", output_file)
