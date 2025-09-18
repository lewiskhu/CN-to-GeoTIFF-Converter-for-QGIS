#!/usr/bin/python

import os
import glob
import pandas as pd
from osgeo import gdal  # use GDAL instead of rasterio

# -----------------------------
# Dataset configuration
# -----------------------------
dataset = 'sif'  # Only 'sif' supported in this script

# GEE Asset paths
gc_bucket_dict = {'sif': "seasonality_data/OCO2_SIF_ANN"}  # Google Cloud Storage bucket folder
gc_bucket = gc_bucket_dict[dataset]

img_coll_dict = {'sif': "users/drewhart/seasonality_data/OCO2_SIF_ANN"}  # GEE Image Collection asset
img_coll = img_coll_dict[dataset]

str_prov_dict = {'sif': ("(string)provider=Oak Ridge National Laboratory "
                         "(ORNL) Distributed Active Archive Center (DAAC)")}
str_prov = str_prov_dict[dataset]

url_dict = {'sif': ('(string)URL=https://daac.ornl.gov/VEGETATION/guides/'
                    'Global_High_Res_SIF_OCO2.html')}
url = url_dict[dataset]

crs_dict = {'sif': 'EPSG:4326'}
crs = crs_dict[dataset]

# -----------------------------
# Local paths (Windows)
# -----------------------------
tif_dir = r"D:\Lewis\Global phenology maping\output_phen\sif"  # Folder with .tif files
csv_file = r"D:\Lewis\Global phenology maping\output_phen\SIF_OCO2_ANN_upload_metadata.csv"  # CSV with start/end dates

# Read CSV
df = pd.read_csv(csv_file)

# -----------------------------
# Process each .tif file
# -----------------------------
tif_files = glob.glob(os.path.join(tif_dir, '*.tif'))
if not tif_files:
    print("No .tif files found in:", tif_dir)
    exit(1)

for f in tif_files:
    # Open raster with GDAL to get nodata value
    ds = gdal.Open(f)
    nodata_val = ds.GetRasterBand(1).GetNoDataValue()

    # File basename
    basename = os.path.splitext(os.path.basename(f))[0]
    row = df[df['id_no'] == basename]

    if row.empty:
        print(f"Metadata not found for {basename}, skipping...")
        continue

    start = str(row['system:time_start'].values[0])
    end = str(row['system:time_end'].values[0])

    # GEE asset name
    asset = os.path.join(img_coll, basename).replace("\\", "/")  # forward slashes for GEE

    # GCS path
    gcs_path = f"gs://{gc_bucket}/{os.path.basename(f)}"

    # Build Earth Engine upload command
    cmd = (f'earthengine upload image --asset_id="{asset}" '
           f'--pyramiding_policy=mean --time_start="{start}" --time_end="{end}" '
           f"--property='{str_prov}' --property='{url}' "
           f'--nodata_value={nodata_val} '
           f'--crs={crs} {gcs_path}')

    print(f"Uploading {basename}...")
    print(cmd)

    # Execute
    os.system(cmd)
    print("="*80)

print("All .tif files processed.")
