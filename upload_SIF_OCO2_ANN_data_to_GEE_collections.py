#!/usr/bin/python

import os
import sys
import glob
import pandas as pd
import rasterio as rio
i    try:
        #    try:
        # F    try:
        # F    try:
        # Create an Earth Engine image with properties
        img = ee.Image(f).set({
            "system:time_start": ee.Date(start).millis(),
            "system:time_end": ee.Date(end).millis(),
            "provider": "Oak Ridge National Laboratory (ORNL) Distributed Active Archive Center (DAAC)",
            "URL": "https://daac.ornl.gov/VEGETATION/guides/Global_High_Res_SIF_OCO2.html",
            "nodata_value": nodata_val
        })
        
        # Create and start an upload task
        task = ee.batch.Export.image.toAsset(
            image=img,
            description=f"Upload {basename}",
            assetId=asset.replace("\\", "/"),
            region=img.geometry(),
            maxPixels=1e13,
            pyramidingPolicy={".*": "MEAN"})rth Engine image
        img = ee.Image(f)
        
        # Set the properties
        img = img.set({
            "system:time_start": ee.Date(start).millis(),
            "system:time_end": ee.Date(end).millis(),
            "provider": "Oak Ridge National Laboratory (ORNL) Distributed Active Archive Center (DAAC)",
            "URL": "https://daac.ornl.gov/VEGETATION/guides/Global_High_Res_SIF_OCO2.html",
            "nodata_value": nodata_val
        })
        
        # Create and start an upload task
        task = ee.batch.Export.image.toAsset(
            image=img,
            description=f"Upload {basename}",
            assetId=asset.replace("\\", "/"),
            region=img.geometry(),
            maxPixels=1e13,
            pyramidingPolicy={".*": "MEAN"})
        img = ee.Image.loadGeoTIFF(f)
        
        # Set the properties
        img = img.set({
            "system:time_start": ee.Date(start).millis(),
            "system:time_end": ee.Date(end).millis(),
            "provider": "Oak Ridge National Laboratory (ORNL) Distributed Active Archive Center (DAAC)",
            "URL": "https://daac.ornl.gov/VEGETATION/guides/Global_High_Res_SIF_OCO2.html",
            "nodata_value": nodata_val
        })
        
        # Create and start an upload task
        task = ee.batch.Export.image.toAsset(
            image=img,
            description=f"Upload {basename}",
            assetId=asset.replace("\\", "/"),
            region=img.geometry(),
            maxPixels=1e13,
            pyramidingPolicy={".*": "MEAN"})
        img = ee.Image.loadGeoTIFF(f)
        
        # Set the properties
        img = img.set({
            "system:time_start": ee.Date(start).millis(),
            "system:time_end": ee.Date(end).millis(),
            "provider": "Oak Ridge National Laboratory (ORNL) Distributed Active Archive Center (DAAC)",
            "URL": "https://daac.ornl.gov/VEGETATION/guides/Global_High_Res_SIF_OCO2.html",
            "nodata_value": nodata_val
        })
        
        # Create the export task
        task = ee.batch.Export.image.toAsset 
            image=img,
            description=f"Upload {basename}",
            assetId=asset.replace("\\", "/"),
            region=img.geometry(),
            maxPixels=1e13,
            pyramidingPolicy={".*": "MEAN"}tialize_ee():
    try:
        ee.Initialize(project='phenologymapping')
    except:
        # If initialization fails, try authenticating first
        ee.Authenticate()
        ee.Initialize(project='phenologymapping')

# Initialize Earth Engine
initialize_ee()
# Get dataset from command line arguments
if len(sys.argv) > 1:

    dataset = sys.argv[1].lower()
else:
    dataset = 'sif'  # default value

assert dataset in ['sif'], 'Valid datasets include: "SIF".'

# Adapted from:
    # https://www.tucson.ars.ag.gov/notebooks/uploading_data_2_gee.html

# Define static variables
gc_bucket_dict = {'sif': "seasonality_data/OCO2_SIF_ANN",
                 }
gc_bucket = gc_bucket_dict[dataset]
print(gc_bucket)

img_coll_dict = {'sif': "users/lewiskhu/PhenologyMapping",
                }
img_coll = img_coll_dict[dataset]

# Create the ImageCollection asset in GEE
#cmd = 'earthengine create collection %s' % img_coll
#print("\n Creating ImageCollection using the following command:\n\t%s\n" % cmd)
#os.system(cmd)

# Get the provider and URL strings
str_prov_dict = {'sif': ("(string)provider=Oak Ridge National "
                         "Laboratory (ORNL) Distributed Active "
                         "Archive Center (DAAC)"),
                }
str_prov = str_prov_dict[dataset]

url_dict = {'sif': ('(string)URL=https://daac.ornl.gov/VEGETATION/guides/'
                    'Global_High_Res_SIF_OCO2.html'),
           }
url = url_dict[dataset]

# set the data directory
data_dir_dict = {'sif': r"D:\Lewis\Global phenology maping\output_phen\sif TIF"}
data_dir = data_dir_dict[dataset]

# read in the DataFrame of file start and end dates
df_file_dict = {'sif': r"D:\Lewis\Global phenology maping\output_phen\SIF_OCO2_ANN_upload_metadata.csv",
               }
df_file = df_file_dict[dataset]
df = pd.read_csv(os.path.join(data_dir, df_file))

# set the CRS
crs_dict = {'sif': 'EPSG:4326',
           }
crs = crs_dict[dataset]

# Get file names to extract date and call ingestion command for each file to be added into an asset as image collection
# Example of filenames used here are from a monthly timeseries: rainfall_US_20140101.tif
for f in glob.glob(os.path.join(data_dir, '*.tif')):
    # get the nodata val
    ds = rio.open(f)
    nodata_val = ds.nodata

    # get file basename and its row of data
    basename = os.path.splitext(os.path.basename(f))[0]
    row = df[df['id_no'] == basename]
    print(f'NOW MOVING: {basename}...\n')
    if row.empty:
        print(f"WARNING: No matching metadata for {basename}, skipping.")
        continue

    # get the start and end timestamps
    start = str(row['system:time_start'].values[0])
    end = str(row['system:time_end'].values[0])

    # create asset name
    asset=os.path.join(img_coll, basename)

    # Upload using Earth Engine Python API
    try:
        # Create and start an upload task
        task = ee.batch.Export.image.toAsset(
            image=ee.Image(f),  # Load the image directly
            description=f"Upload {basename}",
            assetId=asset.replace("\\", "/"),
            region=None,  # Auto-detect from image
            scale=None,  # Use native resolution
            maxPixels=1e13,
            pyramidingPolicy={".*": "MEAN"},  # Apply MEAN pyramiding to all bands
            properties={
                "system:time_start": start,
                "system:time_end": end,
                "provider": "Oak Ridge National Laboratory (ORNL) Distributed Active Archive Center (DAAC)",
                "URL": "https://daac.ornl.gov/VEGETATION/guides/Global_High_Res_SIF_OCO2.html",
                "nodata_value": nodata_val
            }
        )
        
        # Start the task
        task.start()
        print(f"Upload task started for {basename}")
        print(f"Task ID: {task.id}")
    except Exception as e:
        print(f"Error uploading {basename}: {str(e)}")
print('Command:\n--------')
print(f"{'='*80}\n")
