import os
import processing

# === INPUT AND OUTPUT FOLDERS ===
input_folder = r"D:\Lewis\Scientific Articles\Urgent\15671260 erthwardphen_asynch-Publication\Global_SIF_OCO2_MODIS_1863\data"
output_folder = r"D:\Lewis\Scientific Articles\Urgent\sif TIF"

# Make sure output folder exists
os.makedirs(output_folder, exist_ok=True)

# List all .cn files
cn_files = [f for f in os.listdir(input_folder) if f.endswith(".nc")]
total_files = len(cn_files)

# Loop through .cn files
for idx, file in enumerate(cn_files, start=1):
    input_path = os.path.join(input_folder, file)
    output_path = os.path.join(output_folder, os.path.splitext(file)[0] + ".tif")
    
    print(f"[{idx}/{total_files}] Converting: {file} → {os.path.basename(output_path)}")
    
    # Run GDAL translate to convert to GeoTIFF
    processing.run("gdal:translate", {
        'INPUT': input_path,
        'TARGET_CRS': None,          # Keep same CRS as input
        'NODATA': None,
        'COPY_SUBDATASETS': False,
        'OPTIONS': '',
        'EXTRA': '',
        'DATA_TYPE': 0,              # Keep original data type
        'OUTPUT': output_path
    })

print("Conversion finished! All .nc files saved as .tif in:", output_folder)
