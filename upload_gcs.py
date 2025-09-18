from google.cloud import storage
from pathlib import Path
import sys

local = Path(r'D:\Lewis\Global phenology maping\output_phen\sif TIF\sif_ann_201409a.tif')
if not local.exists():
    print('Local file not found:', local, file=sys.stderr)
    sys.exit(2)

bucket_name = 'phenologymapping-seasonality-data'
object_path = 'seasonality_data/sif_ann_201409a.tif'

print('Uploading', local, 'to', f'gs://{bucket_name}/{object_path}')
client = storage.Client(project='phenologymapping')
bucket = client.bucket(bucket_name)
blob = bucket.blob(object_path)
blob.upload_from_filename(str(local))
print('Uploaded:', f'gs://{bucket_name}/{object_path}')
