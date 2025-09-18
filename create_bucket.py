from google.cloud import storage
import sys

bucket_name = 'phenologymapping-seasonality-data'
location = 'US-CENTRAL1'
project = 'phenologymapping'

print('Creating bucket:', bucket_name, 'in project', project)
client = storage.Client(project=project)
try:
    bucket = client.get_bucket(bucket_name)
    print('Bucket already exists:', bucket_name)
except Exception:
    bucket = storage.Bucket(client, bucket_name)
    bucket.location = location
    bucket = client.create_bucket(bucket)
    print('Bucket created:', bucket.name)
