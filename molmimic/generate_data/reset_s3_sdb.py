import boto3
import boto.sdb

conn = boto.sdb.connect_to_region("us-east-1")
while len(conn.get_all_domains()) > 0:
    for d in conn.get_all_domains():
       conn.delete_domain(d.name)

buckets = [bucket for bucket in s3.buckets.all() if not bucket.name.startswith("molmimic")]
for bucket in buckets:
    bucket.objects.all().delete()
    for obj in bucket.objects.all():
        obj.delete()
    for version in bucket.object_versions.iterator():
        version.delete()
    bucket.delete()
