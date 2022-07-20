from __future__ import print_function
import boto3
import boto.sdb
import botocore

conn = boto.sdb.connect_to_region("us-east-1")
while len(conn.get_all_domains()) > 0:
    for d in conn.get_all_domains():
       conn.delete_domain(d.name)

s3 = boto3.resource('s3', "us-east-1", config=botocore.client.Config(signature_version='s3v4'))
buckets = [bucket for bucket in s3.buckets.all() if bucket.name.endswith("--files")]
for bucket in buckets:
    print("Removing", bucket.name)
    bucket.objects.all().delete()
    for obj in bucket.objects.all():
        obj.delete()
    for version in bucket.object_versions.iterator():
        version.delete()
    bucket.delete()

client = boto3.client('ec2', 'us-east-1')
spots = [s for s in client.describe_spot_instance_requests()["SpotInstanceRequests"] if s["State"] in ["open", "active"]]
for i in range(0, len(spots)+1, 50):
    ids_ = [s["InstanceId"] for s in spots[i:i+50]]
    spotids_ = [s["SpotInstanceRequestId"] for s in spots[i:i+50]]
    client.terminate_instances(InstanceIds=ids_)#[spot["InstanceId"]])
    print("Closing", spotids_)
    client.cancel_spot_instance_requests(SpotInstanceRequestIds=spotids_) #[spot["SpotInstanceRequestId"]])

response = client.describe_instances()
for reservation in response["Reservations"]:
    for instance in reservation["Instances"]:
        if instance["InstanceType"] == "t2.small":
            client.terminate_instances(InstanceIds=[instance["InstanceId"]])
