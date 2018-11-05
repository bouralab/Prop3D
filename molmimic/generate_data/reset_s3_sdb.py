import boto3
import boto.sdb
import botocore

conn = boto.sdb.connect_to_region("us-east-1")
while len(conn.get_all_domains()) > 0:
    for d in conn.get_all_domains():
       conn.delete_domain(d.name)

s3 = boto3.resource('s3', "us-east-1", config=botocore.client.Config(signature_version='s3v4'))
buckets = [bucket for bucket in s3.buckets.all() if not bucket.name.startswith("molmimic")]
for bucket in buckets:
    bucket.objects.all().delete()
    for obj in bucket.objects.all():
        obj.delete()
    for version in bucket.object_versions.iterator():
        version.delete()
    bucket.delete()

client = boto3.client('ec2')
response = client.describe_instances()
for reservation in response["Reservations"]:
    for instance in reservation["Instances"]:
        if instance["InstanceType"] == "t2.small":
            client.terminate_instances(InstanceIds=[instance["InstanceId"]])

spots = client.describe_spot_instance_requests()
for spot in spots["SpotInstanceRequests"]:
    if spot["State"] == "open":
        print "Closing", spot["SpotInstanceRequestId"]
        client.cancel_spot_instance_requests(SpotInstanceRequestIds=[spot["SpotInstanceRequestId"]])
