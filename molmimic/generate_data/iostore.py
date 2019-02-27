"""
IOStore class originated here

https://github.com/BD2KGenomics/hgvm-graph-bakeoff-evaluations/blob/master/scripts/toillib.py

and was then here:

https://github.com/cmarkello/toil-lib/blob/master/src/toil_lib/toillib.py

In a perfect world, this would be deprecated and replaced with Toil's stores.

Actually did this here:

https://github.com/glennhickey/toil-vg/tree/issues/110-fix-iostore

But couldn't get Toil's multipart S3 uploader working on large files.  Also,
the toil jobStore interface is a little less clean for our use.

So for now keep as part of toil-vg where it works.  Could also consider merging
into the upstream toil-lib

https://github.com/BD2KGenomics/toil-lib
"""


import sys, os, os.path, json, collections, logging, logging.handlers
import socketserver, struct, socket, threading, tarfile, shutil
import tempfile
import functools
import random
import time
import dateutil
import traceback
import stat
from toil.realtimeLogger import RealtimeLogger
import subprocess
import datetime
import json

# Need stuff for Amazon s3
try:
    import boto3
    import botocore
    have_s3 = True
    # os.environ["AWS_SECRET_KEY"] = "AKIAJCI3AJKBTBRO6RJA"
    # os.environ["AWS_SECURITY_TOKEN"] = "XB81NeUZsPRHLhMkViFSB4BRAnB0vd1J2U2doXte"
except ImportError:
    have_s3 = False
    pass

# We need some stuff in order to have Azure
try:
    import azure
    # Make sure to get the 0.11 BlobService, in case the new azure storage
    # module is also installed.
    from azure.storage.blob import BlobService
    import toil.jobStores.azureJobStore
    have_azure = True
except ImportError:
    have_azure = False
    pass

def robust_makedirs(directory):
    """
    Make a directory when other nodes may be trying to do the same on a shared
    filesystem.

    """

    if not os.path.exists(directory):
        try:
            # Make it if it doesn't exist
            os.makedirs(directory)
        except OSError:
            # If you can't make it, maybe someone else did?
            pass

    # Make sure it exists and is a directory
    assert(os.path.exists(directory) and os.path.isdir(directory))


def write_global_directory(file_store, path, cleanup=False, tee=None, compress=True):
    """
    Write the given directory into the file store, and return an ID that can be
    used to retrieve it. Writes the files in the directory and subdirectories
    into a tar file in the file store.

    Does not preserve the name or permissions of the given directory (only of
    its contents).

    If cleanup is true, directory will be deleted from the file store when this
    job and its follow-ons finish.

    If tee is passed, a tar.gz of the directory contents will be written to that
    filename. The file thus created must not be modified after this function is
    called.

    """

    write_stream_mode = "w"
    if compress:
        write_stream_mode = "w|gz"

    if tee is not None:
        with open(tee, "w") as file_handle:
            # We have a stream, so start taring into it
            with tarfile.open(fileobj=file_handle, mode=write_stream_mode) as tar:
                # Open it for streaming-only write (no seeking)

                # We can't just add the root directory, since then we wouldn't be
                # able to extract it later with an arbitrary name.

                for file_name in os.listdir(path):
                    # Add each file in the directory to the tar, with a relative
                    # path
                    tar.add(os.path.join(path, file_name), arcname=file_name)

        # Save the file on disk to the file store.
        return file_store.writeGlobalFile(tee)
    else:

        with file_store.writeGlobalFileStream(cleanup=cleanup) as (file_handle,
            file_id):
            # We have a stream, so start taring into it
            # TODO: don't duplicate this code.
            with tarfile.open(fileobj=file_handle, mode=write_stream_mode) as tar:
                # Open it for streaming-only write (no seeking)

                # We can't just add the root directory, since then we wouldn't be
                # able to extract it later with an arbitrary name.

                for file_name in os.listdir(path):
                    # Add each file in the directory to the tar, with a relative
                    # path
                    tar.add(os.path.join(path, file_name), arcname=file_name)

            # Spit back the ID to use to retrieve it
            return file_id

def read_global_directory(file_store, directory_id, path):
    """
    Reads a directory with the given tar file id from the global file store and
    recreates it at the given path.

    The given path, if it exists, must be a directory.

    Do not use to extract untrusted directories, since they could sneakily plant
    files anywhere on the filesystem.

    """

    # Make the path
    robust_makedirs(path)

    with file_store.readGlobalFileStream(directory_id) as file_handle:
        # We need to pull files out of this tar stream

        with tarfile.open(fileobj=file_handle, mode="r|*") as tar:
            # Open it for streaming-only read (no seeking)

            # We need to extract the whole thing into that new directory
            tar.extractall(path)


class IOStore(object):
    """
    A class that lets you get your input files and save your output files
    to/from a local filesystem, Amazon S3, or Microsoft Azure storage
    transparently.

    This is the abstract base class; other classes inherit from this and fill in
    the methods.

    """

    def __init__(self):
        """
        Make a new IOStore
        """

        raise NotImplementedError()

    def connect(self):
        """
        Connect to Resource
        """

        raise NotImplementedError()


    def read_input_file(self, input_path, local_path):
        """
        Read an input file from wherever the input comes from and send it to the
        given path.

        If the file at local_path already exists, it is overwritten.

        If the file at local_path already exists and is a directory, behavior is
        undefined.

        """

        raise NotImplementedError()

    def list_input_directory(self, input_path, recursive=False,
        with_times=False):
        """
        Yields each of the subdirectories and files in the given input path.

        If recursive is false, yields files and directories in the given
        directory. If recursive is true, yields all files contained within the
        current directory, recursively, but does not yield folders.

        If with_times is True, yields (name, modification time) pairs instead of
        just names, with modification times represented as datetime objects in
        the GMT timezone. Modification times may be None on objects that do not
        support them.

        Gives relative file/directory names.

        """

        raise NotImplementedError()

    def write_output_file(self, local_path, output_path):
        """
        Save the given local file to the given output path. No output directory
        needs to exist already.

        If the output path already exists, it is overwritten.

        If the output path already exists and is a directory, behavior is
        undefined.

        """

        raise NotImplementedError()

    def write_output_directory(self, local_direcotry, output_directory=""):
        for dirpath, dirnames, filenames in os.walk(local_direcotry):
            for f in filenames:
                fpath = os.path.join(dirpath, f)
                key = fpath[len(direcotry):]
                self.write_output_file(fpath, output_directory+key)

    def exists(self, path):
        """
        Returns true if the given input or output file exists in the store
        already.

        """

        raise NotImplementedError()

    def get_mtime(self, path):
        """
        Returns the modification time of the given gile if it exists, or None
        otherwise.

        """

        raise NotImplementedError()

    def get_size(self, path):
        """
        Returns the size in bytes of the given file if it exists, or None
        otherwise.

        """

        raise NotImplementedError()

    def get_number_of_items(self, path):
        """
        Return the number of items in path if it exits, or None otherwise
        """

        raise NotImplementedError()

    def remove_file(self, path):
        """
        Removes the object from the store
        """

        raise NotImplementedError()


    @staticmethod
    def absolute(store_string):
        """
        Convert a relative path IOStore string to an absolute path one. Leaves
        strings that aren't FileIOStore specifications alone.

        Since new Toil versions change the working directory of SingleMachine
        batch system jobs, we need to have absolute paths passed into jobs.

        Recommended to be used as an argparse type, so that strings can be
        directly be passed to IOStore.get on the nodes.

        """

        if store_string == "":
            return ""
        if store_string[0] == ".":
            # It's a relative ./ path
            return os.path.abspath(store_string)
        if store_string.startswith("file:"):
            # It's a file:-prefixed thing that may be a relative path
            # Normalize the part after "file:" (which is 5 characters)
            return "file:" + os.path.abspath(store_string[5:])

        return store_string


    @staticmethod
    def get(store_string):
        """
        Get a concrete IOStore created from the given connection string.

        Valid formats are just like for a Toil JobStore, except with container
        names being specified on Azure.

        Formats:

        /absolute/filesystem/path

        ./relative/filesystem/path

        file:filesystem/path

        aws:region:bucket (TODO)

        aws:region:bucket/path/prefix (TODO)

        azure:account:container (instead of a container prefix) (gets keys like
        Toil)

        azure:account:container/path/prefix (trailing slash added automatically)

        """

        # Code adapted from toil's common.py loadJobStore()

        if store_string[0] in "/.":
            # Prepend file: tot he path
            store_string = "file:" + store_string

        try:
            # Break off the first colon-separated piece.
            store_type, store_arguments = store_string.split(":", 1)
        except ValueError:
            # They probably forgot the . or /
            raise RuntimeError("Incorrect IO store specification {}. "
                "Local paths must start with . or /".format(store_string))

        if store_type == "file":
            return FileIOStore(store_arguments)
        elif store_type == "aws":
            # Break out the AWS arguments
            region, bucket_name = store_arguments.split(":", 1)

            if "/" in bucket_name:
                # Split the bucket from the path
                bucket_name, path_prefix = bucket_name.split("/", 1)
            else:
                # No path prefix
                path_prefix = ""

            return S3IOStore(region, bucket_name, path_prefix)
        elif store_type == "azure":
            # Break out the Azure arguments.
            account, container = store_arguments.split(":", 1)

            if "/" in container:
                # Split the container from the path
                container, path_prefix = container.split("/", 1)
            else:
                # No path prefix
                path_prefix = ""

            return AzureIOStore(account, container, path_prefix)
        else:
            raise RuntimeError("Unknown IOStore implementation {}".format(
                store_type))



class FileIOStore(IOStore):
    """
    A class that lets you get input from and send output to filesystem files.

    """

    def __init__(self, path_prefix=""):
        """
        Make a new FileIOStore that just treats everything as local paths,
        relative to the given prefix.

        """

        self.path_prefix = path_prefix

    def read_input_file(self, input_path, local_path):
        """
        Get input from the filesystem.
        """

        RealtimeLogger.debug("Loading {} from FileIOStore in {} to {}".format(
            input_path, self.path_prefix, local_path))

        if os.path.exists(local_path):
            # Try deleting the existing item if it already exists
            try:
                os.unlink(local_path)
            except:
                # Don't fail here, fail complaining about the assertion, which
                # will be more informative.
                pass

        # Make sure the path is clear for copying
        assert(not os.path.exists(local_path))

        # Where is the file actually?
        real_path = os.path.abspath(os.path.join(self.path_prefix, input_path))

        if not os.path.exists(real_path):
            RealtimeLogger.error(
                "Can't find {} from FileIOStore in {}!".format(input_path,
                self.path_prefix))
            raise RuntimeError("File {} missing!".format(real_path))

        # Make a temporary file
        temp_handle, temp_path = tempfile.mkstemp(dir=os.path.dirname(local_path))
        os.close(temp_handle)

        # Copy to the temp file
        shutil.copy2(real_path, temp_path)

        # Rename the temp file to the right place, atomically
        RealtimeLogger.info("rename {} -> {}".format(temp_path, local_path))
        os.rename(temp_path, local_path)

        # Look at the file stats
        file_stats = os.stat(real_path)

        if (file_stats.st_uid == os.getuid() and
            file_stats.st_mode & stat.S_IWUSR):
            # We own this file and can write to it. We don't want the user
            # script messing it up through the symlink.

            try:
                # Clear the user write bit, so the user can't accidentally
                # clobber the file in the actual store through the symlink.
                os.chmod(real_path, file_stats.st_mode ^ stat.S_IWUSR)
            except OSError:
                # If something goes wrong here (like us not having permission to
                # change permissions), ignore it.
                pass

    def list_input_directory(self, input_path, recursive=False,
        with_times=False):
        """
        Loop over directories on the filesystem.
        """

        RealtimeLogger.info("Enumerating {} from "
            "FileIOStore in {}".format(input_path, self.path_prefix))

        if not os.path.exists(os.path.join(self.path_prefix, input_path)):
            # Nothing to list over
            return

        if not os.path.isdir(os.path.join(self.path_prefix, input_path)):
            # Can't list a file, only a directory.
            return

        for item in os.listdir(os.path.join(self.path_prefix, input_path)):
            if(recursive and os.path.isdir(os.path.join(self.path_prefix,
                input_path, item))):
                # We're recursing and this is a directory.
                # Recurse on this.
                for subitem in self.list_input_directory(
                    os.path.join(input_path, item), recursive):

                    # Make relative paths include this directory name and yield
                    # them
                    name_to_yield = os.path.join(item, subitem)

                    if with_times:
                        # What is the mtime in seconds since epoch?
                        mtime_epoch_seconds = os.path.getmtime(os.path.join(
                            input_path, item, subitem))
                        # Convert it to datetime

                        yield name_to_yield, mtime_epoch_seconds
                    else:
                        yield name_to_yield
            else:
                # This isn't a directory or we aren't being recursive
                # Just report this individual item.

                if with_times:
                    # What is the mtime in seconds since epoch?
                    mtime_epoch_seconds = os.path.getmtime(os.path.join(
                        input_path, item))

                    yield item, mtime_epoch_seconds
                else:
                    yield item

    def write_output_file(self, local_path, output_path):
        """
        Write output to the filesystem
        """

        RealtimeLogger.debug("Saving {} to FileIOStore in {}".format(
            output_path, self.path_prefix))

        # What's the real output path to write to?
        real_output_path = os.path.join(self.path_prefix, output_path)

        # What directory should this go in?
        parent_dir = os.path.split(real_output_path)[0]

        if parent_dir != "":
            # Make sure the directory it goes in exists.
            robust_makedirs(parent_dir)

        # Make a temporary file
        temp_handle, temp_path = tempfile.mkstemp(dir=self.path_prefix)
        os.close(temp_handle)

        # Copy to the temp file
        shutil.copy2(local_path, temp_path)

        if os.path.exists(real_output_path):
            # At least try to get existing files out of the way first.
            os.unlink(real_output_path)

        # Rename the temp file to the right place, atomically
        os.rename(temp_path, real_output_path)

    def exists(self, path):
        """
        Returns true if the given input or output file exists in the file system
        already.

        """

        return os.path.exists(os.path.join(self.path_prefix, path))

    def get_mtime(self, path):
        """
        Returns the modification time of the given file if it exists, or None
        otherwise.

        """

        if not self.exists(path):
            return None

        # What is the mtime in seconds since epoch?
        mtime_epoch_seconds = os.path.getmtime(os.path.join(self.path_prefix,
            path))
        # Convert it to datetime
        mtime_datetime = datetime.datetime.utcfromtimestamp(
            mtime_epoch_seconds).replace(tzinfo=dateutil.tz.tzutc())

        # Return the modification time, timezoned, in UTC
        return mtime_datetime

    def get_size(self, path):
        """
        Returns the size in bytes of the given file if it exists, or None
        otherwise.

        """

        if not self.exists(path):
            return None

        # Return the size in bytes of the backing file
        return os.stat(os.path.join(self.path_prefix, path)).st_size

    def get_number_of_items(self, path=None):
        """
        Return the number of items in path if it is exist, or None otherwise
        """
        if path is not None and self.exists(path):
            return None
        elif path is None:
            path = self.path_prefix

        return sum(1 for name in os.listdir(path) if os.path.isfile(name))

    def remove_file(self, path):
        try:
            os.remove(path)
        except OSError:
            pass

class BackoffError(RuntimeError):
    """
    Represents an error from running out of retries during exponential back-off.
    """

def backoff_times(retries, base_delay):
    """
    A generator that yields times for random exponential back-off. You have to
    do the exception handling and sleeping yourself. Stops when the retries run
    out.

    """

    # Don't wait at all before the first try
    yield 0

    # What retry are we on?
    try_number = 1

    # Make a delay that increases
    delay = float(base_delay) * 2

    while try_number <= retries:
        # Wait a random amount between 0 and 2^try_number * base_delay
        yield random.uniform(base_delay, delay)
        delay *= 2
        try_number += 1

    # If we get here, we're stopping iteration without succeeding. The caller
    # will probably raise an error.

def backoff(original_function, retries=6, base_delay=10):
    """
    We define a decorator that does randomized exponential back-off up to a
    certain number of retries. Raises BackoffError if the operation doesn't
    succeed after backing off for the specified number of retries (which may be
    float("inf")).

    Unfortunately doesn't really work on generators.
    """

    # Make a new version of the function
    @functools.wraps(original_function)
    def new_function(*args, **kwargs):
        last_error = None
        # Call backoff times, overriding parameters with stuff from kwargs
        for delay in backoff_times(retries=kwargs.get("retries", retries),
            base_delay=kwargs.get("base_delay", base_delay)):
            # Keep looping until it works or our iterator raises a
            # BackoffError
            if delay > 0:
                # We have to wait before trying again
                RealtimeLogger.error("Retry after {} seconds".format(
                    delay))
                time.sleep(delay)
            try:
                return original_function(*args, **kwargs)
            except Exception as last_error:
                # Report the formatted underlying exception with traceback
                RealtimeLogger.error("{} failed due to: {}".format(
                    original_function.__name__,
                    "".join(traceback.format_exception(*sys.exc_info()))))


        # If we get here, the function we're calling never ran through before we
        # ran out of backoff times. Give an error.
        if last_error is None:
            raise BackoffError("Ran out of retries calling {}".format(
                original_function.__name__))
        else:
            raise last_error

    return new_function

class S3IOStore(IOStore):
    """
    A class that lets you get input from and send output to AWS S3 Storage.

    """

    def __init__(self, region, bucket_name, name_prefix=""):
        """
        Make a new S3IOStore that reads from and writes to the given
        container in the given account, adding the given prefix to keys. All
        paths will be interpreted as keys or key prefixes.

        """

        # Make sure s3 libraries actually loaded
        assert(have_s3)

        self.region = region
        self.bucket_name = bucket_name
        self.name_prefix = name_prefix
        self.s3 = None

    def connect(self):
        self.__connect()

    def __connect(self):
        """
        Make sure we have an S3 Bucket connection, and set one up if we don't.
        Creates the S3 bucket if it doesn't exist.
        """

        if self.s3 is None:
            RealtimeLogger.debug("Connecting to bucket {} in region {}".format(
                self.bucket_name, self.region))
            print("Connecting to bucket {} in region {}".format(
                self.bucket_name, self.region))

            # Connect to the s3 bucket service where we keep everything
            self.s3 = boto3.client('s3', self.region, config=
                botocore.client.Config(signature_version='s3v4', retries={"max_attempts":20}))
            self.s3r = boto3.resource('s3', self.region, config=
                botocore.client.Config(signature_version='s3v4', retries={"max_attempts":20}))
            try:
                self.s3.head_bucket(Bucket=self.bucket_name)
            except:
                if self.region == 'us-east-1':
                    self.s3.create_bucket(
                        Bucket=self.bucket_name,
                    )
                else:
                    self.s3.create_bucket(
                        Bucket=self.bucket_name,
                        CreateBucketConfiguration={'LocationConstraint': self.region},
                    )

    #@backoff
    def read_input_file(self, input_path, local_path):
        """
        Get input from S3.
        """

        self.__connect()

        RealtimeLogger.debug("Loading {} from S3IOStore".format(
            input_path))

        # Download the file contents.
        self.s3.download_file(self.bucket_name, os.path.join(self.name_prefix, input_path), local_path)

    @backoff
    def list_input_directory(self, input_path=None, recursive=False,
        with_times=False):
        """
        Yields each of the subdirectories and files in the given input path.

        If recursive is false, yields files and directories in the given
        directory. If recursive is true, yields all files contained within the
        current directory, recursively, but does not yield folders.

        If with_times is True, yields (name, modification time) pairs instead of
        just names, with modification times represented as datetime objects in
        the GMT timezone. Modification times may be None on objects that do not
        support them.

        Gives relative file/directory names.

        """
        if with_times:
            get_output = lambda f: (f.key, f.last_modified)
        else:
            get_output = lambda f: f.key

        bucket = self.s3r.Bucket(self.bucket_name).objects.all() if input_path is None \
            else self.s3r.Bucket(self.bucket_name).objects.filter(Prefix=input_path)

        for obj in bucket:
            yield get_output(obj)

    def write_output_file(self, local_path, output_path):
        """
        Write output to S3.
        """

        self.__connect()

        RealtimeLogger.debug("Saving {} to S3IOStore".format(
            output_path))

        # Download the file contents.
        #self.s3.upload_file(local_path, self.bucket_name, os.path.join(self.name_prefix, output_path))
        self.s3r.Bucket(self.bucket_name).upload_file(local_path, os.path.join(self.name_prefix, output_path))

    @backoff
    def exists(self, path):
        """
        Returns true if the given input or output file exists in the store
        already.

        """

        self.__connect()

        try:
            self.s3.get_object(Bucket=self.bucket_name, Key=path)
            return True
        except (SystemExit, KeyboardInterrupt):
            raise
        except:
            return False

    def get_mtime(self, path):
        """
        Returns the modification time of the given file if it exists, or None
        otherwise.

        """

        raise NotImplementedError()

    def get_size(self, path):
        """
        Returns the size in bytes of the given file if it exists, or None
        otherwise.

        """

        raise NotImplementedError()

    @backoff
    def get_number_of_items(self, path=None):
        """
        Return the number of items in path if it exits, or None otherwise
        """

        return sum(1 for _ in self.list_input_directory(path))

    def remove_file(self, path):
        self.s3r.Object(self.bucket_name, path).delete()


class AzureIOStore(IOStore):
    """
    A class that lets you get input from and send output to Azure Storage.

    """

    def __init__(self, account_name, container_name, name_prefix=""):
        """
        Make a new AzureIOStore that reads from and writes to the given
        container in the given account, adding the given prefix to keys. All
        paths will be interpreted as keys or key prefixes.

        If the name prefix does not end with a trailing slash, and is not empty,
        one will be added automatically.

        Account keys are retrieved from the AZURE_ACCOUNT_KEY environment
        variable or from the ~/.toilAzureCredentials file, as in Toil itself.

        """

        # Make sure azure libraries actually loaded
        assert(have_azure)

        self.account_name = account_name
        self.container_name = container_name
        self.name_prefix = name_prefix

        if self.name_prefix != "" and not self.name_prefix.endswith("/"):
            # Make sure it has the trailing slash required.
            self.name_prefix += "/"

        # Sneak into Toil and use the same keys it uses
        self.account_key = toil.jobStores.azureJobStore._fetchAzureAccountKey(
            self.account_name)

        # This will hold out Azure blob store connection
        self.connection = None

    def __getstate__(self):
        """
        Return the state to use for pickling. We don't want to try and pickle
        an open Azure connection.
        """

        return (self.account_name, self.account_key, self.container_name,
            self.name_prefix)

    def __setstate__(self, state):
        """
        Set up after unpickling.
        """

        self.account_name = state[0]
        self.account_key = state[1]
        self.container_name = state[2]
        self.name_prefix = state[3]

        self.connection = None

    def __connect(self):
        """
        Make sure we have an Azure connection, and set one up if we don't.
        """

        if self.connection is None:
            RealtimeLogger.debug("Connecting to account {}, using "
                "container {} and prefix {}".format(self.account_name,
                self.container_name, self.name_prefix))

            # Connect to the blob service where we keep everything
            self.connection = BlobService(
                account_name=self.account_name, account_key=self.account_key)


    @backoff
    def read_input_file(self, input_path, local_path):
        """
        Get input from Azure.
        """

        self.__connect()


        RealtimeLogger.debug("Loading {} from AzureIOStore".format(
            input_path))

        # Download the blob. This is known to be synchronous, although it can
        # call a callback during the process.
        self.connection.get_blob_to_path(self.container_name,
            self.name_prefix + input_path, local_path)

    def list_input_directory(self, input_path, recursive=False,
        with_times=False):
        """
        Loop over fake /-delimited directories on Azure. The prefix may or may
        not not have a trailing slash; if not, one will be added automatically.

        Returns the names of files and fake directories in the given input fake
        directory, non-recursively.

        If with_times is specified, will yield (name, time) pairs including
        modification times as datetime objects. Times on directories are None.

        """

        self.__connect()

        RealtimeLogger.info("Enumerating {} from AzureIOStore".format(
            input_path))

        # Work out what the directory name to list is
        fake_directory = self.name_prefix + input_path

        if fake_directory != "" and not fake_directory.endswith("/"):
            # We have a nonempty prefix, and we need to end it with a slash
            fake_directory += "/"

        # This will hold the marker that we need to send back to get the next
        # page, if there is one. See <http://stackoverflow.com/a/24303682>
        marker = None

        # This holds the subdirectories we found; we yield each exactly once if
        # we aren't recursing.
        subdirectories = set()

        while True:

            # Get the results from Azure. We don't use delimiter since Azure
            # doesn't seem to provide the placeholder entries it's supposed to.
            result = self.connection.list_blobs(self.container_name,
                prefix=fake_directory, marker=marker)

            RealtimeLogger.info("Found {} files".format(len(result)))

            for blob in result:
                # Yield each result's blob name, but directory names only once

                # Drop the common prefix
                relative_path = blob.name[len(fake_directory):]

                if (not recursive) and "/" in relative_path:
                    # We found a file in a subdirectory, and we aren't supposed
                    # to be recursing.
                    subdirectory, _ = relative_path.split("/", 1)

                    if subdirectory not in subdirectories:
                        # It's a new subdirectory. Yield and remember it
                        subdirectories.add(subdirectory)

                        if with_times:
                            yield subdirectory, None
                        else:
                            yield subdirectory
                else:
                    # We found an actual file
                    if with_times:
                        mtime = blob.properties.last_modified

                        if isinstance(mtime, datetime.datetime):
                            # Make sure we're getting proper localized datetimes
                            # from the new Azure Storage API.
                            assert(mtime.tzinfo is not None and
                                mtime.tzinfo.utcoffset(mtime) is not None)
                        else:
                            # Convert mtime from a string as in the old API.
                            mtime = dateutil.parser.parse(mtime).replace(
                                tzinfo=dateutil.tz.tzutc())

                        yield relative_path, mtime

                    else:
                        yield relative_path

            # Save the marker
            marker = result.next_marker

            if not marker:
                break

    @backoff
    def write_output_file(self, local_path, output_path):
        """
        Write output to Azure. Will create the container if necessary.
        """

        self.__connect()

        RealtimeLogger.debug("Saving {} to AzureIOStore".format(
            output_path))

        try:
            # Make the container
            self.connection.create_container(self.container_name)
        except azure.WindowsAzureConflictError:
            # The container probably already exists
            pass

        # Upload the blob (synchronously)
        # TODO: catch no container error here, make the container, and retry
        self.connection.put_block_blob_from_path(self.container_name,
            self.name_prefix + output_path, local_path)

    @backoff
    def exists(self, path):
        """
        Returns true if the given input or output file exists in Azure already.

        """

        self.__connect()

        marker = None

        while True:

            try:
                # Make the container
                self.connection.create_container(self.container_name)
            except azure.WindowsAzureConflictError:
                # The container probably already exists
                pass

            # Get the results from Azure.
            result = self.connection.list_blobs(self.container_name,
                prefix=self.name_prefix + path, marker=marker)

            for blob in result:
                # Look at each blob

                if blob.name == self.name_prefix + path:
                    # Found it
                    return True

            # Save the marker
            marker = result.next_marker

            if not marker:
                break

        return False


    @backoff
    def get_mtime(self, path):
        """
        Returns the modification time of the given blob if it exists, or None
        otherwise.

        """

        self.__connect()

        marker = None

        while True:

            # Get the results from Azure.
            result = self.connection.list_blobs(self.container_name,
                prefix=self.name_prefix + path, marker=marker)

            for blob in result:
                # Look at each blob

                if blob.name == self.name_prefix + path:
                    # Found it
                    mtime = blob.properties.last_modified

                    if isinstance(mtime, datetime.datetime):
                        # Make sure we're getting proper localized datetimes
                        # from the new Azure Storage API.
                        assert(mtime.tzinfo is not None and
                            mtime.tzinfo.utcoffset(mtime) is not None)
                    else:
                        # Convert mtime from a string as in the old API.
                        mtime = dateutil.parser.parse(mtime).replace(
                            tzinfo=dateutil.tz.tzutc())

                    return mtime

            # Save the marker
            marker = result.next_marker

            if not marker:
                break

        return None

    @backoff
    def get_size(self, path):
        """
        Returns the size in bytes of the given blob if it exists, or None
        otherwise.

        """

        self.__connect()

        marker = None

        while True:

            # Get the results from Azure.
            result = self.connection.list_blobs(self.container_name,
                prefix=self.name_prefix + path, marker=marker)

            for blob in result:
                # Look at each blob

                if blob.name == self.name_prefix + path:
                    # Found it
                    size = blob.properties.content_length

                    return size

            # Save the marker
            marker = result.next_marker

            if not marker:
                break

        return None

    def get_number_of_items(self, path=None):
        """
        Return the number of items in path if it exits, or None otherwise
        """

        return sum(1 for _ in self.list_input_directory(path))
