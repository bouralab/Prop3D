from __future__ import absolute_import

try:
    # In Python 3 we have this quote
    from shlex import quote
except ImportError:
    # But in 2.7 we have this deprecated one
    from pipes import quote

import spython.main

import base64
import requests
import logging
import os
import sys

logger = logging.getLogger(__name__)

def apiSingularityCall(job,
                  image,
                  parameters=None,
                  deferParam=None,
                  bind=None,
                  volumes=None,
                  working_dir=None,
                  containerName=None,
                  runscript=None,
                  entrypoint=None,
                  stream=False,
                  detach=False,
                  auto_remove=False,
                  remove=False,
                  user=None,
                  stdout=None,
                  stderr=False,
                  streamfile=None,
                  timeout=365 * 24 * 60 * 60,
                  **kwargs):
    """
    A toil wrapper for the python singularity API.

    Singularity API Docs: https://singularityhub.github.io/singularity-cli/api/
    Singularity API Code: https://github.com/singularityhub/singularity-cli

    This implements singularity's python API within toil so that calls are run as
    jobs, with the intention that failed/orphaned docker jobs be handled
    appropriately.

    Example of using singularityCall in toil to index a FASTA file with SAMtools:
    def toil_job(job):
        working_dir = job.fileStore.getLocalTempDir()
        path = job.fileStore.readGlobalFile(ref_id,
                                          os.path.join(working_dir, 'ref.fasta')
        parameters = ['faidx', path]
        apiSingularityCall(job,
                      image='quay.io/ucgc_cgl/samtools:latest',
                      working_dir=working_dir,
                      parameters=parameters)

    Note that when run with detatch=False, or with detatch=True and stdout=True
    or stderr=True, this is a blocking call. When run with detatch=True and
    without output capture, the container is started and returned without
    waiting for it to finish.

    :param toil.Job.job job: The Job instance for the calling function.
    :param str image: Name of the Singularity image to be used. Most docker
                      images are supported.
                     (e.g. 'quay.io/ucsc_cgl/samtools:latest')
    :param list[str] parameters: A list of string elements.  If there are
                                 multiple elements, these will be joined with
                                 spaces.  This handling of multiple elements
                                 provides backwards compatibility with previous
                                 versions which called docker using
                                 subprocess.check_call().
                                 **If list of lists: list[list[str]], then treat
                                 as successive commands chained with pipe.
    :param str working_dir: The working directory.
    :param int deferParam: Action to take on the container upon job completion.
           FORGO (0) leaves the container untouched and running.
           STOP (1) Sends SIGTERM, then SIGKILL if necessary to the container.
           RM (2) Immediately send SIGKILL to the container. This is the default
           behavior if defer is set to None.
    :param str name: The name/ID of the container.
    :param str runscript: Prepends commands sent to the container.
    :param str entrypoint: See runscript.
    :param bool detach: Run the container in detached mode. (equivalent to '-d')
    :param bool stdout: Return logs from STDOUT when detach=False (default: True).
                        Block and capture stdout to a file when detach=True
                        (default: False). Output capture defaults to output.log,
                        and can be specified with the "streamfile" kwarg.
    :param bool stderr: Return logs from STDERR when detach=False (default: False).
                        Block and capture stderr to a file when detach=True
                        (default: False). Output capture defaults to output.log,
                        and can be specified with the "streamfile" kwarg.
    :param str streamfile: Collect container output to this file if detach=True and
                        stderr and/or stdout are True. Defaults to "output.log".
    :param dict log_config: Specify the logs to return from the container.  See:
                      https://docker-py.readthedocs.io/en/stable/containers.html
    :param bool remove: Remove the container on exit or not.
    :param str user: The container will be run with the privileges of
                     the user specified.  Can be an actual name, such
                     as 'root' or 'lifeisaboutfishtacos', or it can be
                     the uid or gid of the user ('0' is root; '1000' is
                     an example of a less privileged uid or gid), or a
                     complement of the uid:gid (RECOMMENDED), such as
                     '0:0' (root user : root group) or '1000:1000'
                     (some other user : some other user group).
    :param environment: Allows one to set environment variables inside of the
                        container, such as:
    :param int timeout: Use the given timeout in seconds for interactions with
                        the Docker daemon. Note that the underlying docker module is
                        not always able to abort ongoing reads and writes in order
                        to respect the timeout. Defaults to 1 year (i.e. wait
                        essentially indefinitely).
    :param kwargs: Additional keyword arguments supplied to the docker API's
                   run command.  The list is 75 keywords total, for examples
                   and full documentation see:
                   https://docker-py.readthedocs.io/en/stable/containers.html

    :returns: Returns the standard output/standard error text, as requested, when
              detatch=False. Returns the underlying
              docker.models.containers.Container object from the Docker API when
              detatch=True.
    """
    runscript = runscript or entrypoint
    stream = stream or detach

    options = []

    # make certain that files have the correct permissions
    if user is not None:
        options += ["--userns", user]
    else:
        #Add to --userns?
        thisUser = os.getuid()
        thisGroup = os.getgid()
        if user is None:
            user = str(thisUser) + ":" + str(thisGroup)

    if containerName is None:
        containerName = getContainerName(job)

    if working_dir is None:
        working_dir = os.getcwd()
        options += ["--workdir", working_dir]

    if bind is None and volumes is None:
        bind = ["{}:/data".format(working_dir)]
    else:
        if isinstance(bind, str):
            bind = [bind]
        elif isinstance(bind, tuple):
            bind = list(tuple)
        elif volumes is None and not isinstance(bind, list):
            raise RuntimeError("Bind paths must a single string, tuple or list")
        else:
            bind = []

        if isinstance(volumes, dict):
            bind += ["{}:{}".format(host_path, container_path["bind"]) for \
                host_path, container_path in volumes.items()]

    if len(bind) == 0:
        bind = None

    if parameters is None:
        parameters = []

    # If 'parameters' is a list of lists, treat each list as a separate command
    # and chain with pipes.
    if len(parameters) > 0 and isinstance(parameters[0], list):
        if isinstance(runscript, str):
            runscript = [runscript]
        else: #runscript is None:
            runscript = ['/bin/bash', '-c']
        chain_params = \
            [' '.join((quote(arg) for arg in command)) \
             for command in parameters]
        command = ' | '.join(chain_params)
        pipe_prefix = "set -eo pipefail && "
        command = runscript+[pipe_prefix + command]
        logger.debug("Calling singularity with: " + repr(command))

    # If 'parameters' is a normal list, join all elements into a single string
    # element, quoting and escaping each element.
    # Example: ['echo','the Oread'] becomes: ["echo 'the Oread'"]
    # Note that this is still a list, and the docker API prefers this as best
    # practice:
    # http://docker-py.readthedocs.io/en/stable/containers.html
    elif len(parameters) > 0 and isinstance(parameters, list) and runscript is not None:
        if isinstance(runscript, str):
            runscript = [runscript]
        # else: #runscript is None:
        #     pass
            #runscript = ['/bin/bash', '-c']
        command = ' '.join((quote(arg) for arg in parameters))
        logger.debug("Calling singularity with: " + repr(command))

    # If the 'parameters' lists are empty, they are respecified as None, which
    # tells the API to simply create and run the container
    else:
        runscript = None
        command = None

    working_dir = os.path.abspath(working_dir)

    # # Ensure the user has passed a valid value for deferParam
    # assert deferParam in (None, FORGO, STOP, RM), \
    #     'Please provide a valid value for deferParam.'
    #
    # client = spython.main.get_client()
    #
    # if deferParam == STOP:
    #     job.defer(singularityStop, containerName)
    #
    # if deferParam == FORGO:
    #     remove = False
    # elif deferParam == RM:
    #     remove = True
    #     job.defer(singularityKill, containerName)
    # elif remove is True:
    #     job.defer(singularityKill, containerName)

    # if runscript is None:
    #     client_command = client.run
    # else:
    #     client_command = client.execute

    client = spython.main.get_client()

    try:
        if runscript is None:
            out = client.run(image=image,
                       args=parameters,
                       stream=False,
                       bind=bind,
                       **kwargs)
        else:
            out = client.execute(image=image,
                                 command=command,
                                 stream=False,
                                 bind=bind,
                                 **kwargs)
    except Exception as e:
        logger.error("Singularity had non-zero exit. Error: {} \n " +
            "Check your command: {}".format(e, repr(command)))
        raise

    if not stream:
        return out

    if stdout or stderr:
        if streamfile is None:
            streamfile = 'output.log'
        for line in out:
            with open(streamfile, 'w') as f:
                f.write(line)

    # If we didn't capture output, the caller will need iterate on
    # the container to know when it is done. Even if we did capture
    # output, the caller needs the container to get at the exit code.
    return out

def pullSingularityImage(image, pull_folder=None):
    if pull_folder is None:
        pull_folder = os.environ.get("HOME", "")

    base_image = os.path.basename(image)+".simg"

    for dir in (pull_folder, os.getcwd(), os.environ.get("HOME")):
        path = os.path.join(dir, base_image)
        if os.path.isfile(path):
            return path

    client = spython.main.get_client()
    client.pull(image, name=base_image, pull_folder=pull_folder)

    return os.path.join(pull_folder, base_image)


def singularityKill(container_name, gentleKill=False, timeout=365 * 24 * 60 * 60):
    """
    Immediately kills a container.  Equivalent to "docker kill":
    https://docs.docker.com/engine/reference/commandline/kill/
    :param container_name: Name of the container being killed.
    :param client: The docker API client object to call.
    :param int timeout: Use the given timeout in seconds for interactions with
                        the Docker daemon. Note that the underlying docker module is
                        not always able to abort ongoing reads and writes in order
                        to respect the timeout. Defaults to 1 year (i.e. wait
                        essentially indefinitely).
    """
    return


def singularityStop(container_name):
    """
    Gracefully kills a container.  Equivalent to "docker stop":
    https://docs.docker.com/engine/reference/commandline/stop/
    :param container_name: Name of the container being stopped.
    :param client: The docker API client object to call.
    """
    return #dockerKill(container_name, gentleKill=True)


def containerIsRunning(container_name, timeout=365 * 24 * 60 * 60):
    """
    Checks whether the container is running or not.
    :param container_name: Name of the container being checked.
    :returns: True if status is 'running', False if status is anything else,
    and None if the container does not exist.
    :param int timeout: Use the given timeout in seconds for interactions with
                        the Docker daemon. Note that the underlying docker module is
                        not always able to abort ongoing reads and writes in order
                        to respect the timeout. Defaults to 1 year (i.e. wait
                        essentially indefinitely).
    """
    return

def getContainerName(job):
    """Create a random string including the job name, and return it."""
    return '--'.join([str(job),
                      base64.b64encode(os.urandom(9), b'-_').decode('utf-8')])\
                      .replace("'", '').replace('"', '').replace('_', '')
