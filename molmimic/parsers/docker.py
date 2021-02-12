from __future__ import absolute_import

import os
import sys
import subprocess

if "USE_SINGULARITY" in os.environ:
    USE_SINGULARITY = os.environ["USE_SINGULARITY"][0].lower()=="t"
else:
    USE_SINGULARITY = subprocess.check_output(["singularity"]).starts_with(
        "USAGE: singularity")

if not USE_SINGULARITY:
    from toil.lib._docker import apiSingularityCall, singularityKill, \
        singularityStop, containerIsRunning
else:
    from toil.lib._docker import apiDockerCall, dockerKill, dockerStop, \
        containerIsRunning

def dockerCheckOutput(*args, **kwargs):
    raise RuntimeError("dockerCheckOutput() using subprocess.check_output() has been removed, "
                       "please switch to apiDockerCall().")

def dockerCall(*args, **kwargs):
    raise RuntimeError("dockerCall() using subprocess.check_output() has been removed, "
                       "please switch to apiDockerCall().")

def subprocessDockerCall(*args, **kwargs):
    raise RuntimeError("subprocessDockerCall() has been removed, "
                       "please switch to apiDockerCall().")

def apiDockerCall(*args, **kwargs):
    if USE_SINGULARITY:
        return apiSingularityCall(*args, **kwargs)
    else:
        return apiDockerCall(*args, **kwargs)

def dockerKill(container_name, gentleKill=False, timeout=365 * 24 * 60 * 60):
    if USE_SINGULARITY:
        return singularityKill(container_name, gentleKill, timeout)
    else:
        return dockerKill(container_name, gentleKill, timeout)

def dockerStop(container_name):
    if USE_SINGULARITY:
        return singularityStop(container_name)
    else:
        return dockerStop(container_name)

def containerIsRunning(container_name, timeout=365 * 24 * 60 * 60):
    if USE_SINGULARITY:
        return singularityStop(container_name, timeout)
    else:
        return dockerStop(container_name, timeout)

def getContainerName(job):
    """Create a random string including the job name, and return it."""
    return '--'.join([str(job),
                      base64.b64encode(os.urandom(9), b'-_').decode('utf-8')])\
                      .replace("'", '').replace('"', '').replace('_', '')
