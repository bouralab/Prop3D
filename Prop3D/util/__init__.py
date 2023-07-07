import os
import sys
import re
import subprocess
import requests
import shutil
from contextlib import contextmanager

def atof(text, use_int=False):
    converter = int if use_int else float
    try:
        retval = converter(text)
    except ValueError:
        retval = text
    return retval

def natural_keys(text, use_int=False):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    float regex comes from https://stackoverflow.com/a/12643073/190597
    '''
    if isinstance(text, (list,tuple)):
        if len(text)==3 and isinstance(text[0], str) and isinstance(text[0], int) and isinstance(text[2], str):
            return text
        assert 0, "{} must be str".format(text)
    key =  [ atof(c, use_int=use_int) for c in re.split(r'[+-]?([0-9]+(?:[.][0-9]*)?|[.][0-9]+)', str(text)) ]
    assert len(key)==3, "key:{}; text:{}".format(key, text)
    key[0] = ' '.join(key[0].split())
    key[2] = ' '.join(key[2].split())
    if len(key[0]) == "":
        key[0] = " "
    if len(key[2]) == "":
        key[2] = " "
    return tuple(key)

def to_int(s):
    if isinstance(s, int):
        return s
    elif isinstance(s, str):
        return int("".join([d for d in s if d.isdigit()]))
    else:
        raise RuntimeError("Must be an an in or string")

def number_of_lines(path):
    numlines = 0
    with open(path) as f:
        numlines = sum([1 for line in f])
    return numlines

def SubprocessChain(commands, output=subprocess.PIPE):
    if len(commands) > 2:
        prev_proc = subprocess.Popen(
            commands[0],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            env=os.environ)
        for cmd in commands[1:-1]:
            proc = subprocess.Popen(
                cmd,
                stdin=prev_proc.stdout,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                env=os.environ)
            prev_proc = proc
        final_proc = subprocess.Popen(
            commands[-1],
            stdin=prev_proc.stdout,
            stdout=output,
            stderr=subprocess.PIPE,
            env=os.environ)
        return final_proc.communicate()
    elif len(commands) == 2:
        prev_proc = subprocess.Popen(
            commands[0],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            env=os.environ)
        final_proc = subprocess.Popen(
            commands[1],
            stdin=prev_proc.stdout,
            stdout=output,
            stderr=subprocess.PIPE,
            env=os.environ)
        return final_proc.communicate()
    elif len(commands) == 1:
        final_proc = subprocess.Popen(
            commands[0],
            stdout=output,
            stderr=subprocess.PIPE,
            env=os.environ)
    else:
        raise RuntimeError
    return final_proc.communicate()

def safe_remove(files, warn=False):
    if isinstance(files, str):
        files = [files]

    success = []
    for f in files:
        if isinstance(f, str) and os.path.exists(f):
            try:
                os.remove(f)
                success.append(True)
            except (OSError, FileNotFoundError) as e:
                success.append(False)
                pass
        else:
            success.append(False)

    from toil.realtimeLogger import RealtimeLogger
    
    if isinstance(files, str):
        return success[0]

    return success

def safe_call(job, func, input, *args, **kwds):
    try:
        return func(job, input, *args, **kwds)
    except (SystemExit, KeyboardInterrupt):
        raise
    except:
        raise

def download_file(url, local_filename=None):
    if local_filename is None:
        local_filename = url.split('/')[-1]
    with requests.get(url, stream=True) as r:
        with open(local_filename, 'wb') as f:
            r.raw.read = functools.partial(response.raw.read, decode_content=True)
            shutil.copyfileobj(r.raw, f)

    return local_filename

@contextmanager
def silence_stdout():
    new_target = open(os.devnull, "w")
    old_target, sys.stdout = sys.stdout, new_target
    try:
        yield new_target
    finally:
        sys.stdout = old_target

@contextmanager
def silence_stderr():
    new_target = open(os.devnull, "w")
    old_target, sys.stderr = sys.stderr, new_target
    try:
        yield new_target
    finally:
        sys.stderr = old_target

def getcwd():
    """Get current working directory. Uses shell to avoid drive mapping names from python"""
    return subprocess.check_output("pwd", shell=True).decode("utf-8").strip()

def izip_missing(iterA, iterB, **kwds):
    """Iterate through two iterables, while making sure they are in the same
    order. If there are missing values, you can skip the value entirely or
    return only the iterator with the value and a special fill value.
    
    Parameters
    ----------
    iterA : the first iterator
    iterB : the second iterator
    key : function that returns items to compare. Must return strings, ints, or
        an object with the __lt__, __gt__, and __eq__ methods. Optional.
    fillvalue : The value to return if the item is missing. Optional.
    
    Returns
    -------
    A : item from first iterator, or fillValue
    B : item from second iterator, or fillValue
    """
    #Get the comparison functions
    key = kwds.get("key", lambda x: x)
    keyA = kwds.get("keyA", key)
    keyB = kwds.get("keyB", key)

    useMissing = "fillvalue" in kwds
    fillvalue = kwds.get("fillvalue")

    verbose = kwds.get("verbose", False)

    #Start both iterators
    A = next(iterA)
    B = next(iterB)
    try:
        while True:
            if keyA(A) == keyB(B):
                if verbose:
                    print(keyA(A), "==", keyB(B))
                yield A, B
                A = next(iterA)
                B = next(iterB)
            elif keyA(A) < keyB(B):
                if verbose:
                    print(keyA(A), "<", keyB(B))
                if useMissing:
                    yield(A, fillvalue)
                A = next(iterA)
            elif keyA(A) > keyB(B):
                if verbose:
                    print(keyA(A), ">", keyB(B))
                if useMissing:
                    yield fillvalue, B
                B = next(iterB)
            else:
                raise RuntimeError("Invalid comparator")
    except StopIteration:
        pass

def reset_ip():
    import requests
    import boto.ec2

    conn = boto.ec2.connect_to_region("us-east-1")

    try:
        response = requests.get('http://169.254.169.254/latest/meta-data/instance-id')
    except ConnectionError:
        return

    instance_id = response.text

    reservations = conn.get_all_instances(filters={'instance-id' : instance_id})
    instance = reservations[0].instances[0]

    old_address = instance.ip_address


    RealtimeLogger.info("Changing {} IP from {}".format(instance_id, old_address))


    conn.disassociate_address(old_address)

    new_address = conn.allocate_address().public_ip

    RealtimeLogger.info("Changing {} IP from {} to {}".format(instance_id, old_address, new_address))

    conn.associate_address(instance_id, new_address)

    RealtimeLogger.info("Changed {} IP from {} to {}".format(instance_id, old_address, new_address))

def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')
    
def str2boolorlist(v):
    if isinstance(v, bool):
       return v
    elif isinstance(v, str) and v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif isinstance(v, str) and v.lower() in ('no', 'false', 'f', 'n', '0'):
        return None
    if len(v) == 0:
        return True
    elif len(v) == 1 and v[0].lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif len(v) == 1 and v[0].lower() in ('no', 'false', 'f', 'n', '0'):
        return None
    return v