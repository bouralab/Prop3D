import os
import sys
import subprocess
import tempfile
import shutil

from toil.realtimeLogger import RealtimeLogger
from molmimic.util import silence_stdout, silence_stderr

CONTAINER_PATH = os.environ.get("CONTAINER_PATH", os.environ["HOME"])

FORCE_LOCAL = os.environ.get("FORCE_LOCAL", "false")[0].lower()=="t"
USE_SINGULARITY = os.environ.get("USE_SINGULARITY", "true")[0].lower()=="t"
USE_DOCKER = os.environ.get("USE_DOCKER", "false")[0].lower()=="t"

if FORCE_LOCAL:
    USE_SINGULARITY = False
    USE_DOCKER = False
elif USE_SINGULARITY:
    USE_DOCKER = False
    try:
        rc = subprocess.check_output(["singularity"]).decode("utf-8")
        USE_SINGULARITY = rc.startswith("USAGE: singularity")
    except FileNotFoundError:
        USE_SINGULARITY = False
    FORCE_LOCAL = not USE_SINGULARITY
elif USE_DOCKER:
    USE_SINGULARITY = False
    try:
        rc = subprocess.check_output(["docker"]).decode("utf-8")
        USE_DOCKER = rc.strip().startswith("Usage:	docker")
    except FileNotFoundError:
        USE_DOCKER = False
    FORCE_LOCAL = not USE_DOCKER

if USE_SINGULARITY:
    from molmimic.parsers.singularity import apiSingularityCall as containerCall
    from molmimic.parsers.singularity import singularityKill as containerKill
    from molmimic.parsers.singularity import singularityStop as containerStop
    from molmimic.parsers.singularity import containerIsRunning
    from molmimic.parsers.singularity import pullSingularityImage as pullContainer
elif USE_DOCKER:
    from toil.lib.docker import apiDockerCall as containerCall
    from toil.lib.docker import dockerKill as containerKill
    from toil.lib.docker import dockerStop as containerStop
    from toil.lib.docker import containerIsRunning

if USE_DOCKER or FORCE_LOCAL:
    #Docker automatically caches image, no need to save
    def pullContainer(image, pull_folder=""): return image

class Container(object):
    """

    Arguments: A list of 2- or 3-tuples countaining the parameter name, rule,
        and optioanlly a string formatting rule.

    Rules: A function or custom method name as string. Predefined rules:
        'str': No action, stay as string
        'path:in': Define a path to be used as input. Paths can be modified for
            use in contianer.
        'path:out': Define a path to be used as output. Paths created in the
            container can be copied to the local machine
        'path:in:stdin': Define a path to be used as input. For contianerized it
            will use the 'in' and LOCAL commands will use stdin.

        Custom methods must return the updated value. If None, value will be
        removed from parameter list.
    """

    IMAGE = None
    LOCAL = None
    PARAMETERS = []
    RULES = {}
    ENTRYPOINT = None
    RETURN_FILES = False

    def __init__(self, job=None, return_files=False, force_local=False, fallback_local=False, work_dir=None):
        assert (self.IMAGE, self.LOCAL).count(None) <= 1, "Must define container or local path"

        if self.LOCAL is not None:
            assert isinstance(self.LOCAL, list)

        self.work_dir = work_dir if work_dir is None else os.getcwd()
        self.force_local = (FORCE_LOCAL or force_local or self.IMAGE is None) and self.LOCAL is not None
        self.fallback_local = fallback_local or self.LOCAL is not None
        self.return_files = return_files or self.RETURN_FILES
        self.job = job
        self.rules = {
            str: lambda k, v: str(v),
            "str": lambda k, v: str(v),
            None: lambda k, v: v,
            "None": lambda k, v: v,
            "path:in": self.format_in_path,
            "path:out": self.format_out_path,
            "path:in:stdin": "format_in_std_in",
            "store_true": "store_true"}
        self.rules.update(self.RULES)

        self.parameters = []
        self.param_names = []
        self.param_funcs = {}
        self.params_to_update = {}
        self.parameter_formatters = {}

        self.optional_parameters = {}
        self.optional_param_names = []
        self.optional_param_funcs = {}
        self.optional_params_to_update = {}
        self.optional_parameter_formatters = {}
        self.optional_parameters_defaults = {}

        self.number_of_parameters = 0
        self.number_of_optional_parameters = 0

        self.change_paths = {}
        self.files_to_remove = []

        self.is_local = False
        self.stdin = None

        for i, p in enumerate(self.PARAMETERS):
            if isinstance(p, str):
                self.parameters.append(p)
            elif isinstance(p, (list, tuple)):
                if len(p) == 1 and isinstance(p[0], str):
                    #ignore, treat as string
                    self.parameters.append(p)
                elif len(p) in [2,3]:
                    key, rule = p[0], p[1]

                    optional = False
                    default = None
                    if key.startswith(":"):
                        #Optional value
                        key_parts = key[1:].split(":")
                        key = key_parts[0]
                        default = key_parts[1] if len(key_parts) > 1 else None

                        if rule == "store_true":
                            default = None

                        optional = True

                    if len(p) == 3:
                        if isinstance(p[2], str) and len(p[2])>0:
                            formatter = p[2]
                        else:
                            #Formatter is probably None, which means the parameter
                            #is saved for other methods to use, but is not passed
                            #to the program
                            formatter = None
                    else:
                        formatter = "{}"

                    self.parameter_formatters[key] = formatter

                    if rule in self.rules:
                        if callable(self.rules[rule]):
                            func = self.rules[rule]
                        elif hasattr(self, self.rules[rule]):
                            func = getattr(self, self.rules[rule])
                        elif hasattr(self, rule):
                            func = getattr(self, rule)
                        else:
                            assert 0, rule

                        if not optional:
                            if default is not None:
                                self.parameters.append(formatter.format(default))
                            else:
                                self.parameters.append(None)
                            self.param_names.append(key)
                            self.param_funcs[i] = func
                            self.params_to_update[key] = i
                            self.number_of_parameters += 1
                        else:
                            if default is not None:
                                self.parameters.append(formatter.format(default))
                            else:
                                self.parameters.append(None)
                            self.optional_param_names.append(key)
                            self.optional_param_funcs[i] = func
                            self.optional_params_to_update[key] = i
                            self.number_of_optional_parameters += 1
                            self.optional_parameters_defaults[key] = default
                    else:
                        raise RuntimeError("invalid parameters: {}".format(p))

    def __call__(self, *args, **kwds):
        if self.force_local:
            return self.local(*args, **kwds)

        assert self.job is not None

        parameters = self.format_parameters(args, kwds)
        RealtimeLogger.info("Run {} with params {}".format(self.IMAGE, parameters))

        image = pullContainer(self.IMAGE, pull_folder=CONTAINER_PATH)

        try:
            with silence_stdout(), silence_stderr():
                out = containerCall(
                    self.job,
                    image=image,
                    entrypoint=self.ENTRYPOINT,
                    working_dir="/data",
                    volumes={self.work_dir:{"bind":"/data", "mode":"rw"}},
                    parameters=parameters,
                    return_result=True)
        except (SystemExit, KeyboardInterrupt):
            raise
        except:
            if self.fallback_local:
                import traceback as tb
                RealtimeLogger.error(tb.format_exc())
                return self.local(*args, **kwds)
            raise

        message = "".join(out["message"])

        if out["return_code"]:
            self.clean()
            raise RuntimeError(message)

        try:
            out_files = self.check_output()
        except AssertionError:
            print("captured output")
            self.stdout = out
            self.clean()
            raise

        if self.return_files:
            self.stdout = message
            self.clean()
            return out_files

        return message

    def local(self, *args, **kwds):
        self.is_local = True
        parameters = self.format_parameters(args, kwds)

        env = self.set_local_env()

        try:
            if self.stdin is not None:
                with open(self.stdin) as f:
                    out = subprocess.check_output(self.LOCAL+parameters,
                        env=env, stdin=f)
            else:
                out = subprocess.check_output(self.LOCAL+parameters, env=env)
        except (SystemExit, KeyboardInterrupt):
            raise
        except:
            raise
        self.is_local = False
        self.stdin = None
        self.clean()

    def set_local_env(self):
        return None

    def format_parameters(self, args, kwds):
        if len(args)+len(kwds) == self.number_of_parameters:
            pass

        # if len(args) == self.number_of_parameters and len(kwds) == self.number_of_optional_parameters":
        #     #Correct number of args, all default options
        #     kwds = zip(self.param_names, args).update(kwds)
        if len(args) == self.number_of_parameters:
            #Correct number of args, maybe optionsals optionals
            kwds.update(dict(zip(self.param_names, args)))
        elif len(args) == 0 and len(kwds) >= self.number_of_parameters:
            #Correct number of args, optionals may be present
            pass
        else:
            #Less args than kwds, make sure missing args are in kwds
            params = dict(zip(self.param_names[:len(args)], args))
            need_params = self.param_names[len(args):]
            if len(need_params) == 0 or all(k in kwds for k in need_params):
                kwds = params.update(kwds)
            else:
                raise RuntimeError

        parameters = self.parameters[:]
        for k, v in kwds.items():
            try:
                idx = self.params_to_update[k]
                param_func = self.param_funcs[idx]
            except KeyError as e:
                try:
                    idx = self.optional_params_to_update[k]
                    param_func = self.optional_param_funcs[idx]
                except KeyError as e:
                    #Invalid parameter, skip
                    RealtimeLogger.info("ignoring parameter {}".format(k))
                    continue

            val = param_func(k, v)

            if val is None:
                RealtimeLogger.info("ignoring parameter {} since val is None".format(k))
                continue

            try:
                formatter = self.parameter_formatters[k]
                if formatter is not None:
                    parameters[idx] = formatter.format(val)
                else:
                    parameters[idx] = None
            except KeyError:
                try:
                    formatter = self.optional_parameter_formatters[k]
                    if formatter is not None:
                        parameters[idx] = formatter.format(val)
                    else:
                        parameters[idx] = None
                except KeyError:
                    parameters[idx] = val

        return [p for p in parameters if p is not None]

    def format_in_path(self, name, path):
        if self.is_local:
            if not os.path.isfile(path):
                raise RuntimeError("{} is not found".format(path))
            return os.path.abspath(path)
        else:
            new_path = os.path.join("/data", os.path.basename(path))
            local_path = os.path.abspath(os.path.join(self.work_dir, os.path.basename(path)))
            if not os.path.abspath(os.path.dirname(path)) == os.path.abspath(self.work_dir):
                shutil.copy(path, local_path)
            return new_path

    def format_out_path(self, name, path):
        path = os.path.abspath(path)

        if not os.path.dirname(path):
            raise RuntimeError("parent dir of {} is not found".format(path))

        if self.is_local:
            self.change_paths[path] = path
            return os.path.abspath(path)
        else:
            new_path = os.path.join("/data", os.path.basename(path))
            fix_path = os.path.join(self.work_dir, os.path.basename(path))
            self.change_paths[fix_path] = path
            return new_path

    def format_in_std_in(self, name, path):
        if not self.is_local:
            return self.format_in_path(name, path)

        self.stdin = path

        return None

    def store_true(self, name, value):
        if value:
            return "--{}".format(name)

    def check_output(self):
        out_files = []
        for out_path, change_path in self.change_paths.items():
            assert os.path.isfile(out_path), "Outfile '{}' not found. ".format(
                out_path) + "Could not change to '{}': {}".format(
                    change_path, os.listdir(self.work_dir))
            if os.path.abspath(out_path) != os.path.abspath(change_path):
                shutil.move(out_path, change_path)
            out_files.append(change_path)

        if len(out_files) == 1:
            return out_files[0]

        return out_files

    def clean(self):
        while len(self.files_to_remove) > 0:
            try:
                os.remove(self.files_to_remove.pop())
            except OSError:
                pass

    def tempfile(self):
        tfile = tempfile.NamedTemporaryFile(dir=self.work_dir, delete=False)
        tfile.close()
        return tfile.name

    def kill(self):
        containerKill()

    def stop(self):
        containerStop()

    def isRunning(self):
        containerIsRunning()
