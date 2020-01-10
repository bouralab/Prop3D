import os
import sys
import subprocess

if "USE_SINGULARITY" in os.environ:
    USE_SINGULARITY = os.environ["USE_SINGULARITY"][0].lower()=="t"
else:
    USE_SINGULARITY = subprocess.check_output(["singularity"]).starts_with(
        "USAGE: singularity")

if not USE_SINGULARITY:
    from molmimic.parsers.singularity import apiSingularityCall as containerCall
    from molmimic.parsers.singularity import singularityKill as containerKill
    from molmimic.parsers.singularity import singularityStop as containerStop
    from molmimic.parsers.singularity import containerIsRunning
else:
    from toil.lib.docker import apiDockerCall as containerCall
    from toil.lib.docker import dockerKill as containerKill
    from toil.lib.docker import dockerStop as containerStop
    from toil.lib.docker import containerIsRunning

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
    RETURN_FILES = False

    def __init__(self, job=None, return_files=False, force_local=False, fallback_local=False, work_dir=None):
        assert (self.IMAGE, self.LOCAL).count(None) != 2, "Must define conatiner or local path"

        if self.LOCAL is not None:
            assert isinstance(self.LOCAL, list)

        self.work_dir = work_dir if work_dir is None else os.getcwd()
        self.force_local = force_local or self.IMAGE is None
        self.fallback_local = fallback_local
        self.return_files = return_files or self.RETURN_FILES
        self.job = job
        self.rules = {
            "str": str,
            "path:in": self.format_in_path,
            "path:out": self.format_out_path,
            "path:in:stdin": "format_in_std_in",
            "store_true": "store_true"}.update(self.RULES)

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

        self.number_of_parameters = 0
        self.number_of_optional_parameters = 0

        self.change_paths = {}

        self.is_local = False
        self.stdin = None

        for i, p in enumerate(self.PARAMETERS):
            if isinstance(p, str):
                self.paramters.append(p)
            elif isinstance(p, (list, tuple)):
                if len(p) == 1 and isinstance(p[0], str):
                    #ignore, treat as string
                    paramters.append(p)
                elif len(p) in [2,3] and all(isinstance(_p, str) for _p in p):
                    key, rule = p[0], p[1]

                    optional = False
                    if key.startswith(":"):
                        #Optional value
                        key, default = key[1:].split(":")
                        if len(default)==0 and rule == "store_true":
                            default = None
                        else:
                            raise RuntimeError("invalid optional parameters: {}".format(p))

                        optional = True

                    if rule in self.rules:
                        if callable(self.rules[rule]):
                            func = self.rules[rule]
                        elif hasattr(self, rule):
                            func = getattr(self, rule)



                        if not optional:
                            self.parameters.append(None)
                            self.param_names.append(key)
                            self.param_funcs[i] = func
                            self.params_to_update[key] = i
                            self.number_of_parameters += 1
                        else:
                            self.parameters.append(default)
                            self.optional_param_names.append(key)
                            self.optional_param_funcs[i] = func
                            self.optional_params_to_update[key] = i
                            self.number_of_optional_parameters += 1
                            self.optional_parameters_defaults[key] == default

                        if len(p) == 3 and:
                            formatter = p[3]
                            self.parameter_formatters[key] = formatter
                    else:
                        raise RuntimeError("invalid parameters: {}".format(p))

    def __call__(self, *args, **kwds):
        if self.force_local:
            return self.local(parameters)

        assert self.job is not None

        parameters = self.format_parameters(self, args, kwds)

        try:
            out = containerCall(job,
                          image=self.IMAGE,
                          working_dir="/data",
                          volumes={self.work_dir:{"bind":"/data", "mode":"rw"}},
                          parameters=parameters)
        except (SystemExit, KeyboardInterrupt):
            raise
        except:
            if self.fallback_local:
                import traceback as tb
                RealtimeLogger.error(tb.format_exc())
                return self.local(parameters)
            raise

        out_files = self.check_output()

        if self.return_files:
            return out_files

        return out

    def local(self, *args, **kwds):
        self.is_local = True
        parameters = self.format_parameters(self, args, kwds)

        env = self.set_local_env()

        try:
            if self.stdin is not None:
                with open(self.stdin) as f:
                    out = subprocess.check_output(self.LOCAL+parameters, env=env, stdin=f)
            else:
                out = subprocess.check_output(self.LOCAL+parameters, env=env)
        except (SystemExit, KeyboardInterrupt):
            raise
        except:
            raise
        self.is_local = False
        self.stdin = None

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
            kwds = zip(self.param_names, args).update(kwds)
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

        # assert len(args)+len(kwds) == self.number_of_parameters, "Must use paramters: {}".format(", ".join(self.param_names))
        # if len(kwds) == self.number_of_parameters:
        #     pass
        # elif len(args) == self.number_of_parameters:
        #     kwds = zip(self.param_names, args)
        # else:
        #     kwds = zip(self.param_names, args).update(kwds)

        self.change_paths = {}

        parameters = self.parameters[:]
        for k, v in kwds.items():
            idx = self.params_to_update[k]
            val = self.param_funcs[k](k, v)

            if val is None:
                continue

            try:
                formatter = self.parameter_formatters[k]
                if formatter is not None:
                    parameters[idx] = formatter.format(val)
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
        if not os.path.dirname(path):
            raise RuntimeError("parent dir of {} is not found".format(path))

        if self.is_local:
            return os.path.abspath(path)
        else:
            new_path = os.path.join("/data", os.path.basename(path))
            fix_path = os.path.join(self.work_dir, os.path.basename(path))
            self.change_paths[fix_path] = new_path
            return new_path

    def format_in_std_in(self, name, path):
        if not self.is_local:
            return self.format_in_path(path)

        self.stdin = path

        return None

    def store_true(self, name, value):
        if value:
            return "--{}".format(name)

    def check_output(self):
        out_files = []
        for out_path, change_path in self.change_paths.items():
            assert os.path.isfile(out_path), "Outfile not found: {}".format(os.listdir(work_dir))
            shutil.move(out_path, change_path)
            out_files.append(change_path)

        if len(out_files) == 1:
            return out_files[0]

        return out_files


    def kill(self):
        containerKill()

    def stop(self):
        containerStop()

    def isRunning(self):
        containerIsRunning()
