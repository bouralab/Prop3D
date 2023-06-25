import os
import json
from shutil import copyfileobj
from functools import partial
import urllib.request
import requests
from contextlib import closing

#from toil.realtimeLogger import RealtimeLogger

class RealtimeLogger:
    #staticmethod
    def info(*args, **kwds):
        print(args, kwds)

class WebService(object):
    def __init__(self, base_url, store, work_dir=None, download=True, clean=True, max_attempts=2):
        self.work_dir = os.getcwd() if work_dir is None else work_dir
        self.base_url = base_url
        self.store = store
        self._download = True
        self._clean = clean
        self.max_attempts = max_attempts
        #self.get = memory.cache(self.get)
        self.files = {}

    def __get__(self, key):
        return self.get(key)

    def __del__(self):
        self.clean()

    def clean(self):
        if not hasattr(self, 'files'):
            return 
            
        for key, (fname, should_delete) in self.files.items():
            if should_delete:
                try:
                    os.remove(fname)
                except (OSError, FileNotFoundError):
                    pass

    def fix_key(self, key):
        return key

    def extension(self, key=None):
        return ""

    def parse(self, file_path, key):
        with open(file_path) as fh:
            return fh.read()

    def check_line(self, key, line, attempts):
        """Check line to test for failure. Return True if should rerun.
        Subclass to add functionality."""
        return False

    def __get__(self, key):
        return self.get(key)

    def get(self, key, attempts=None, last_source=None):
        if attempts is None:
            attempts = self.max_attempts

        if isinstance(key, (list, tuple)):
            key = "/".join(self.fix_key(str(k)) for k in key)
        elif isinstance(key, str):
            key = self.fix_key(key)
        else:
            raise KeyError(key)

        store_key = "{}{}".format(key, self.extension(key))
        fname = os.path.join(self.work_dir, "{}-{}{}".format(
            self.store.store_name.replace("/", "-"), key.replace("/", "-"),
            self.extension(key)))

        #Check if file should be download or get from store
        if os.path.isfile(fname) and not last_source=="local":
            #File already exists or previosly downloaded in this session
            RealtimeLogger.info("API read from file")
            source = "local"
            should_remove = False #If previosly downloaded in this session, it will not remove
        elif self.store.exists(key) and not last_source=="IOStore":
            RealtimeLogger.info("API get from store")
            self.store.read_input_file(key, fname)
            source = "IOStore"
            should_remove = True
        else:
            should_remove = True
            if self._download:
                rc = self.download(key, fname)
                source = "webservice"
                if not rc:
                    RealtimeLogger.info("DOWNLOAD API -- FAILED")
                    #Error no EPPIC file
                    raise KeyError("Key '{}' does not exist in web server".format(key))
            else:
                RealtimeLogger.info("DOWNLOAD API -- NO DOWNLOAD")
                raise KeyError("Key '{}' does not exist".format(key))

        #Check if file contains data -- important because some files contain error messages
        try:
            with open(fname) as f:
                pass
        except IOError as e:
            #Might be empty
            try:
                os.remove(fname)
            except OSError:
                pass

            RealtimeLogger.info("Failed reading, {} bc {}".format(fname, e))

            if attempts > 0:
                return self.get(key, attempts=attempts-1, last_source=source)
            else:
                raise KeyError("Key '{}' is an invalid file".format(key))
        RealtimeLogger.info("Donwlaod step 5")
        try:
            result = self.parse(fname, key)
        except (SystemExit, KeyboardInterrupt) as e:
            raise
        except ValueError as e:
            rerun = False
            RealtimeLogger.info("Donwlaod step pre6")
            try:
                with open(fname) as f:
                    for line in f:
                        RealtimeLogger.info(f"check line {line}")
                        rerun = self.check_line(key, line, attempts)
            except Exception:
                rerun = True

            try:
                os.remove(fname)
            except (OSError, FileNotFoundError):
                pass

            RealtimeLogger.info("Donwlaod step 6, will rerun={}; attempts={}".format(rerun, attempts))

            if rerun and attempts > 0:
                return self.get(key, attempts=attempts-1, last_source=source)
            else:
                RealtimeLogger.info("Not restarting")
                raise KeyError("Key '{}' is an invalid file".format(key))
        except Exception as e:
            RealtimeLogger.info("API Failed parsing json ({}): {}".format(type(e), e))
            raise KeyError("Key '{}' Not found; {} is an invalid file {} {}".format(key, fname, type(e), e))

        if key not in self.files:
            self.files[key] = (fname, should_remove)

        if should_remove and self._clean:
            self.clean()

        return result

    def download(self, key, fname):
        RealtimeLogger.info("API download from web service")

        url = "{}{}".format(self.base_url, key)

        RealtimeLogger.info("Downloading {}".format(url))

        try:
            if url.startswith("http"):
                with requests.get(url, stream=True) as r, open(fname, 'wb') as f:
                    RealtimeLogger.info("Donwlaod step 1")
                    r.raw.read = partial(r.raw.read, decode_content=True)
                    copyfileobj(r.raw, f)
                    RealtimeLogger.info("Donwlaod step 2")
            elif url.startswith("ftp"):
                with closing(urllib.request.urlopen(url)) as r, open(fname, 'wb') as f:
                    copyfileobj(r, f)
            else:
                raise RuntimeError(f"Invalid protocol for {url}")

        except Exception as e:
            RealtimeLogger.info("API Download Error {}: {}".format(type(e), e))
            return False

        RealtimeLogger.info("Donwlaod step 3")

        if self.should_add_to_store(key, fname):
            #Save to store
            self.store.write_output_file(fname, key)

        RealtimeLogger.info("Donwlaod step 4")

        return True

    def should_add_to_store(self, key, fname):
        return True

class JSONApi(WebService):
    def __init__(self, *args, **kwds):
        super().__init__(*args, **kwds)

    def extension(self, key):
        return ".json"

    def parse(self, file_path, key):
        with open(file_path) as fh:
            return json.load(fh)
