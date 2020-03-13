import sys
import os
import yaml
import inspect
import shutil
import logging
import subprocess
import pkg_resources
import warnings
from datetime import datetime
from scgid.config import Config
from scgid.dependencies import Dependencies
from scgid.reuse import ReusableOutputManager
from scgid.modcomm import logger_name_gen, LoggingEntity, ErrorHandler, get_head, get_root
from scgid.error import ConfigError, ArgumentError, ModuleError, check_result, Ok
from scgid.parsers import PathStore
from scgid.library import get_logging_config

class Module (object):
    def __init__(self, call, name = None, parent=None, loglevel=logging.INFO):
        self.call = call
        self.name = call
        self.root = get_root()
        self.wd = None
        self.argparser = None
        self.parsed_args = None
        self.raw_args = sys.argv
        if name is not None:
            self.name = name

        self.logger = None
        self.loglevel = loglevel

        self.start_logging()

        self.config = Config()
        self.config.logfile_path = None
        self.config.runid = self.root.start_ts

        self.set_module_logging_level(self.loglevel)

        self.config.load_yaml()

    def log_config (self):
        self.logger.info(f"The logfile for this run is located at `{os.path.join(self.config.rundir, 'logfiles', self.config.logfile_path)}`.")
        self.simplelogger.info (f"\n`scgid {type(self).__name__.lower()}` was invoked with the following call:\n{'-'*50}")
        self.simplelogger.info(f"{' '.join(self.raw_args)}\n")
        self.simplelogger.info(f"The full SCGid configuration for run {self.config.runid} is as follows:\n{'-'*50}\n{self.config}\n{'-'*50}\n")

        self.logger.info(f"Initilization and configuration complete.")
        self.logger.info(f"Starting `scgid {type(self).__name__.lower()}`")

    # Defined by subclasses
    def generate_argparser(self):
        return None

    def check_path_args(self):
        warnings.warn("Using PathStore (versus PathAction) in argparser nullifies that need for this function to be run at start of module.", DeprecationWarning)
        for arg in [v for v in self.argparser.__dict__["_actions"] if isinstance(v, PathAction)]:
            path = getattr(self.parsed_args, arg.dest)
            if not os.path.isfile(path):
                return ArgumentError(f"Value given for `{'|'.join(arg.option_strings)}` does not exist: {path}")
            else:
                pass
        return None

    @check_result
    def log_to_rundir(self, name):
        loggers = [logging.getLogger(x) for x in ["SCGid", "data"]]
        
        if not os.path.isdir(os.path.join(self.config.get("rundir"), "logfiles")):
            os.mkdir(os.path.join(self.config.get("rundir"), "logfiles"))

        for l in loggers:
            for h in l.handlers:
                h.close()
            l.handlers = []

        now = datetime.now().strftime("%Y-%m-%d_%H:%M:%S")

        logfile_path = os.path.join(self.config.get("rundir"), "logfiles", f"SCGid.{name}.{self.root.start_ts}.log")
        self.config.logfile_path = logfile_path
        
        try:
            shutil.move(os.path.join(self.config.get("rundir"), "../", "scgid.log.tmp"), logfile_path)
        except FileNotFoundError:
            pass

        cfg = get_logging_config( logfile_path )

        logging.config.dictConfig( cfg )

        return Ok(None)
    
    def start_logging(self):
        self.logger = logging.getLogger( logger_name_gen() )
        self.simplelogger = logging.getLogger("data")

        '''
        logfiles_path = os.path.join(
                os.getcwd(),
                self.config.get("prefix"),
                "logfiles"
                )
        '''
        #if not os.path.isdir(logfiles_path):
            
    def set_module_logging_level(self, loglevel = logging.INFO):
        ident = type(self).__name__
        downstream_loggers = [logging.getLogger(name) for name in logging.root.manager.loggerDict if ident in name]
        for dsl in downstream_loggers:
            dsl.setLevel(loglevel)

    def set_rundir (self, prefix):
        rundir = "{}{}".format(prefix, self.config.OUTPUTSUFFIX)
        try:
            os.chdir(rundir)
        except:
            os.mkdir(rundir)
            os.chdir(rundir)
            self.logger.info("Creating directory `%s`", os.getcwd())
        
        self.config.rundir = os.getcwd()

    def setwd(self, name):
        name = name.split(".")[1]

        try:
            os.chdir(name)
        except:
            os.mkdir(name)
            os.chdir(name)
            self.logger.info("Creating directory `%s`", os.getcwd())
        if os.path.isdir('temp'):
            shutil.rmtree('temp')
        os.mkdir('temp')
        os.chdir('temp')
        self.logger.info("Entering temporary working directory: `%s`", os.getcwd())
        self.wd = os.getcwd()

    def migrate_temp_dir(self):
        for output_file in os.listdir('.'):
            shutil.move(output_file, f"../{output_file}")
        os.chdir("../")
        os.rmdir("temp")
    
    def resetwd(self):
        os.chdir("../../")

    # This method takes a dictionary or arguments typically read from a .yaml file (argdict) generated by SCGidPipeline.create_options_file()
    # and a reference to the module's instance of argparse.ArgumentParser (parser - generated by Module.generate_argparser()). Passing this reference
    # makes it easier to centralize this script in scgid.module and avoid code duplication.
    def translate_argdict(self, argdict, parser):
        translated_argdict = {}
        for a in parser._actions:

            # Skip positional arguments
            if len(a.option_strings) == 0:
                continue

            # Grab the long-form option_string, with the `--` before it. This is how it is referenced in the object returned by argparse.ArgumentParser.parsed_args()    
            key = a.option_strings[-1]
            assert key.startswith("--"), "Missing `--blah` argument option string defined in argparse.ArgumentParser associated with this module."
            key = key.replace("-","")


            # There should only be one key for each arugment in the argdict...
            # This should be impossible because of`iter( set(itertools.chain(gct_vars, kmers_vars, codons_vars)) )` in SCGidPipeline.create_options_file()
            if key in translated_argdict:
                raise KeyError("Key duplication in argdict.")
                sys.exit(1)

            else:
                if a.metavar in argdict:
                    value = argdict[a.metavar]
                    
                    # Set whitespace values, str(None) values, and values of zero length to python None. Skip for integer arguments
                    if not isinstance(value, int) and value is not None:
                        if (
                            value.isspace() or 
                            len(value) == 0 or 
                            value == "None"
                            ):

                            value = None

                    # Assert that all required arguments are not None
                    if (
                        a.required and
                        value is None
                        ):

                        return ArgumentError(f"Required option {'|'.join(a.option_strings)} is missing from options file.")
                    
                    # Add value to the argdict, making sure to add the absolute path if needed
                    if (
                        value is not None and 
                        isinstance(a, PathStore)
                        ):

                        translated_argdict[key] = os.path.abspath(value)

                    else:    
                        translated_argdict[key] = value

                else:

                    # This warning will print if a metavar in the argdict generated from argparse.ArgumentParser._actions is not present in argdict. 
                    # This happens with `--help`.
                    if self.root is not None:
                        self.root.logger.debug(f"Problem translating key `{key}`. The key does not occur in argdict")
        
        return translated_argdict

    def cli_invocation(self):
        set_opts = {k.strip():v.strip() for k,v in self.translated_args.items() if v is not None}
        optstr = [f"--{k} {v}" for k,v in set_opts.items()]
        optstr.insert(0, type(self).__name__.lower())
        optstr.insert(0, "scgid")
        return optstr
            

    '''
    def try_catch_error(self, op):
        if Error.iserror(op):
            self.handle_error(op)
        else:
            return op
    def handle_error(self, Err):
        ErrType = type(Err).__name__
        if Err.is_fatal():
            print "[FATAL {}] {}".format(ErrType, Err)
            sys.exit(Err.exitcode)
        else:
            print "[{}] {}".format(ErrType, Err)
            if not Err.recover():
                print "[FATAL {}] {}".format(ErrType, Err)
                sys.exit(Err.exitcode)
    '''

            