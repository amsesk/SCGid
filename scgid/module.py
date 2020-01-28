import sys
import os
import yaml
import inspect
import shutil
import logging
import subprocess
import pkg_resources
from scgid.config import Config
from scgid.dependencies import Dependencies
from scgid.reuse import ReusableOutputManager
from scgid.modcomm import logger_name_gen, LoggingEntity, ErrorHandler, get_head
from scgid.error import ConfigError, ArgumentError, ModuleError
from scgid.parsers import PathAction

class Module (object):
    def __init__(self, call, name = None, parent=None):
        self.call = call
        self.name = call
        self.caller = inspect.stack()[1][0].f_locals["self"]
        self.wd = None
        self.argparser = None
        self.parsed_args = None
        if name is not None:
            self.name = name
        self.logger = None

        #Set by self.initialize()
        self.config = Config()
    
    def generate_argparser(self):
        pass

    def check_path_args(self):
        warnings.warn("Using PathStore (versus PathAction) in argparser nullifies that need for this function to be run at start of module.", DeprecationWarning)
        for arg in [v for v in self.argparser.__dict__["_actions"] if isinstance(v, PathAction)]:
            path = getattr(self.parsed_args, arg.dest)
            if not os.path.isfile(path):
                return ArgumentError(f"Value given for `{'|'.join(arg.option_strings)}` does not exist: {path}")
            else:
                pass
        return None
    
    def start_logging(self):
        self.logger = logging.getLogger( logger_name_gen() )
        self.simplelogger = logging.getLogger("SCGid.unfmt")

    def setwd (self, name, prefix):
        name = name.split(".")[1]
        rundir = "{}{}".format(prefix, self.config.OUTPUTSUFFIX)
        try:
            os.chdir(rundir)
        except:
            os.mkdir(rundir)
            os.chdir(rundir)
            self.logger.info("Creating directory `%s`", os.getcwd())
        
        self.config.rundir = os.getcwd()
        
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

    def translate_argdict(self, argdict, parser):
        modclass = (type(self).__name__)
        #actions = eval(modclass).generate_argparser()._actions
        actions = parser._actions
        translated_argdict = {}
        for a in actions:
            if len(a.option_strings) == 0:
                continue
            elif len(a.option_strings) == 2:
                key = a.option_strings[1]
            else:
                key = a.option_strings[0]
            assert key.startswith("--"), "Bad argument key"
            key = key.replace("-","")

            if key in translated_argdict:
                raise KeyError("Key duplication in options file.")
                sys.exit(1)

            else:
                if a.metavar in argdict:
                    value = argdict[a.metavar]
                    if value.isspace() or len(value) == 0 or value == "None":
                        value = None
                    translated_argdict[key] = value

                else:
                    print(f"Problem translating key `{key}`. The key does not occur in argdict")
        
        return translated_argdict

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

            