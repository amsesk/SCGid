import sys
import os
import yaml
import inspect
import shutil
import logging
from scripts.dependencies import Dependencies
from scripts.reuse import ReusableOutputManager
from scripts.modcomm import logger_name_gen, LoggingEntity
from scripts.error import MissingConfigError, MalformedConfigError

class Module (object):
    def __init__(self, call, name = None, parent=None):
        self.call = call
        self.name = call
        self.caller = inspect.stack()[1][0].f_locals["self"]
        self.wd = None
        if name is not None:
            self.name = name
        self.logger = None

        #Set by self.initialize()
        self.config = Config()
    def start_logging(self):
        self.logger = logging.getLogger(
            logger_name_gen()
        )

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
            self.logger.info("Creating directory `%s`", name)
            os.mkdir(name)
            os.chdir(name)
            self.logger.info("Creating directory `%s`", os.getcwd())
        if os.path.isdir('temp'):
            shutil.rmtree('temp')
        os.mkdir('temp')
        os.chdir('temp')
        self.logger.info("Entering working directory: `%s`", os.getcwd())
        self.wd = os.getcwd()

    def initialize(self):
        pass

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

class Config(LoggingEntity):
    def __init__(self):
        self.SCGID_SCRIPTS = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
        self.SCGID = os.path.dirname(self.SCGID_SCRIPTS)
        self.OUTPUTSUFFIX = "_scgid_output"
        self.dependencies = Dependencies()
        self.reusable = ReusableOutputManager()
        self.__dict__.update(self.load_yaml(self.SCGID))
        self.rundir = None

    def __repr__(self):
        return "\n".join(["{}: {}".format(key, setting) for key,setting in self.__dict__.items()])

    def get(self, value):
        return getattr(self, value)

    def load_cmdline(self, parsed_args):
        for unset in [x for x,p in parsed_args.__dict__.items() if p is None]:
            try:
                setattr(parsed_args, unset, getattr(self, "DEFAULT_{}".format(unset)))
            except:
                pass
        self.__dict__.update( vars(parsed_args) )

    def load_yaml(self, HOME):
        loc = os.path.join(HOME, "config.yaml")
        if os.path.isfile( "{}.local".format(loc) ):
            loc = "{}.local".format(loc)
        try:
            with open(loc, 'r') as cfg:
                return yaml.load(cfg, Loader=yaml.BaseLoader)
        except IOError:
            raise MissingConfigError("Unable to locate config.yaml")
        except:
            raise MalformedConfigError("Malformed config.yaml")

            