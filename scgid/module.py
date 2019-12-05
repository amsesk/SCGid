import sys
import os
import yaml
import inspect
import shutil
import logging
import subprocess
from scgid.dependencies import Dependencies
from scgid.reuse import ReusableOutputManager
from scgid.modcomm import logger_name_gen, LoggingEntity, ErrorHandler
from scgid.error import ConfigError

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
    
    def generate_argparser(self):
        pass
    
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

class Config(LoggingEntity, ErrorHandler):
    def __init__(self):
        self.logger = logging.getLogger( logger_name_gen() )
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
        try: 
            attr = getattr(self, value)
            return attr

        except KeyError:
            raise KeyError(f"Member `{value}` not in Config")
        

    def load_cmdline(self, parsed_args, case_args={}):
        for unset in [x for x,p in parsed_args.__dict__.items() if p is None]:
            try:
                setattr(parsed_args, unset, getattr(self, unset))
            except:
                self.logger.debug(f"No loaded value from config.yaml to set None-type argument `{unset}`")
        
        self.__dict__.update( vars(parsed_args) )
    
    def load_argdict(self, argdict):
        self.__dict__.update(argdict)

    def load_yaml(self, HOME):
        loc = os.path.join(HOME, "config.yaml")
        if os.path.isfile( "{}.local".format(loc) ):
            loc = "{}.local".format(loc)
        try:
            with open(loc, 'r') as cfg:
                return yaml.load(cfg, Loader=yaml.BaseLoader)
        except IOError:
            return ConfigError("Unable to locate configuration file")
        except:
            return ConfigError("Malformed configuration file")
    
    def check_case_args (self, mapper_dict):
        for option_arg, options in mapper_dict.items():
            choice = options[self.get(option_arg)]
            if not choice["bool"]:
                self.logger.warning(choice["warning"])

    def check_esom_path (self):
        # Check to see if esomtrn exists in configured `esom_path`
        if not os.path.isfile(
            os.path.join(
                self.get("esom_path"),
                "esomtrn"
            )):
            # If not, try to find it in environment and update config
            p = subprocess.Popen( ["which", "esomtrn"], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
            path, _ = p.communicate()
            if p.returncode:
                return ConfigError("ESOM path unset in $PATH or config.yaml")

            else:
                path = os.path.dirname(path).decode("utf-8")
            
                # Updating config here
                setattr(self, "esom_path", path)

            