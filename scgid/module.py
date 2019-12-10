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

            