import sys
import os
import yaml
from lib import SimpleNamespace
import inspect
import shutil
import Error
import logging
from Depends import Dependencies

class Module (object):
    def __init__(self, call, name = None, parent=None):
        self.call = call
        self.name = call
        self.parent = parent
        self.OUTPUTSUFFIX = "_scgid_output"
        if name is not None:
            self.name = name
        self.log = None

        #Set by self.initialize()
        self.config = Config()
    def start_logging(self):
        self.log = logging.getLogger("{}.{}".format(
            self.parent.__name__,
            ".".join(self.config.get("mod"))
            ))

    def setwd (self, name, prefix):
        rundir = "{}{}".format(prefix, self.OUTPUTSUFFIX)
        try:
            os.chdir(rundir)
        except:
            os.mkdir(rundir)
            os.chdir(rundir)
            self.log.info("Creating directory `%s`", os.getcwd())
        try:
            os.chdir(name)
        except:
            self.log.info("Creating directory `%s`", name)
            os.mkdir(name)
            os.chdir(name)
            self.log.info("Creating directory `%s`", os.getcwd())
        if os.path.isdir('temp'):
            shutil.rmtree('temp')
        os.mkdir('temp')
        os.chdir('temp')
        self.log.info("Entering working directory: `%s`", os.getcwd())

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

class Config(object):
    def __init__(self):
        self.SCGIDBIN = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
        self.SCGID = os.path.dirname(self.SCGIDBIN)
        self.dependencies = Dependencies()
        self.__dict__.update(self.load_yaml(self.SCGID))

    def __repr__(self):
        return "\n".join(["{}: {}".format(key, setting) for key,setting in self.__dict__.iteritems()])

    def get(self, value):
        return getattr(self, value)

    def load_cmdline(self, parsed_args):
        for unset in [x for x,p in parsed_args.__dict__.iteritems() if p is None]:
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
                return yaml.load(cfg)
        except IOError:
            raise Error.MissingConfigError("Unable to locate config.yaml")
        except:
            raise Error.MalformedConfigError("Malformed config.yaml")
class Kmers(Module):
    def __init__(self):
        super(Kmers, self).__init__("kmers")
class Train(Kmers):
    def __init__(self):
        super(Train, self).__init__()
            