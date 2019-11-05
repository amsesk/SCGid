import sys
import os
import yaml
from types import SimpleNamespace
import inspect
import Error
from Depends import Dependencies

class Module ():
    def __init__(self, call, pargs, name = None):
        self.call = call
        self.name = call
        if name is not None:
            self.name = name
        self.pargs = pargs
        self.bin_dir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
        self.home = os.path.dirname(self.bin_dir)

        #Set by self.initialize()
        self.dependencies = Dependencies(self.pargs)
        self.config = None

    def initialize(self):
        self.config = self.get_config()
        self.try_catch_error( self.dependencies.check() )

    def get_config(self):
        with open( os.path.join(self.home, "config.yaml") ) as config:
            return yaml.load(config)

    def set_dependencies(self, *deps):
        self.dependencies = Dependencies(self.pargs, *deps)

    def try_catch_error(self, op):
        if Error.iserror(op):
            self.die_on_error(op)

    def die_on_error(self, Err):
        print "[FATAL {}] {}".format(type(Err).__name__, Err)
        sys.exit(Err.exitcode)