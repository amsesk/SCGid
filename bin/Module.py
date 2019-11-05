import Error
import sys

class Module ():
    def __init__(self, call, dependencies, name = None):
        self.call = call
        self.name = call
        if name is not None:
            self.name = name
        self.dependencies = dependencies
    def initialize(self):
        self.try_catch_error( self.dependencies.check() )
    def try_catch_error(self, op):
        if Error.iserror(op):
            self.die_on_error(op)
    def die_on_error(self, Err):
        print "[FATAL {}] {}".format(type(Err).__name__, Err)
        sys.exit(Err.exitcode)