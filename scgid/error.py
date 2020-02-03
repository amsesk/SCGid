import sys
import traceback
import warnings
from scgid.modcomm import get_error_handler

def check_result(func):
    def wrapper(*args, **kwargs):
        ret = func(*args, **kwargs)
        if is_ok(ret):
            return ret.values

        elif isinstance(ret, FatalError):
            return ret

        else:
            warnings.warn(f"Value returned by `{func.__name__}` is not a result. This behavior is being deprecated.", DeprecationWarning)
            return ret 

    return wrapper

class Ok():
    def __init__(self, values = None):
        self.values = values

    def unpack(self):
        return self.values
        

def is_ok(retval):
    if isinstance(retval, Ok):
        return True

    else:
        return False

class FatalError(Exception):
    def __init__(self):
        self.__name__ = type(self).__name__

        # Get call stack reference to the entity that typically inherits scgid.module.Module that will be printing errors
        self.handler = get_error_handler()

        # Exit status (self.errno) set to 1 generally, but should be reset in child superclasses
        self.errno = 1

        # Message text is set to something useful in terminal child class
        self.msg = "Fatal error occurred. If you are seeing this, there is an unset self.msg in called child class that needs to be set. Please report this."
    
    def catch(self):
        self.handler.logger.critical(f"[{self.__name__}] {self.msg}")
        sys.exit(self.errno)

class ConfigError(FatalError):
    def __init__(self, msg):
        super().__init__()
        self.msg = msg
        self.errno = 2
        self.catch()

class MissingDependencyError(FatalError):
    def __init__(self, msg):
        super().__init__()
        self.msg = msg
        self.errno = 3
        self.catch()

class ErrorClassNotImplemented(FatalError):
    def __init__(self):
        super().__init__()
        self.errno = 55
        self.msg = "You are seeing this because a novel error has occured. Please report this as an issue at `https://www.github.com/amsesk/SCGid.git`"
        self.catch()

class ModuleError(FatalError):
    def __init__(self):
        super().__init__()
        self.errno = 4
        if type(self).__name__ == "ModuleError":
            self.catch()
        
        # Child error class should be catching itself
        else:
            pass

class ArgumentError(FatalError):
    def __init__(self, msg):
        super().__init__()
        self.msg = msg
        self.errno = 5
        self.catch()

class InternalError(FatalError):
    def __init__(self, msg):
        super().__init__()
        self.errno = 99
        self.msg = msg
        self.catch()

        
