import sys
import traceback
from scgid.modcomm import get_error_handler

class FatalError(Exception):
    def __init__(self, msg):
        self.__name__ = type(self).__name__
        self.errno = 1
        self.handler = get_error_handler()
    
    def catch(self):
        self.handler.logger.critical((f"[{self.__name__}] {self}"))
        sys.exit(self.errno)

class ConfigError(FatalError):
    def __init__(self, msg):
        super().__init__(msg)
        self.errno = 2
        self.catch()

class MissingDependencyError(FatalError):
    def __init__(self, msg):
        super().__init__(msg)
        self.errno = 3
        self.catch()

class ModuleError(FatalError):
    def __init__(self, msg):
        super().__init__(msg)
        self.errno = 4
        self.catch()

class InternalError(FatalError):
    def __init__(self, msg):
        super().__init__(msg)
        self.errno = 99
        self.catch()

        
