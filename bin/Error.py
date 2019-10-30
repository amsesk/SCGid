import sys
import traceback

class Error(Exception):
    pass
class MissingDependencyError(Error):
    def __init__(self, unavail):
        self.unavail = unavail
        self.exitcode = 7
    def __str__(self):
        return "Required dependency missing from environment: {}".format( ", ".join(self.unavail) )

class InternalError(Error):
    def __init__(self):
        self.exitcode = 1
    def __str__(self):
        return ""

def iserror(res):
    if isinstance(res, Error):
        return True
        
