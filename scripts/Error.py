import sys
import traceback

class Error(Exception):
    def __init__(self):
        self.fatal = False
    def is_fatal(self):
        return self.fatal
    def recover(self):
        return None

class MalformedConfigError(Error):
    def __init__(self, message):
        super(MalformedConfigError, self).__init__()
        self.message = message
        self.exitcode = 2
        self.fatal = False
    def recover(self):
        print "Recover failed"
        return False
    def __str__(self):
        return self.message
class MissingConfigError(Error):
    def __init__(self, message):
        super(MissingConfigError, self).__init__()
        self.message = message
        self.exitcode = 2
        self.fatal = False
    def recover(self):
        print "Recover failed"
        return False
    def __str__(self):
        return self.message

class MissingDependencyError(Error):
    def __init__(self, unavail):
        super(MissingDependencyError, self).__init__()
        self.unavail = unavail
        self.exitcode = 7
        self.fatal = True
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
        
