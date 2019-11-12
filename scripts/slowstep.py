# Class to handle slow steps that need only be done once per SCGid run
import re
import inspect
class ReusableOutput(object):
    def __init__(self, arg, pattern, generate):
        self.arg = arg
        self.parent = None
        self.re_pattern = re.compile(pattern)
        self.needs_doing = self.is_present(self.parent.prefix)
        self.generate = generate

        print inspect.stack()

    def is_present(self, prefix):
        pass

    def generate(self):
        return self.genfunc()