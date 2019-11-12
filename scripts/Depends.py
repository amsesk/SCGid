from collections import namedtuple
from Error import MissingDependencyError
import logging
import itertools
import os

CaseDependencyCouplet = namedtuple("CaseDependencyCouplet", ["argid", "value"])
class Dependency(object):
    def __init__(self, cmd):
        self.cmd = cmd
        self.available = None

    def is_available(self, path_contents):
        if self.cmd in path_contents:
            return True
        else:
            return False

class ConstDependency(Dependency):
    def __init__(self, cmd):
        super(ConstDependency, self).__init__(cmd)

class CaseDependency(Dependency):
    def __init__(self, cmd, argid, value):
        super(CaseDependency, self).__init__(cmd)
        self.couplet = CaseDependencyCouplet(argid, value)


class Dependencies():
    def __init__(self):
        self.deps = ()
        self.log = logging.getLogger(__name__)

    def __repr__(self):
        return ", ".join([x.cmd for x in self.deps])

    def populate (self, *deps):
        self.deps = deps

    def check(self, parsed_args):
        avail = [x for x in os.environ["PATH"].split(":") if os.path.isdir(x)]
        avail = list(itertools.chain.from_iterable([os.listdir(x) for x in avail]))
        for d in self.deps:
            d.available = d.is_available(avail)

        missing = [x.cmd for x in self.deps if isinstance(x, ConstDependency) and not x.available] + [x.cmd for x in self.deps if isinstance(x, CaseDependency) and parsed_args.__dict__[x.couplet.argid] == x.couplet.value and not x.available]
        if len(missing) > 0:
            self.log.critical(MissingDependencyError(missing))
            raise MissingDependencyError(missing)
        else:
            return 0
        
'''
class Depends():
    def __init__(self, deptree=None, *args):
        self.cmds = {x: None for x in args}

    def check(self):
        avail = list(itertools.chain.from_iterable([os.listdir(x) for x in os.environ["PATH"].split(":")]))
        for c, a in self.cmds.iteritems():
            if c in avail:
                self.cmds[c] = True
            else:
                self.cmds[c] = False
        unavail = {c: a for c, a in self.cmds.iteritems() if not a}
        if len(unavail) > 0:
            return Error.MissingDependencyError(unavail)
        else:
            return 0

    def add(self, *args):
        self.cmds.update({x: None for x in args})
'''