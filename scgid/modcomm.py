import inspect
import logging
import sys
import os

# Log naming via callstack
class Root:
    pass
    
class LoggingEntity(object):
    pass

class Head:
    pass

class ErrorHandler:
    pass

def get_caller ():
    for f in inspect.stack():
        try:
            return f[0].f_locals["self"]
        except:
            pass # Not a module, with a `self`
def get_callstack ():
    curr = inspect.currentframe()
    callstack = [curr]
    while True:
        if curr.f_back is None:
            break
        else:
            curr = curr.f_back
            callstack.append(curr)
    return callstack
def callstack_clsfilt (callstack):
    cls_callstack = []
    for l in callstack:
        try:
            _ = l.f_locals["self"]
            cls_callstack.append(l)
        except:
            pass
    return cls_callstack
def logger_name_gen():
    callstack = callstack_clsfilt( 
            get_callstack()
        )
    cls_stack = [x.f_locals["self"] for x in callstack]
    cls_stack_dedup = []

    # Dedepulicate list
    for i in cls_stack:
        if i not in cls_stack_dedup:
            cls_stack_dedup.append(i)

    return ".".join(
        #[type(x).__name__ for x in cls_stack_dedup][::-1]
        [type(x).__name__ for x in cls_stack_dedup if isinstance(x, LoggingEntity)][::-1]
    )
def get_head():
    callstack = callstack_clsfilt( 
            get_callstack()
        )
    cls_stack = [x.f_locals["self"] for x in callstack]
    cls_stack_dedup = []

    # Dedepulicate list
    for i in cls_stack:
        if i not in cls_stack_dedup:
            cls_stack_dedup.append(i)

    return [x for x in cls_stack_dedup if isinstance(x, Head)][0]

def get_root():
    callstack = callstack_clsfilt( 
            get_callstack()
        )
    cls_stack = [x.f_locals["self"] for x in callstack]
    cls_stack_dedup = []

    # Dedepulicate list
    for i in cls_stack:
        if i not in cls_stack_dedup:
            cls_stack_dedup.append(i)
    
    # try/except block to allow testing of individual module components
    try:
        return [x for x in cls_stack_dedup if isinstance(x, Root)][0]
    
    except IndexError:
        return None

def get_error_handler():
    callstack = callstack_clsfilt( 
            get_callstack()
        )
    cls_stack = [x.f_locals["self"] for x in callstack]
    cls_stack_dedup = []

    # Dedepulicate list
    for i in cls_stack:
        if i not in cls_stack_dedup:
            cls_stack_dedup.append(i)

    return [x for x in cls_stack_dedup if isinstance(x, ErrorHandler)][0]

def pkgloc():
    source = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
    home = os.path.dirname(source)
    return home, source

# Handlers
class ExitOnExceptionHandler(logging.StreamHandler):
    def emit(self, record):
        if record.levelno == logging.CRITICAL:
            pass