import inspect
import logging

# Log naming via callstack
class LoggingEntity(object):
    def __init__(self):
        pass

class Head:
    def __init__(self):
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

# Handlers
class ExitOnExceptionHandler(logging.StreamHandler):
    def emit(self, record):
        if record.levelno == logging.CRITICAL:
            raise SystemExit(-1)