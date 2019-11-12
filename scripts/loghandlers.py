import logging

class ExitOnExceptionHandler(logging.StreamHandler):
    def emit(self, record):
        if record.levelno == logging.CRITICAL:
            raise SystemExit(-1)