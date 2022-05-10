import subprocess
import re
import sys
import os
import pickle
import random
import itertools as it
import logging
from io import StringIO
from collections import OrderedDict
from scgid.sequence import DNASequence, DNASequenceCollection, AASequence, AASequenceCollection

SPDB_OS_REGEXP_PATTERN = "OS=(.+?)(?:OX=.+|GN=.+|PE=.+|SV=.+|$)"
ERASE_LINE = '\x1b[2K'
CURSOR_UP_ONE = '\x1b[1A'
output_cols = {
        'GREEN':'\033[0;32m',
        'RED':'\033[0;31m',
        'RESET':'\033[0m',
        'YELLOW':'\033[93m'
        }
'''
class spinner():
    def __init__(self):
        self.sequence = ['|','/','-','\\','|','/','-','\\']
        self.current = 0
        self.end = 7
    def __iter__(self):
        return self
    def next(self):
        if self.current == self.end:
            self.current = 0
            ow_last_stdout(self.sequence[self.end])
            return str()
        else:
            self.current += 1
            ow_last_stdout(self.sequence[self.current-1])
            return str()
'''
def get_logging_config(logfile_path = "scgid.log.tmp"):
    LOGGING_CONFIG = {
        "version": 1,

        "formatters": {
            "verbose": {
                "format": "%(asctime)s %(name)-12s %(levelname)-8s %(message)s",
            },

            "simple": {
                "format": "%(message)s",
            },
        },

        'handlers': {
            'console_verbose': {
                'level':'INFO',
                'class':'logging.StreamHandler',
                'formatter': 'verbose',
                'stream': 'ext://sys.stdout',
            },
            'file_verbose': {
                'level': 'INFO',
                'class': 'logging.FileHandler',
                'formatter': 'verbose',
                'filename': logfile_path,
                'mode': 'a',
            },
            'console_simple': {
                'level':'INFO',
                'class':'logging.StreamHandler',
                'formatter': 'simple',
                'stream': 'ext://sys.stdout',
            },
            'file_simple': {
                'level': 'INFO',
                'class': 'logging.FileHandler',
                'formatter': 'simple',
                'filename': logfile_path,
                'mode': 'a',
            },
        },

        "loggers": {

            "SCGid": {
                "level": "DEBUG",
                "handlers": ['console_verbose', 'file_verbose'],
                },

            "data": {
                "level": "DEBUG",
                "handlers": ['console_simple', 'file_simple'],

                },
            },
        }

    return LOGGING_CONFIG

class bcolors:
    MAGENTA = '\033[95m'
    CYAN = '\033[96m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def ow_last_stdout(string):
    sys.stdout.write(CURSOR_UP_ONE)
    sys.stdout.write(ERASE_LINE)
    sys.stdout.write(string)


def report_outcome(msg, color, outcome):
    ow_last_stdout('\r'+msg + output_cols[color] + outcome + output_cols['RESET'] + '\n')

def is_fasta(f, strict=False, verbose = False):
    with open(f, 'r') as openf:
        lines = [l.strip() for l in openf.readlines()]

        fail = False
        has_headers = True
        made_up_of_seqs = True
        just_had_header = False
        idx = 0

        #lines = go_through_list(lines)
        for line in lines:
            if not just_had_header:
                if re.match('^>.+',line) is not None:
                    has_headers = True
                    just_had_header = True
                else:
                    has_headers = False
            elif just_had_header:
                if strict:
                    if re.match('^[-atcguATCGUgalmfwkqespvicyhrndtbzxuGALMFWKQESPVICYHRNDTBZXU*]+$',line) is not None:
                        made_up_of_seqs = True
                    else:
                        made_up_of_seqs = False
                else:
                    if re.match('^[-A-Za-z*]+$',line) is not None:
                        made_up_of_seqs = True
                    else:
                        made_up_of_seqs = False
                try:
                    if re.match('^>.+',lines[idx+1]) is not None:
                        just_had_header = False
                except:
                    break
            if not has_headers or not made_up_of_seqs:
                fail = True
                if verbose:
                    print( f"Deviation from fasta format on line #{idx+1}" )
            idx+=1
        if not fail:
            return True
        else:
            return False

def gff3_to_fasta(gff3, outname):
    prots = []
    recording = False
    for line in open(gff3).readlines():
        line = line.strip()
        if recording is False:
            s = re.search("[#][ ]protein[ ]sequence[ ][=][ ]",line)
            if s is not None:
                recording = True
                s = re.search("[#][ ]protein[ ]sequence[ ][=][ ][[]([A-Za-z]+)",line)

                # Case where augustus outputs 0 length protein
                if s is None:
                    recording = False
                    continue
                else:
                    protein_sequence = s.group(1)

                    # Case where protein sequence is entirely on the first line
                    s = re.search("[#][ ]protein[ ]sequence[ ][=][ ][[]([A-Za-z]+)[]]$",line)
                    if s is not None:
                        recording = False

                    # Case where protein sequence is spread over multiple lines
                    else:
                        pass
            else:
                s = re.search("[#][ ]end[ ]gene[ ]([A-Za-z0-9_.]+)",line)
                if s is not None:
                    prots.append(AASequence(s.group(1), protein_sequence))
                    protein_sequence = ""

        elif recording is True:
            s = re.search("[#][ ]([A-Za-z]+)[]]$", line)
            if s is not None:
                recording = False
                protein_sequence += s.group(1)
            else:
                s = re.search("[#][ ]([a-zA-Z]+)",line)
                protein_sequence += s.group(1)
    with open(outname,'w') as fasta:
        for i in prots:
            fasta.write(i.to_fasta()+"\n")

def log_command_line_errors (err, log_inst):
    if err != '':
        print(err)
        return 0
        err = err.split('\n')
        for i in err:
            if len(i) > 0:
                log_inst.critical(i)

def log_command_line_stdout (out, log_inst):
    for l in [x.strip() for x in StringIO(out.decode("utf-8")).readlines()] + [""]: ### Add extra line to deal with weird ESOM OOM stdout that cuts off
        log_inst.info(l)
    return 0

def subprocessC (args):
    p = subprocess.call(args)
    if p != 0:
        sys.exit(1)

def subprocessP (args, log_inst, log_stdout=False):
    p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out,err = p.communicate()
    if log_stdout:
        assert isinstance(log_inst,list), "Internal error. In order to log stdout, a list of logs must be passed to subprocessP"
        log_command_line_stdout(out,log_inst[1])
        if p.returncode != 0:
            log_command_line_errors(err, log_inst[0])
            sys.exit(1)
    else:
        if p.returncode != 0:
            log_command_line_errors(err, log_inst)
            sys.exit(1)
    return out

def subprocessT (args):
    p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()

def file_grep (pattern, file, mode='first'):
    if mode not in ['first','multiple']:
        raise NameError("Invalid option, "+mode)
    result_list = []
    with open(file,'r') as f:
        for line in f.readlines():
            if re.search(pattern,line) is not None:
                if mode == 'first':
                    return line.strip()
                elif mode == 'multiple':
                    result_list.append(line.strip())
    return result_list

def random_color():
    r = random.randrange(0,255)
    g = random.randrange(0,255)
    b = random.randrange(0,255)
    return [r,g,b]

def random_colors(n, minimum_distance=40):
    pal = []
    while len(pal) < n:
        newcolor = random_color()
        if any([abs(newcolor[0]-x[0])<minimum_distance for x in pal]) and any([abs(newcolor[1]-x[1])<minimum_distance for x in pal]) and any([abs(newcolor[2]-x[2])<minimum_distance for x in pal]):
            continue
        else:
            pal.append(newcolor)
    return pal

def alltaxtab2dict (alltaxtab): ## takes a per_line generator as arg (like file_yield_lines)
    db = {}
    for line in open(alltaxtab, 'r'):
        line = line.replace("'","")
        spl = [x.strip() for x in line.split('\t')]
        if len(spl[2]) == 0 or len(spl[8]) == 0:
            continue
        db[spl[2]] = spl[8]
    return db

def spdb_grab_os(spdb): ## takes a per_line generator as arg (like file_yield_lines)
    pattern = re.compile("^>.+{}".format(SPDB_OS_REGEXP_PATTERN))
    all_sp_os = []
    for line in open(spdb, 'r'):
        sp_os = None
        s = re.search(pattern, line)
        if s is not None:
            sp_os = s.group(1).strip()
            sp_os = sp_os.replace("'","")
            sp_os = sp_os.replace("#","")
            all_sp_os.append(sp_os)
    return list(set(all_sp_os))

# This function takes a 2-nested dictionary and flattens it into a new dictionary
# Primary keys are discarded
# Duplicated secondary keys will overwrite eachother
def flatten_dict(nested_dict):
    flattened_dict = {}
    for _,d in nested_dict.items():
        for k,v in d.items():
            flattened_dict[k] = v
    return flattened_dict
