import sys
import pandas as pd
import numpy as np
import operator
import logging
import threading
import re
from sequence import *
import os
import inspect
import cPickle as pickle
import subprocess
import settings
from time import localtime, strftime, sleep
from ete3 import Tree, TreeStyle, NodeStyle, NCBITaxa, TextFace ## to avoid seg fault on flux, change back once fixed
import random
#path_to_this_file = inspect.getfile(inspect.currentframe())
#this_file = os.path.split(path_to_this_file)[1]

### SOME VARIABLES ###
SPDB_OS_REGEXP_PATTERN = "OS=(.+?)(?:OX=.+|GN=.+|PE=.+|SV=.+|$)"
ERASE_LINE = '\x1b[2K'
CURSOR_UP_ONE = '\x1b[1A'
output_cols = {
        'GREEN':'\033[0;32m',
        'RED':'\033[0;31m',
        'RESET':'\033[0m'       
        }

### GENERAL FUNCTIONS FOR ALL MODULES ###
#%%
def get_run_opts ():
    opts_fname = "scgid.opts"
    opts = {}
    if not os.path.isfile(opts_fname):
        with open(opts_fname, 'w') as f:
            f.write("run_wd = "+os.getcwd()+"\n")        
    with open(opts_fname, 'r') as f:
        for o in f.readlines():
            spl = map(str.strip, o.split("="))
            opts[spl[0]] = spl[1]
        return opts

#%%
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

#%%
class color_changing_formatter(logging.Formatter):
    def __init__(self, fmt="%(levelno)s: %(msg)s",datefmt="%Y-%m-%d @ %H:%M:%S"):
        logging.Formatter.__init__(self, fmt, datefmt)

    def format(self, record):
        orig_format = self._fmt
        if record.levelno == 20:
            self._fmt = bcolors.CYAN+'[%(asctime)s]'+bcolors.ENDC+' %(message)s'
        elif record.levelno == 30:
            self._fmt = bcolors.WARNING+'[%(asctime)s]'+bcolors.ENDC+' %(message)s'
        elif record.levelno in [40,50]:
            self._fmt = bcolors.FAIL+'[%(asctime)s]'+bcolors.ENDC+' %(message)s'
        result = logging.Formatter.format(self, record)
        self._fmt = orig_format
        return result

#%%
def logger1(run_wd, module):
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    fh = logging.FileHandler(os.path.join(run_wd,'scgid.log'))
    fh.setLevel(logging.INFO)
    ch_formatter = color_changing_formatter()
    tolog_formatter = logging.Formatter('[%(asctime)s]['+module+'] %(message)s','%y-%m-%d @ %H:%M:%S')
    ch.setFormatter(ch_formatter)
    fh.setFormatter(tolog_formatter)
    logger.addHandler(ch)
    logger.addHandler(fh)
    
    return logger    
#%%
def start_logging(module, args):
    #from log import logger
    run_opts = get_run_opts()
    logger = logger1(run_opts['run_wd'], module)
    now = strftime("%B %d, %Y at %H:%M:%S", localtime())
    logger.info("XXXXXXXXXXXXXXXX NEW CALL TO <SCGID "+module.upper()+"> ON "+now.upper()+" XXXXXXXXXXXXXXXX")
    logger.info(' '.join(["scgid",module,' '.join(args[1:])]))
    return logger
                
#%%
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
#%%
def wrapper (func, ret):
    if type(ret) == dict:
        ret.update(func())
    elif type(ret) == list:
        ret.append(func())
    elif ret is None:
        func()

def do_long_wait(func, rt):
    ret_types = {
            'dict': dict(),
            'list': list(),
            'none': None
            }
    ret = ret_types[rt]
    t = threading.Thread(target=wrapper, args=(func, ret))
    t.start()
    spin = spinner()
    print str() #move cursor down to not erase previous line
    while t.is_alive():
        print next(spin)
        t.join(0.2)
    print CURSOR_UP_ONE,ERASE_LINE  #erase last spinner char
    if ret is None:
        return 0
    else:
     return ret

#%%
def pickle_loader (pklFile):
    try:
        while True:
            yield pickle.load(pklFile)
    except EOFError:
            pass
#%%
def pkl_fasta_in_out (org_fname, seq_type = "nucl", contig_info = True):
    #pkl_fname = org_fname.split("/")[-1]+'.pkl'
    pkl_fname = "{}.pkl".format(org_fname)
    #if os.path.isfile(os.getcwd()+"/"+pkl_fname):
    if os.path.isfile(pkl_fname) and os.stat(org_fname).st_ctime < os.stat(pkl_fname).st_ctime:
        #print "LOADING FROM PKL"
        objlist = []
        with open(pkl_fname,'r') as input:
            for seq in pickle_loader(input):
                objlist.append(seq)
    else:
        #print "MAKING NEW PKL"
        objlist = readFasta(org_fname, seq_type, contig_info)
        with open(pkl_fname,'w') as output:
            for seq in objlist:
                pickle.dump(seq, output, pickle.HIGHEST_PROTOCOL)
    return objlist


#%%
def ow_last_stdout(string):
    sys.stdout.write(CURSOR_UP_ONE)
    sys.stdout.write(ERASE_LINE)
    sys.stdout.write(string)
#%%
def report_outcome(msg, color, outcome):
    ow_last_stdout('\r'+msg + output_cols[color] + outcome + output_cols['RESET'] + '\n')

#%%
### FTP related functions ###

def ftp_retr_progress (block, lfile, tsize):
    ret = open(lfile,'ab').write(block)
    csize = os.path.getsize(lfile)
    prog = (csize*50)/(tsize)
    sys.stdout.write(CURSOR_UP_ONE)
    sys.stdout.write(ERASE_LINE)
    sys.stdout.write("["+"#"*(prog)+" "*(50-prog)+"] "+str(prog*2)+"%\n")
    return(ret)

def ftp_retr_and_report (ftp_inst, sname, rfile, lfile):
    tsize = ftp_inst.size(rfile)
    print "> Retrieving %s from %s ...\n[" % (rfile,sname) + " "*50 + "] 0%"
    ftp_inst.retrbinary("RETR %s" % (rfile), lambda block: ftp_retr_progress(block, lfile, tsize))

### gc_cov functions ###

#%%
def get_numeric_range(pd_series):
    range = (pd_series.min(),pd_series.max())
    return range

#%%
def calc_gc_window_1tailed (info_table, target_taxa, inc_factor = 0.01, plot = False):
    assert type(info_table) is pd.core.frame.DataFrame, "Input info_table to parse_infotable() must be a pandas dataframe."

    ### Get summary stats needed for determining gc window from info_table ###
    target = info_table.loc[info_table['parse_lineage'] == 'target']
    target_mean_gc = np.mean(target.gc)
    target_std_gc = np.std(target.gc)
    target_gc_range = (target.gc.min(),target.gc.max())
    all_target = float(target.shape[0])
    all_nontarget = float(info_table.shape[0] - all_target)

    ## Account for case where there is no nontarget (should be really rare), avoid division by zero error
    if all_nontarget == 0:
        all_nontarget = 1

    increment = float(target_std_gc)*inc_factor
    points = {
            "target": [],
            "nontarget": []
            }

    while (target_mean_gc - increment >= target_gc_range[0]) and (target_mean_gc + increment <= target_gc_range[1]):
        sel = (target_mean_gc - increment, target_mean_gc + increment)
        window = info_table[(info_table['gc'] >= sel[0]) & (info_table['gc'] <= sel[1])]
        target_in_window = float(window.loc[window['parse_lineage'] == 'target'].shape[0])
        nontarget_in_window = float(window.shape[0] - target_in_window)
        points["target"].append(target_in_window/all_target)
        points["nontarget"].append(nontarget_in_window/all_nontarget)
        increment += target_std_gc*inc_factor


    diff = map(operator.sub,points["target"],points["nontarget"])
    tradeoffs = map(operator.mul,points["target"],diff)
    points["tradeoffs"] = tradeoffs

    num_maxes = [i for i in points["tradeoffs"] if i == max(points["tradeoffs"])]
    if len(num_maxes) > 1:
        ## get max tradeoff at highest step (ie for asymtotic tradeoff curves)
        idx = 0
        indices_of_max = []
        for val in points["tradeoffs"]:
            if val == max(points["tradeoffs"]):
                indices_of_max.append(idx)
            idx+=1
        best_step_idx = max(indices_of_max)
    else:
        best_step_idx = points["tradeoffs"].index(max(points["tradeoffs"]))

    opt_tail = target_std_gc * inc_factor * (best_step_idx+1)
    final_window = (target_mean_gc - opt_tail, target_mean_gc + opt_tail)

    if plot is True:
        print "Final window:",final_window[0],"-->",final_window[1]
        return points

    else:
        return final_window
#%%
def calc_gc_window_2tailed (info_table, target_taxa, inc_factor = 0.01, plot = False):
    target = info_table.loc[info_table['parse_lineage'] == 'target']
    target_mean_gc = np.mean(target.gc)
    target_std_gc = np.std(target.gc)
    target_gc_range = (target.gc.min(),target.gc.max())
    all_target = float(target.shape[0])
    all_nontarget = float(info_table.shape[0] - all_target)

    ## Account for case where there is no nontarget (should be really rare), avoid division by zero error
    if all_nontarget == 0:
        all_nontarget = 1

    increment = float(target_std_gc)*inc_factor
    points = {
            "neg":{
                    "target": [],
                    "nontarget": []
                    },
            "pos":{
                    "target": [],
                    "nontarget": []
                    }
            }
    while target_mean_gc - increment >= target_gc_range[0]:
        sel = (target_mean_gc - increment, target_mean_gc)
        window = info_table[(info_table['gc'] >= sel[0]) & (info_table['gc'] <= sel[1])]
        target_in_window = float(window.loc[window['parse_lineage'] == 'target'].shape[0])
        nontarget_in_window = float(window.shape[0] - target_in_window)
        points["neg"]["target"].append(target_in_window/all_target)
        points["neg"]["nontarget"].append(nontarget_in_window/all_nontarget)
        increment += target_std_gc*inc_factor

    ### Reset increment ###
    increment = float(target_std_gc)*inc_factor

    while target_mean_gc + increment <= target_gc_range[1]:
        sel = (target_mean_gc, target_mean_gc + increment)
        window = info_table[(info_table['gc'] >= sel[0]) & (info_table['gc'] <= sel[1])]
        target_in_window = float(window.loc[window['parse_lineage'] == 'target'].shape[0])
        nontarget_in_window = float(window.shape[0] - target_in_window)
        points["pos"]["target"].append(target_in_window/all_target)
        points["pos"]["nontarget"].append(nontarget_in_window/all_nontarget)
        increment += target_std_gc*inc_factor


    diff_neg = map(operator.sub,points["neg"]["target"],points["neg"]["nontarget"])
    diff_pos = map(operator.sub,points["pos"]["target"],points["pos"]["nontarget"])

    tradeoffs_neg = map(operator.mul,points["neg"]["target"],diff_neg)
    tradeoffs_pos = map(operator.mul,points["pos"]["target"],diff_pos)

    points["neg"]["tradeoffs"] = tradeoffs_neg
    points["pos"]["tradeoffs"] = tradeoffs_pos

    num_maxes_neg = [i for i in points["neg"]["tradeoffs"] if i == max(points["neg"]["tradeoffs"])]
    num_maxes_pos = [i for i in points["pos"]["tradeoffs"] if i == max(points["pos"]["tradeoffs"])]

    ## make sure that the following FOR loop is necessary in the negative direction (ie asymptotic tradeoff curves)

    if len(num_maxes_neg) > 1:
        ## get max tradeoff at highest step in NEGATIVE direction (ie for asymtotic tradeoff curves)
        idx = 0
        indices_of_max_neg = []
        for val in points["neg"]["tradeoffs"]:
            if val == max(points["neg"]["tradeoffs"]):
                indices_of_max_neg.append(idx)
            idx+=1
        best_step_idx_neg = max(indices_of_max_neg)
    else:
        best_step_idx_neg = points["neg"]["tradeoffs"].index(max(points["neg"]["tradeoffs"]))

    if len(num_maxes_pos) > 1:
        ## get max tradeoff at highest step in POSITIVE direction (ie for asymtotic tradeoff curves)
        idx = 0
        indices_of_max_pos = []
        for val in points["pos"]["tradeoffs"]:
            if val == max(points["pos"]["tradeoffs"]):
                indices_of_max_pos.append(idx)
            idx+=1
        best_step_idx_pos = max(indices_of_max_pos)
    else:
        best_step_idx_pos = points["pos"]["tradeoffs"].index(max(points["pos"]["tradeoffs"]))


    opt_tail_neg = target_std_gc * inc_factor * (best_step_idx_neg+1)
    opt_tail_pos = target_std_gc * inc_factor * (best_step_idx_pos+1)

    opt_window_neg = target_mean_gc - opt_tail_neg
    opt_window_pos = target_mean_gc + opt_tail_pos

    final_window = (opt_window_neg, opt_window_pos)

    if plot is True:
        print "Final window:",opt_window_neg,"-->",opt_window_pos
        return points

    else:
        return final_window
#%%

### When should there NOT be a coverage cut-off and how do I code that? ###

def calc_coverage_window_1tailed (info_table, target_taxa, inc_factor = 0.01, plot = False):
    target = target = info_table.loc[info_table['parse_lineage'] == 'target']
    target_mean_cov = np.mean(target.coverage)
    target_std_cov = np.std(target.coverage)
    target_cov_range = (target.coverage.min(),target.coverage.max())
    all_target = float(target.shape[0])
    all_nontarget = float(info_table.shape[0] - all_target)
    increment = float(target_std_cov) * inc_factor

    ## Account for case where there is no nontarget (should be really rare), avoid division by zero error
    if all_nontarget == 0:
        all_nontarget = 1

    points = {
            "target": [],
            "nontarget": []
            }
    while (target_mean_cov - increment >= target_cov_range[0]) and (target_mean_cov + increment <= target_cov_range[1]):
        sel = (target_mean_cov - increment, target_mean_cov + increment)
        window = info_table[(info_table['coverage'] >= sel[0]) & (info_table['coverage'] <= sel[1])]
        target_in_window = float(window.loc[window['parse_lineage'] == 'target'].shape[0])
        nontarget_in_window = float(window.shape[0] - target_in_window)
        points["target"].append(target_in_window/all_target)
        points["nontarget"].append(nontarget_in_window/all_nontarget)
        increment += target_std_cov*inc_factor

    diff = map(operator.sub,points["target"],points["nontarget"])
    tradeoffs = map(operator.mul,points["target"],diff)
    points["tradeoffs"] = tradeoffs

    num_maxes = [i for i in points["tradeoffs"] if i == max(points["tradeoffs"])]

    if len(num_maxes) > 1:
        ## get max tradeoff at highest step (ie for asymtotic tradeoff curves)
        idx = 0
        indices_of_max = []
        for val in points["tradeoffs"]:
            if val == max(points["tradeoffs"]):
                indices_of_max.append(idx)
            idx+=1
        best_step_idx = max(indices_of_max)
    else:
        best_step_idx = points["tradeoffs"].index(max(points["tradeoffs"]))

    opt_tail = target_std_cov * inc_factor * (best_step_idx+1)
    final_window = (target_mean_cov - opt_tail, target_mean_cov + opt_tail)

    if plot is True:
        print "Final window:",final_window[0],"-->",final_window[1]
        return points

    else:
        return final_window

#%%
def calc_coverage_window_2tailed (info_table, target_taxa, inc_factor = 0.01, plot = False):

    target = info_table.loc[info_table['parse_lineage'] == 'target']
    target_mean_cov = np.mean(target.coverage)
    target_std_cov = np.std(target.coverage)
    target_cov_range = (target.coverage.min(),target.coverage.max())
    increment = float(target_std_cov)*inc_factor
    all_target = float(target.shape[0])
    all_nontarget = float(info_table.shape[0] - all_target)

    ## Account for case where there is no nontarget (should be really rare), avoid division by zero error
    if all_nontarget == 0:
        all_nontarget = 1

    points = {
            "neg":{
                    "target": [],
                    "nontarget": []
                    },
            "pos":{
                    "target": [],
                    "nontarget": []
                    }
            }

    while target_mean_cov - increment >= target_cov_range[0]:
        sel = (target_mean_cov - increment, target_mean_cov)
        window = info_table[(info_table['coverage'] >= sel[0]) & (info_table['coverage'] <= sel[1])]
        target_in_window = float(window.loc[window['parse_lineage'] == 'target'].shape[0])
        nontarget_in_window = float(window.shape[0] - target_in_window)
        points["neg"]["target"].append(target_in_window/all_target)
        points["neg"]["nontarget"].append(nontarget_in_window/all_nontarget)
        increment += target_std_cov*inc_factor

    ### Reset increment ###
    increment = float(target_std_cov)*inc_factor

    while target_mean_cov + increment <= target_cov_range[1]:
        sel = (target_mean_cov, target_mean_cov + increment)
        window = info_table[(info_table['coverage'] >= sel[0]) & (info_table['coverage'] <= sel[1])]
        target_in_window = float(window.loc[window['parse_lineage'] == 'target'].shape[0])
        nontarget_in_window = float(window.shape[0] - target_in_window)
        ## Account for case where there is no nontarget (should be really rare), avoid division by zero error
        if all_nontarget == 0:
            all_nontarget = 1

        points["pos"]["target"].append(target_in_window/all_target)
        points["pos"]["nontarget"].append(nontarget_in_window/all_nontarget)


        increment += target_std_cov*inc_factor


    diff_neg = map(operator.sub,points["neg"]["target"],points["neg"]["nontarget"])
    diff_pos = map(operator.sub,points["pos"]["target"],points["pos"]["nontarget"])
    #print len(points["pos"]["target"])
    #print len(diff_pos)

    tradeoffs_neg = map(operator.mul,points["neg"]["target"],diff_neg)
    tradeoffs_pos = map(operator.mul,points["pos"]["target"],diff_pos)

    points["neg"]["tradeoffs"] = tradeoffs_neg
    points["pos"]["tradeoffs"] = tradeoffs_pos

    ## how many points represent the maximum trade-off?
    num_maxes_neg = [i for i in points["neg"]["tradeoffs"] if i == max(points["neg"]["tradeoffs"])]
    num_maxes_pos = [i for i in points["pos"]["tradeoffs"] if i == max(points["pos"]["tradeoffs"])]

    ## make sure that the following FOR loop is necessary in the negative direction (ie asymptotic tradeoff curves)

    if len(num_maxes_neg) > 1:
        ## get max tradeoff at highest step in NEGATIVE direction (ie for asymptotic tradeoff curves)
        idx = 0
        indices_of_max_neg = []
        for val in points["neg"]["tradeoffs"]:
            if val == max(points["neg"]["tradeoffs"]):
                indices_of_max_neg.append(idx)
            idx+=1
        best_step_idx_neg = max(indices_of_max_neg)

    else:
        best_step_idx_neg = points["neg"]["tradeoffs"].index(max(points["neg"]["tradeoffs"]))

    ## get max tradeoff at highest step in POSITIVE direction (ie for asymtotic tradeoff curves)
    if len(num_maxes_pos) > 1:
        ## get max tradeoff at highest step in POSITIVE direction (ie for asymptotic tradeoff curves)
        idx = 0
        indices_of_max_pos = []
        for val in points["pos"]["tradeoffs"]:
            if val == max(points["pos"]["tradeoffs"]):
                indices_of_max_pos.append(idx)
            idx+=1
        best_step_idx_pos = max(indices_of_max_pos)

    else:
        best_step_idx_pos = points["pos"]["tradeoffs"].index(max(points["pos"]["tradeoffs"]))
    #print best_step_idx_pos

    #best_step_idx_neg = points["neg"]["tradeoffs"].index(max(points["neg"]["tradeoffs"]))
    #best_step_idx_pos = points["pos"]["tradeoffs"].index(max(points["pos"]["tradeoffs"]))

    opt_tail_neg = target_std_cov * inc_factor * (best_step_idx_neg+1)
    opt_tail_pos = target_std_cov * inc_factor * (best_step_idx_pos+1)

    #print opt_tail_pos

    opt_window_neg = target_mean_cov - opt_tail_neg
    opt_window_pos = target_mean_cov + opt_tail_pos

    #print opt_window_pos

    final_window = (opt_window_neg, opt_window_pos)

    if plot is True:
        print "Final window:",opt_window_neg,"-->",opt_window_pos
        return points

    else:
        return final_window
#%%
def get_window_table (info_table, window, param):
    assert type(info_table) is pd.core.frame.DataFrame, "Input info_table to parse_infotable() must be a pandas dataframe."
    assert type(window) is tuple, "Arg window must be a tuple."
    assert param in info_table.columns, "Specified parameter is not present in input info_table."

    #print window[0],window[1]
    window_table = info_table[(info_table[param] >= window[0]) & (info_table[param] <= window[1])]
    return window_table

#%%

def generate_windows(info_table, target_taxa, target):

    windows = []

    gc1noCov = {}
    gc2noCov = {}
    cov1noGc = {}
    cov2noGc = {}

    gc1cov1 = {}
    gc1cov2 = {}
    gc2cov1 = {}
    gc2cov2 = {}

    cov1gc1 = {}
    cov1gc2 = {}
    cov2gc1 = {}
    cov2gc2 = {}

    gc1 = calc_gc_window_1tailed (info_table,target_taxa,inc_factor=0.01,plot=False)
    gc2 = calc_gc_window_2tailed (info_table,target_taxa,inc_factor=0.01,plot=False)
    cov1 = calc_coverage_window_1tailed(info_table,target_taxa,inc_factor=0.001,plot=False)
    cov2 = calc_coverage_window_2tailed(info_table,target_taxa,inc_factor=0.001,plot=False)
    #status.task("(4/12)...")

    gc1cov1['gc'] = gc1
    gc1cov2['gc'] = gc1
    gc2cov1['gc'] = gc2
    gc2cov2['gc'] = gc2

    gc1noCov['gc'] = gc1
    gc1noCov['cov'] = get_numeric_range(target.coverage)
    gc2noCov['gc'] = gc2
    gc2noCov['cov'] = get_numeric_range(target.coverage)

    cov1noGc['gc'] = get_numeric_range(target.gc)
    cov1noGc['cov'] = cov1
    cov2noGc['gc'] = get_numeric_range(target.gc)
    cov2noGc['cov'] = cov2

    cov1gc1['cov'] = cov1
    cov1gc2['cov'] = cov1
    cov2gc1['cov'] = cov2
    cov2gc2['cov'] = cov2

    gc1 = get_window_table(info_table,gc1,'gc')
    gc2 = get_window_table(info_table,gc2,'gc')
    cov1 = get_window_table(info_table,cov1,'coverage')
    cov2 = get_window_table(info_table,cov2,'coverage')

    gc1cov1['cov'] = calc_coverage_window_1tailed(gc1,target_taxa,inc_factor=0.01,plot=False)
    gc1cov2['cov'] = calc_coverage_window_2tailed(gc1,target_taxa,inc_factor=0.01,plot=False)

    gc2cov1['cov'] = calc_coverage_window_1tailed(gc2,target_taxa,inc_factor=0.01,plot=False)
    gc2cov2['cov'] = calc_coverage_window_2tailed(gc2,target_taxa,inc_factor=0.01,plot=False)

    #status.task("(8/12)...")

    cov1gc1['gc'] = calc_gc_window_1tailed(cov1,target_taxa,inc_factor=0.01,plot=False)
    cov1gc2['gc'] = calc_gc_window_2tailed(cov1,target_taxa,inc_factor=0.01,plot=False)

    cov2gc1['gc'] = calc_gc_window_1tailed(cov2,target_taxa,inc_factor=0.01,plot=False)
    cov2gc2['gc'] = calc_gc_window_2tailed(cov2,target_taxa,inc_factor=0.01,plot=False)

    #status.success("(12/12)...")

    windows.append(gc1cov1)
    windows.append(gc1cov2)
    windows.append(gc2cov1)
    windows.append(gc2cov2)
    windows.append(cov1gc1)
    windows.append(cov1gc2)
    windows.append(cov2gc1)
    windows.append(cov2gc2)
    windows.append(gc1noCov)
    windows.append(gc2noCov)
    windows.append(cov1noGc)
    windows.append(cov2noGc)

    return windows
#%%

class gc_cov_window ():
    def __init__ (self, lower_gc=0.00, upper_gc=0.00, lower_cov=0.00, upper_cov=0.00,
                  tailedness="", tp=0.00, ntp=0.00):
        self.lower_gc = lower_gc
        self.upper_gc = upper_gc
        self.lower_cov = lower_cov
        self.upper_cov = upper_cov
        self.tailedness = tailedness
        self.tp = tp
        self.ntp = ntp
#%%
def node_num_only (name):
    spl = name.split('_')
    return '_'.join(spl[0:2])

#%%
def gff3_to_fasta(gff3, outname):
    prots = []
    recording = False
    for line in open(gff3).readlines():
        if recording is False:
            s = re.search("[#][ ]protein[ ]sequence[ ][=][ ]",line)
            if s is not None:
                recording = True
                s = re.search("[#][ ]protein[ ]sequence[ ][=][ ][[]([A-Za-z]+)",line)
                protein_sequence = s.group(1)
        elif recording is True:
            s = re.search("[#][ ]end[ ]gene[ ]([A-Za-z0-9_.]+)",line)
            if s is not None:
                recording = False
                prots.append(sequence(s.group(1), protein_sequence, seq_type = "prot", contig_info = False))
                protein_sequence = ""
            else:
                s = re.search("[#][ ]([a-zA-Z]+)",line)
                protein_sequence += s.group(1)
    with open(outname,'w') as fasta:
        for i in prots:
            fasta.write(i.outFasta()+"\n")
#%%
### RSCU functions ###

#%%
def get_blank_codon_counts ():
    blank = {
        'UUU': 0, 'UUC': 0, 'UUA': 0, 'UUG': 0,
        'CUU': 0, 'CUC': 0, 'CUA': 0, 'CUG': 0,
        'AUU': 0, 'AUC': 0, 'AUA': 0, 'AUG': 0,
        'GUU': 0, 'GUC': 0, 'GUA': 0, 'GUG': 0,
        'UCU': 0, 'UCC': 0, 'UCA': 0, 'UCG': 0,
        'AGU': 0, 'AGC': 0, 'CCU': 0, 'CCC': 0,
        'CCA': 0, 'CCG': 0, 'ACU': 0, 'ACC': 0,
        'ACA': 0, 'ACG': 0, 'GCU': 0, 'GCC': 0,
        'GCA': 0, 'GCG': 0, 'UAU': 0, 'UAC': 0,
        'UAA': 0, 'UAG': 0, 'UGA': 0, 'CAU': 0,
        'CAC': 0, 'CAA': 0, 'CAG': 0, 'AAU': 0,
        'AAC': 0, 'AAA': 0, 'AAG': 0, 'GAU': 0,
        'GAC': 0, 'GAA': 0, 'GAG': 0, 'UGU': 0,
        'UGC': 0, 'UGG': 0, 'CGU': 0, 'CGC': 0,
        'CGA': 0, 'CGG': 0, 'AGA': 0, 'AGG': 0,
        'GGU': 0, 'GGC': 0, 'GGA': 0, 'GGG': 0
        }
    return dict(blank)
#%%
def get_codon_table ():
    codon_table = {
        'Phe': ['UUU','UUC'],
        'Leu': ['UUA','UUG','CUU','CUC','CUA','CUG'],
        'Ile': ['AUU','AUC','AUA'],
        'Met': ['AUG'],
        'Val': ['GUU','GUC','GUA','GUG'],
        'Ser': ['UCU','UCC','UCA','UCG','AGU','AGC'],
        'Pro': ['CCU','CCC','CCA','CCG'],
        'Thr': ['ACU','ACC','ACA','ACG'],
        'Ala': ['GCU','GCC','GCA','GCG'],
        'Tyr': ['UAU','UAC'],
        'STOP': ['UAA','UAG','UGA'],
        'His': ['CAU','CAC'],
        'Gln': ['CAA','CAG'],
        'Asn': ['AAU','AAC'],
        'Lys': ['AAA','AAG'],
        'Asp': ['GAU','GAC'],
        'Glu': ['GAA','GAG'],
        'Cys': ['UGU','UGC'],
        'Trp': ['UGG'],
        'Arg': ['CGU','CGC','CGA','CGG','AGA','AGG'],
        'Gly': ['GGU','GGC','GGA','GGG']
        }
    return dict(codon_table)
#%%
def transcribe(nucl_seq):
    transcript = ""
    trans = {
            'A':'U',
            'T':'A',
            'G':'C',
            'C':'G',
            'N':'N'
            }
    for letter in nucl_seq:
        transcript+=trans[letter]
    return transcript
#%%
def complement(string):
    out = ""
    comp = {
            'A':'U',
            'U':'A',
            'G':'C',
            'C':'G',
            'N':'N'
            }
    for letter in string:
        out+=comp[letter]
    return out

#%%
def translate_dna_to_protein (string):
    dna_to_aa = {'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
    out = ""
    for pos in range(0,len(string),3):
        out+= dna_to_aa[string[pos:pos+3]]
    return out

#%%
def best_blast_hit (tabular, bitcol=7):
    best = {}
    for line in open(tabular).readlines():
        spl = line.split("\t")
        spl = map(str.strip,spl)
        label = spl[0]
        bit = float(spl[bitcol])
        if label in best.keys():
            if best[label][bitcol] < bit:
                best[label] = spl
        else:
            best[label] = spl

    return best
#%%
def extract_cds_gff3 (gff3, nucl):
    nucl_dict = {}
    for obj in nucl:
        #''' shortnames
        node = obj.label.split('_')[0:2]
        node = '_'.join(node)
        #node = node.replace("NODE_","N")
        #'''
        #node = obj.label
        nucl_dict[node] = obj.sequence
    cds_only = {}
    for line in open(gff3).readlines():
        if line[0] == "#":
            continue
        if line.split('\t')[2] == "CDS":
            spl = line.split("\t")
            #''' shortnames
            node = spl[0].split('_')[0:2]
            node = '_'.join(node)
            #node = node.replace("NODE_","N")
            #'''
            #node = spl[0]
            s = re.search("[.](g[0-9]+)[.]",spl[8])
            pid = s.group(1)
            if node in cds_only.keys():
                if pid in cds_only[node].keys():
                    cds_only[node][pid].append(spl[1:len(spl)-1])
                else:
                    cds_only[node][pid] = [spl[1:len(spl)-1]]
            else:
                cds_only[node] = {
                        pid: [spl[1:len(spl)-1]]
                        }
    cds_cat = []
    for node,pids in cds_only.iteritems():
        node_cds_cat = ""
        #print node
        for prot,cds_chunks in pids.iteritems():
            gene_cds = ""
            for cds in cds_chunks:
                start = cds[2]
                stop = cds[3]
                gene_cds += nucl_dict[node][int(start)-1:int(stop)]
            if len(gene_cds) % 3 == 0:
                node_cds_cat += gene_cds
        cds_cat.append(sequence(node, node_cds_cat, seq_type='nucl', contig_info=False))
    return cds_cat
#%%
def remove_small_sequences (list_of_sequence_objects, min_len):
    list_to_output = []
    for sequence_object in list_of_sequence_objects:
        if len(sequence_object.sequence) >= int(min_len):
            list_to_output.append(sequence_object)
    return list_to_output

#%%
def annotate_tips_prot (tree, target_taxa, info_table):
    ## tip styling by annotation ##
    tar = NodeStyle()
    tar["vt_line_color"] = "blue"
    tar["hz_line_color"] = "blue"
    tar["hz_line_width"] = 25
    tar["shape"] = "circle"
    
    ntar = NodeStyle()
    ntar["vt_line_color"] = "red"
    ntar["hz_line_color"] = "red"
    ntar["hz_line_width"] = 25
    ntar["shape"] = "circle"
    
    no_class = NodeStyle()
    no_class["vt_line_color"] = "black"
    no_class["hz_line_color"] = "black"
    no_class["shape"] = "circle"
    no_class["hz_line_width"] = 25
    
    info_table.decide_taxonomy()
    '''
    #Reformat info table for tip annotation
    grouped = info_table.groupby('contig')
    
    reformed = grouped.agg({'pid': lambda x: ','.join(x),
                            'evalue': lambda x: ','.join(str(e) for e in x),
                            'parse_lineage': lambda x: ','.join(x)})
    
    reformed.parse_lineage = reformed.parse_lineage.apply(str.split,args=',')
    reformed.pid = reformed.pid.apply(str.split,args=',')
    reformed.evalue = reformed.evalue.apply(str.split,args=',')
    reformed = reformed.reset_index()
    
    to_add = []
    to_exclude = []
    
    for i in reformed.itertuples():
        #i[0] = index, i[1] = contig, i[2] = pid, i[3] = evalue, i[4] parsed_lineage_lists
        num_t = i.parse_lineage.count('target')
        num_nt = i.parse_lineage.count('nontarget')
        shortname = i.contig.split("_")[0:2]
        shortname = shortname[0][0]+shortname[1]
        if num_t > num_nt:
            to_add.append(shortname)
        elif num_t == num_nt:
            evalues = [float(e) for e in i.evalue]
            maxes = [x for x in evalues if x == min(evalues)] ## best evalue means lowest, ie min()
            #print maxes
            if len(maxes) == 1:
                #print i.evalue
                best_idx = evalues.index(min(evalues))
                if i.parse_lineage[best_idx] == 'target':
                    to_add.append(shortname)
                else:
                    to_exclude.append(shortname)
        else:
            to_exclude.append(shortname)
    '''
    ## Apply the variable formatting to annotated tips ##
    for n in tree.traverse():
        if n.is_leaf():
            n.add_feature("annotation","unclassified")
            if n.name in info_table.keep:
                n.annotation = "target"
                n.set_style(tar)
                n.add_face(TextFace(n.name,fsize=500),column=0)
            elif n.name in info_table.dump:
                n.annotation = "nontarget"
                n.set_style(ntar)
                
            else:
                n.set_style(no_class)
    #sys.exit()
    return tree
            
    
    
def annotate_tips(tree, target_taxa, path_to_tids):

    ## tip styling by annotation ##
    tar = NodeStyle()
    tar["shape"] = "circle"
    tar["size"] = 35
    tar["fgcolor"] = "blue"

    ntar = NodeStyle()
    ntar["shape"] = "circle"
    ntar["size"] = 35
    ntar["fgcolor"] = "red"

    no_class = NodeStyle()
    no_class["shape"] = 'circle'
    no_class["size"] = 35
    no_class["fgcolor"] = "black"

    for l in tree.get_leaves():
        l.add_feature("annotation","unclassified")

    ## get lineage information from ncbi taxdb ##
    #to_log = ""
    #sys.stdout = to_log
    ncbi = NCBITaxa()
    ids = {}
    for i in open(path_to_tids).readlines():
        cols = i.split('\t')
        name = '_'.join(cols[0].split('_')[0:2])
        tids = cols[1].strip()
        if len(tids.split(';')) > 1:
            ids[name] = {}
            for tid in tids.split(';'):
                num = int(tid)
                if name in ids.keys():
                    ids[name].update(ncbi.get_taxid_translator(ncbi.get_lineage(num)))
                else:
                    ids[name] = ncbi.get_taxid_translator(ncbi.get_lineage(num))
        else:
            num = int(tids)
            ids[name] = ncbi.get_taxid_translator(ncbi.get_lineage(num))

    for node,l in ids.iteritems():
        formatted = []
        for tid, taxon in l.iteritems():
            formatted.append(taxon)
        with open("ncbi_taxonomy",'a') as tax:
            tax.write(node+"\t")
            tax.write(';'.join(formatted)+'\n')
        ## Annotate the tip that got a hit ##
        if any(i in target_taxa['target'] for i in formatted) is True and any(i in target_taxa['exception'] for i in formatted) is False:
            for this_one in tree.iter_search_nodes(name=node):
                this_one.annotation = 'target'
        else:
            for this_one in tree.iter_search_nodes(name=node):
                this_one.annotation = 'nontarget'
    #logger.warning(to_log)
    #sys.stdout = sys.__stdout__

    ## Apply the variable formatting to annotated tips ##
    for n in tree.traverse():
        if n.is_leaf():
            if n.annotation is "target":
                n.set_style(tar)
            elif n.annotation is "nontarget":
                n.set_style(ntar)
            elif n.annotation is "unclassified":
                n.set_style(no_class)

    return tree

#%%
def specify_Xmx_esom(path_to_esomstart, mem):

    with open(path_to_esomstart,'r') as f:
        pos = 0
        trunc_point = 0
        for line in f.readlines():
            if line.split(' ')[0] == 'java':
                trunc_point = pos
                break
            pos += len(line)

    with open(path_to_esomstart,'ar') as f:
        f.seek(trunc_point)
        f.truncate()
        new_line = "java -cp \"$CP\" -Xmx"+mem+" $MAIN \"$@\""
        f.write(new_line)
#%%
def update_ESOM_HOME(path_to_esomstart, new_home):
    with open(path_to_esomstart, 'r') as f:
        match = 'ESOM_HOME[=]["]'
        reprint = []
        for line in f.readlines():
            if re.search(match,line) is not None:
                current_path = re.search('[\"]([^\"]+)[\"]',line).group(1)
                new_home_line = line.replace(current_path,new_home)
                reprint.append(new_home_line)
            else:
                reprint.append(line)
    with open(path_to_esomstart,'w') as f:
        f.write(''.join(reprint))
#%%
def file_grep (pattern, file, mode='first'):
    if mode not in ['first','multiple']:
        raise NameError("Invalid option, "+mode)
    result_list = []
    with open(file,'r') as f:
        for line in f.readlines():
            if re.search(pattern,line) is not None:
                if mode is 'first':
                    return line.strip()
                elif mode is 'multiple':
                    result_list.append(line.strip())
    return result_list


#%%
def replace_line(file, old, new):
    head = ""
    tail = ""
    line_start = 0
    line_end = 0
    with open(file,'r') as f:
        pos = 0
        trunc_point = 0
        present = False
        for line in f.readlines():
            if line.strip() == old:
                line_start = pos
                line_end = pos + len(line)
                present = True
                break
            pos += len(line)
        if present is False:
            raise ValueError("The line you want to replace is not present in the file.")
        f.seek(0)
        head = f.read(line_start)
        f.seek(line_end)
        tail = f.read()
    with open(file, 'w') as f:
        f.write(head)
        f.write(new+'\n')
        f.write(tail)

#%%
def replace_line_by_pattern (file, pattern, new):
    head = ""
    tail = ""
    line_start = 0
    line_end = 0
    with open(file,'r') as f:
        pos = 0
        trunc_point = 0
        present = False
        for line in f.readlines():
            if re.match(pattern, line) is not None:
                line_start = pos
                line_end = pos + len(line)
                present = True
                break
            pos += len(line)
        if present is False:
            raise ValueError("The line you want to replace is not present in the file.")
        f.seek(0)
        head = f.read(line_start)
        f.seek(line_end)
        tail = f.read()
    with open(file, 'w') as f:
        f.write(head)
        f.write(new+'\n')
        f.write(tail)

#%%
def get_line(file, to_get):
    with open(file,'r') as f:
        for line in f.readlines():
            if line.strip() == to_get:
                return to_get
#%%
def log_command_line_errors (err, log_inst):
    if err != '':
        err = err.split('\n')
        for i in err:
            if len(i) > 0:
                log_inst.critical(i)

#%%
def subprocessC (args):
    p = subprocess.call(args)
    if p != 0:
        sys.exit(1)

#%%
def subprocessP (args, log_inst):
    p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out,err = p.communicate() 
    if p.returncode != 0:
        log_command_line_errors(err, log_inst)
        sys.exit(1)
    return out

#%%
def subprocessT (args):
    p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()

#%%

def is_fasta(f, strict=False, verbose = False):
    lines = open(f).readlines()
    lines = map(str.strip,lines)

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
                print "Deviation from fasta format on line #"+str(idx+1)
        idx+=1
    if not fail:
        return True
    else:
        return False

#%%
def random_color():
    r = random.randrange(0,255,40)
    g = random.randrange(0,255,40)
    b = random.randrange(0,255,40)
    return [r,g,b]

#%%
class go_through_list:
    def __init__(self, l):
        self.current = 0
        self.end = len(l)
        self.thelist = l

    def __iter__ (self):
        return self

    def next(self):
        if self.current == self.end:
            raise StopIteration
        else:
            self.current+=1
            return self.thelist[self.current-1]

#%%
def read_lines_from_large_file (large_file, n=500):
    last_pos = -1
    p = 0
    stop = False
    with open(large_file,'r') as f:
        while True:
            for i in range(p,n):
                next_line = f.readline().strip()
                if f.tell() == last_pos:
                    stop = True
                    break
                last_pos = f.tell()
                yield next_line
                p += 1
            if stop:
                break
            n += n
#%%

def file_yield_lines(large_file):
    with open(large_file, 'r') as lf:
        for line in lf:
            yield line
            
def spdb_grab_os(spdb_line_gen): ## takes a per_line generator as arg (like file_yield_lines)
    pattern = re.compile("^>.+{}".format(SPDB_OS_REGEXP_PATTERN))
    all_sp_os = []
    for line in spdb_line_gen:
        sp_os = None
        s = re.search(pattern, line)
        if s is not None:
            sp_os = s.group(1).strip()
            sp_os = sp_os.replace("'","")
            sp_os = sp_os.replace("#","")
            all_sp_os.append(sp_os)
    return list(set(all_sp_os))

def alltaxtab2dict (alltaxtab_line_gen): ## takes a per_line generator as arg (like file_yield_lines)
    db = {}
    for line in alltaxtab_line_gen:
        line = line.replace("'","")
        spl = line.split('\t')
        spl = map(str.strip, spl)
        if len(spl[2]) == 0 or len(spl[8]) == 0:
            continue
        db[spl[2]] = spl[8]
    return db
#%%
def proc_sam_output_unaligned (lines, pe1, pe2, orph):
    for line in lines:
        spl = line.split('\t')
        if len(spl) < 4:
            continue
        else:
            if spl[1] == '77':
                with open(pe1,'a') as out:
                    out.write('@'+spl[0]+'\n'+spl[9]+'\n\n'+spl[10]+'\n')
            elif spl[1] == '141':
                with open(pe2,'a') as out:
                    out.write('@'+spl[0]+'\n'+spl[9]+'\n\n'+spl[10]+'\n')
            elif spl[1] == '4':
                with open(orph,'a') as out:
                    out.write('@'+spl[0]+'\n'+spl[9]+'\n'+spl[10]+'\n')
            else:
                continue
#%%
def parse_spdb_blastout(sp_fasta, blastout, log_inst=None):
    path_to_spdb = os.path.abspath(sp_fasta)
    sp_fasta = pkl_fasta_in_out(sp_fasta,seq_type="prot",contig_info=False)
    ids = {}

    ldict = []
    for entry in sp_fasta:
        spl = entry.label.split(" ",1)
        newrow = {
            'accession': spl[0],
            'description': spl[1]
            }
        ldict.append(newrow)
    sp_fasta = None #free

    lib = pd.DataFrame(ldict)
    lib = lib.set_index('accession')
    
    output = []
    with open(blastout, 'r') as b:
        for line in b:
            spl = line.split("\t")
            acc = line.split("\t")[1]
            contig = line.split("\t")[0]
            evalue = line.split("\t")[10]
        
            try:
                hit = lib.loc[acc].description
                output.append("{}\t{}\t{}\t{}".format(contig, acc, hit, evalue))
            except:
                if log_inst is not None:
                    log_inst.critical("Accession ({}) not found in swissprot-style database, {}. It is likely that you specified a different version of that database than that used to BLAST against.".format(acc, path_to_spdb))
                    raise ValueError
                else:
                    print "ERROR"
        
        return output
 

