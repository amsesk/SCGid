#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  8 09:28:00 2018

@author: aimzez
"""

### This library is meant to serve the initiate.py script, prior to install of some package dependencies

#%% imports
import sys
import os
import ftplib
import re

#%% some variables
ERASE_LINE = '\x1b[2K'
CURSOR_UP_ONE = '\x1b[1A'
output_cols = {
        'GREEN':'\033[0;32m',
        'RED':'\033[0;31m',
        'RESET':'\033[0m'       
        }
#%%
def ow_last_stdout(string):
    sys.stdout.write(CURSOR_UP_ONE)
    sys.stdout.write(ERASE_LINE)
    sys.stdout.write(string)
    
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
def report_outcome(msg, color, outcome):
    sys.stdout.write(CURSOR_UP_ONE)
    #sys.stdout.flush()
    sys.stdout.write(ERASE_LINE)
    sys.stdout.write('\r'+msg + output_cols[color] + outcome + output_cols['RESET'] + '\n')
    
#%% FTP-related functions

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