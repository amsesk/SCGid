#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 19:58:39 2017

@author: kevinamses
"""
import sys
import subprocess
import os
import inspect
from lib import file_grep, update_ESOM_HOME
import settings
import re

esom_help = """
----- SCGID KMERS -----

Usage: scgid kmers <task> [args...] --> try scgid kmers <task> -h|--help to see task-specific arguments

Please specify a task to continue:

    train       Train the ESOM topology
    annotate    Annotate the trained ESOM toplogy with taxonomic information
    extract     Pull-out the target region by class number and generate final draft genome
"""

## Help screen ##
if len(sys.argv) == 1 or sys.argv[1] in ['-h','--help']:
    print esom_help
    sys.exit()
elif sys.argv[1] in ['-h','--help']:
    print esom_help
    sys.exit()

bin_dir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
pkg_home = os.path.dirname(bin_dir)

esom_bin = os.path.join(settings.esom_path,'bin')

## Make sure that esom path is set correctly in the esomstart and esomtrn ##
esomstart = os.path.join(settings.esom_path,"bin","esomstart")
esomstart_declare = file_grep("ESOM_HOME=",esomstart)
current_path = re.search('[\"]([^\"]+)[\"]',esomstart_declare).group(1)
if current_path is not settings.esom_path:
    update_ESOM_HOME (esomstart,settings.esom_path)

esomtrn = os.path.join(settings.esom_path,"bin","esomtrn")
esomtrn_declare = file_grep("ESOM_HOME=",esomtrn)
current_path = re.search('[\"]([^\"]+)[\"]',esomtrn_declare).group(1)
if current_path is not settings.esom_path:
    update_ESOM_HOME (esomtrn,settings.esom_path)

if sys.argv[1] == "train":
    arguments = sys.argv[2:]
    call = os.path.join(pkg_home,'bin','train_esom.py')
    arguments.insert(0,call)
    py = sys.executable
    arguments.insert(0,py)
    subprocess.call(arguments)

elif sys.argv[1] == "annotate":
    arguments = sys.argv[2:]
    call = os.path.join(pkg_home,'bin','annotate_esom.py')
    arguments.insert(0,call)
    py = sys.executable
    arguments.insert(0,py)
    subprocess.call(arguments)

elif sys.argv[1] == "extract":
    arguments = sys.argv[2:]
    call = os.path.join(pkg_home,'bin','extract_esom.py')
    arguments.insert(0,call)
    py = sys.executable
    arguments.insert(0,py)
    subprocess.call(arguments)
