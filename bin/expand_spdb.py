#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  8 14:11:10 2018

@author: kevinamses
"""
#%% imports
import sys
import os
import argparse
import settings
import re
import ast
import inspect
from itertools import izip
from lib import subprocessC, start_logging, replace_line_by_pattern, file_yield_lines

bin_dir =os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
this_module = "spexpand"

#%% argument handling
parser = argparse.ArgumentParser()
parser.add_argument('-db','--spdb', metavar="swissprot_fasta",action="store", default=settings.path_to_spdb, required=False, help = "A file in FASTA format containing the swissprot database to be used. FASTA headers must be in standard swissprot format [... OS= ... GN= | PE= | OX= | EOL ]. See README for more information and examples.")
parser.add_argument('-t','--taxdb', metavar='taxonomy_db', action='store', default=settings.path_to_taxdb, required=False, help='The path to the taxonomy database associated with -db|--spdb')
parser.add_argument('-p','--proteins', metavar='proteins_to_add', action='store', required=True, help='A file in FASTA format that contains the protein sequences that you would like to add to the swissprot-style database.')
parser.add_argument('-l','--lineages', metavar='lineages_to_add', action='store', required=True, help='A two-column, tab-separated list of all unique OSs present in -p|--proteins (column 1) and semicolon-separated lineage information (column2).')
parser.add_argument('-o','--output', metavar='output', action='store', required=True, help='The output path for your expanded swissprot-style database.')
parser.add_argument('-d','--defaults', action='store_true', required=False, help='Include this flag if you would like scgid to update package settings to reflect the newly created databases as your defaults.')
args = parser.parse_args()

spdb_start = os.path.abspath(args.spdb)
taxdb_start = os.path.abspath(args.taxdb)
seqs_to_add = os.path.abspath(args.proteins)
lin_to_add = os.path.abspath(args.lineages)
spdb_end = os.path.abspath(args.output)
taxdb_end = spdb_end+'.taxdb'

outdir = os.path.split(os.path.abspath(args.output))[0]

#%% cd to output dir and initiate logging
try:
    os.chdir(outdir)
except:
    os.mkdir(outdir)
    os.chdir(outdir)
logs = start_logging(this_module, args, sys.argv)
logger = logs[0]

#%% Some checks and balances

if not os.path.isfile(spdb_start):
    logger.critical(args.spdb+" does not exist!")
    raise IOError
if not os.path.isfile(taxdb_start):
    logger.critical(args.taxdb+" does not exist!")
    raise IOError
if not os.path.isfile(seqs_to_add):
    logger.critical(args.proteins+" does not exist!")
    raise IOError
if not os.path.isfile(lin_to_add):
    logger.critical(args.lineages+" does not exist!")
    raise IOError

#%% First thing to do it concatenate spdb_start and seqs_to_add into spdb_out... grabbing all species designations (sp_os) along the way
logger.info("Concatenating new protein sequences and current swissprot-style database into expanded database.")
with open(spdb_end, 'w') as end:
    for fi in [spdb_start, seqs_to_add]:
        for line in file_yield_lines(fi):
            end.write(line)
logger.info("Expanded swissprot-style database written to: %s" % (spdb_end))

#%% Process user-specfied lineages into a dictionary

logger.info("Processing your supplied additional lineage data for taxdb and reading current taxdb...")
new_lin_dict = {}
for line in file_yield_lines(lin_to_add):
    spl = line.split('\t')
    if len(spl) != 2 or len(spl[1].split(";")) < 2:
        logger.critical("Lineage file is not formatted correctly. Must be a two-column, tab-separated list with species designation in column 1 and semicolon-separated lineage information in column 2.")
        raise ValueError
    else:
        new_os = spl[0].strip()
        new_lin = spl[1].strip()
        new_lin_dict[new_os] = new_lin

#%% Load the current taxonomy database and combine with user-specified lineages into new taxdb
with open(taxdb_start) as t:
    db_dict = t.read()
    taxdb_new = ast.literal_eval(db_dict)

taxdb_new.update(new_lin_dict)

with open(taxdb_end, 'w') as t:
    i=1
    t.write("{\n")
    for sp_os, lineage in taxdb_new.iteritems():
        t.write("'%s': '%s'" % (sp_os, lineage))
        if i != len(taxdb_new):
            t.write(",\n")
        else:
            t.write("\n}")
        i+=1

logger.info("Expanded taxonomy database written to %s" %(taxdb_end))
#%% Update settings to reflect newly-created databases as defaults

if args.defaults:
    logger.info("Updated settings.py to reflect these databases as your new defaults.")
    replace_line_by_pattern(os.path.join(bin_dir,"settings.py"), "path_to_spdb=", "path_to_spdb=\"%s\"" % (spdb_end))
    replace_line_by_pattern(os.path.join(bin_dir,"settings.py"), "path_to_taxdb=", "path_to_taxdb=\"%s\"" % (taxdb_end))

logger.info("scgid spexpand DONE")

#%%














