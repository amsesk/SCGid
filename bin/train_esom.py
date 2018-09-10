#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 16:31:07 2017

@author: kevinamses
"""

import sys
import argparse
import os
import subprocess
import inspect
import math
import re
import settings
from lib import start_logging, specify_Xmx_esom, subprocessP, subprocessC, do_long_wait

#%%

## Argument Handling ##
parser = argparse.ArgumentParser()
parser.add_argument('-n','--nucl', metavar = "contig_fasta", action="store",required=True, help = "A FASTA file containing the nucleotide assembly. (MANDATORY)")
parser.add_argument('-f','--prefix', metavar = 'prefix_for_output', required=False, default='scgid', help="The prefix that you would like to be used for all output files. DEFAULT = scgid")

## print_tetramer_freqs options
parser.add_argument('-m','--mintig', metavar = "minimum_contig_size", action="store",required=False, default="1000", help = "Contig size cutoff (in nucleotides) for inclusion in the ESOM training. Default = 1000 bp")
parser.add_argument('-w','--window', metavar = "window_size", action="store",required=False, default="1000", help = "Size of the window in which kmer frequencies are calculated. Default = 1000 bp")
parser.add_argument('-k','--kmer', metavar = "kmer_size", action="store",required=False, default="4", help = "Kmer size for which frequencies will be calculated. Default = 4")

## esomtrn options
parser.add_argument('-r','--rows', metavar = "rows_in_map", action="store",required=False, help = "The number of rows to be present in the output ESOM. Default = 5.5x the number of neurons")
parser.add_argument('-c','--cols', metavar = "columns_in_map", action="store",required=False, help = "The number of columns to be present in the output ESOM. Default = 5.5x the number of neurons")

parser.add_argument('-sr','--start_radius', metavar = "start_radius", action="store",required=False, default = '50', help = "Start radius for the ESOM. (MANDATORY)")
parser.add_argument('-e','--epochs', metavar = "training_epochs", action="store",required=False, default = "20", help = "Number of epochs to train over. Default = 20")
parser.add_argument('--Xmx', metavar = "available_memory", action="store",required=False, default = "512m", help = "Set memoray available to train the ESOM. Specicy as such: X megabytes = Xm, X gigabytes = Xg")

args =  parser.parse_args()
esom_bin = os.path.join(settings.esom_path,'bin')
bin_dir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
pkg_home = os.path.dirname(bin_dir)

#%%
this_module = "kmers"
prefix = args.prefix
assembly = os.path.abspath(args.nucl)
file_prefix = os.path.split(assembly)[1]


#%%
try:
    os.chdir(args.prefix+'_scgid_output')
except:
    os.mkdir(args.prefix+'_scgid_output')
    os.chdir(args.prefix+'_scgid_output')

logger = start_logging(this_module, sys.argv)

try:
    os.chdir('esom')
except:
    os.mkdir('esom')
    os.chdir('esom')

#%% Let's get started
logger.info("Starting kmer-based genome selection process.")
logger.info("Calculating 4-mer frequencies across scaffolds and generating ESOM input files.")
## print tetramer freqs ##
arguments = ['cp',assembly,file_prefix]
subprocess.call(arguments)
call = os.path.join(bin_dir,'print_tetramer_freqs_deg_filter_esom_VD.pl')
arguments = ['perl',call,'-s',file_prefix,'-m',args.mintig,'-w',args.window,'-k',args.kmer]
logger.info(' '.join(arguments))
do_long_wait(lambda: subprocessP(arguments, logger), 'none')

## If rows and cols aren't provided, make a square map that is big enough for number of neurons #
if args.rows is None and args.cols is None:
    logger.info("Nothing specified for -r|--rows or -c|--cols. Using defult dimensions for map.")
    nodes = len(open(file_prefix+".lrn").readlines()) - 4 ## number of lines in .lrn file minus 4 header lines
    neurons = float(nodes*5.5)
    dim = float(math.sqrt(neurons))
    rows = str(round(dim)).split('.')[0]
    cols = str(round(dim)).split('.')[0]
else:
    logger.info("Using user-specified dimensions for map: "+arg.rows+"x"+args.cols)
    rows = args.rows
    cols = args.cols

## Train ESOM ##
out = file_prefix+"."+rows+"x"+cols+"e"+args.epochs+".wts"
b = file_prefix+"."+rows+"x"+cols+"e"+args.epochs+".bm"

#Alter hard-coded memory-availability in esomstart to match user-specified resources
specify_Xmx_esom(os.path.join(esom_bin,"esomstart"),args.Xmx)

##Train the ESOM
logger.info("Starting ESOM training process")
train_args = [os.path.join(esom_bin,"esomtrn"), "--permute", "--out", out,"-b", b,"--cls", file_prefix+".cls", "--lrn",file_prefix+".lrn", "--algorithm", "kbatch", "--rows", rows, "--columns", cols, "-bmc", "6", "--start-radius", args.start_radius, "--epochs", args.epochs, "-k", "0.15", "--bmsearch", "standard", "--dist", "euc"]

logger.info(' '.join(train_args))

do_long_wait(lambda: subprocessP(train_args, logger), 'none')

logger.info("ESOM trained successfully, time to annotate.")

## DONE
