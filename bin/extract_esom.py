#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 20:01:11 2017

@author: kevinamses
"""
import sys
import argparse
import os
import subprocess
import inspect
from lib import *

bin_dir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
pkg_home = os.path.dirname(bin_dir)

parser = argparse.ArgumentParser()
parser.add_argument('-n','--nucl', metavar = "contig_fasta", action="store",required=True, help = "A FASTA file containing the nucleotide assembly. (MANDATORY)")
parser.add_argument('-c', '--cls', metavar = 'class_file', action = 'store', required = True, help = 'The .cls output file from "esom train".')
parser.add_argument('-nf', '--names_file', metavar = 'names_file', action = 'store', required = False, default = None, help = 'The .names output file from "esom train".')
parser.add_argument('-cid', '--classnum', metavar = 'class_number', action = 'store', required = True, help = "The class number of interest. That is, the class that represents the selection of target contigs you made in esomana.")
parser.add_argument('-f','--prefix', metavar = 'prefix_for_output', required=False, default='scgid', help="The prefix that you would like to be used for all output files. DEFAULT = scgid")
parser.add_argument('-l','--loyal', metavar = 'loyalty_threshold', required=False, default='51', help="The loyalty threshold for keeping a contig based on where its various windows end-up in ESOM. DEFAULT = 51")
args = parser.parse_args()

output_dir = args.prefix+'_scgid_output'

nucl = os.path.abspath(args.nucl)
prefix = args.prefix
cls_file = os.path.abspath(args.cls)

#%% Look for files in default locations if nothing given for -nf|--names_file
if args.names_file is None:
    names_file = os.path.join(output_dir,'esom',os.path.split(nucl)[1]+'.names')
    if not os.path.isfile(names_file):
        raise IOError("Nothing given for -nf|--names_file and .names file not found in default location.")

else:
    names_file = os.path.abspath(args.names_file)
    if not os.path.isfile(names_file):
        raise IOError("No such .names file, "+names_file)

#%% navigate to head of working directory
try:
    os.chdir(output_dir)
except:
    os.mkdir(output_dir)
    os.chdir(output_dir)

#%% navigate to esom output directory
try:
    os.chdir('esom')
except:
    os.mkdir('esom')
    os.chdir('esom')

#%% generate fasta for esom class of interest

cmd = ['perl',os.path.join(bin_dir,'getClassFasta.pl'),'-cls',cls_file,'-names',names_file,'-num',args.classnum,'-fasta',nucl,'-loyal',args.loyal]
print ' '.join(cmd)
mv = ['mv',args.classnum+'.fasta',prefix+'_esom_final_genome.fasta']
#rm_conf = ['rm',args.classnum+'.conf']
subprocessC(cmd)
subprocessC(mv)
#subprocess.call(rm_conf)
