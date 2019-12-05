#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  5 17:58:57 2018

@author: kevinamses
"""

import argparse
import os
import subprocess
import sys
from lib import logger1, log_command_line_errors, read_lines_from_large_file, subprocessP, subprocessC, subprocessT

parser = argparse.ArgumentParser()

parser.add_argument('-nt','--nontarget', metavar = "nontarget_fasta", action="store",required=True, help = "A FASTA file containing the nontarget contigs output by the blob, rscu, esom, or consensus methods. (MANDATORY)")
parser.add_argument('-pe1','--paired1', metavar = 'raw_reads_directory', required=True, help="The path to the directory containing the raw reads used in assembly.")
parser.add_argument('-pe2','--paired2', metavar = 'raw_reads_directory', required=True, help="The path to the directory containing the raw reads used in assembly.")
parser.add_argument('-u','--unpaired', metavar = 'raw_reads_directory', required=True, help="The path to the directory containing the raw reads used in assembly.")
parser.add_argument('-f','--prefix', metavar = 'prefix_for_output', required=False, default='scgid', help="The prefix that you would like to be used for all output files. DEFAULT = scgid")
parser.add_argument('--cpus', metavar = 'cores', action = 'store', required = False, default = "1", help = "The number of cores available for bowtie2 to use for alignment.")

args =  parser.parse_args()

ntfasta = os.path.abspath(args.nontarget)
prefix = args.prefix
this_module = "purify"


## get absolute paths of reads ##
pe1 = []
pe2 = []
u = []
for f in args.paired1.split(','):
    pe1.append(os.path.abspath(f))
for f in args.paired2.split(','):
    pe2.append(os.path.abspath(f))
for f in args.unpaired.split(','):
    u.append(os.path.abspath(f))

## Get into output directory ##
try:
    os.chdir(args.prefix+'_scgid_output')
except:
    os.mkdir(args.prefix+'_scgid_output')
    os.chdir(args.prefix+'_scgid_output')

logger = start_logging(this_module, sys.argv)

try:
    os.chdir('purify')
except:
    os.mkdir('purify')
    os.chdir('purify')

## check for bowtie2 availability
try:
    subprocessT(['bowtie2','-help'])
except:
    logger.critical("Bowtie2 not available in PATH. Exitting...")
    sys.exit(1)

## Build the bowtie database ##
logger.info("Building bowtie database")
cmd = ['bowtie2-build', ntfasta, prefix+'_nontarget_ref']
logger.info(' '.join(cmd))
p = subprocessP(cmd, logger)


## Align raw reads to the bowtie database ##
logger.info("Aligning raw reads to bowtie database")
cmd = ['bowtie2', '-x', prefix+'_nontarget_ref', '-1', ','.join(pe1), '-2', ','.join(pe2), '-U', ','.join(u), '-S', prefix+'_nontarget_read_alignment.sam', '-q', '--end-to-end', '--sensitive','--threads',args.cpus]
logger.info(' '.join(cmd))
p = subprocessP(cmd, logger)

## Identify unaligned (ie target) reads and print to fastq ##
logger.info("Processing SAM alignment and printing target reads to fastq format")

if os.path.isfile(prefix+'_purified_pe1.fastq'):
    os.remove(prefix+'_purified_pe1.fastq')
if os.path.isfile(prefix+'_purified_pe2.fastq'):
    os.remove(prefix+'_purified_pe2.fastq')
if os.path.isfile(prefix+'_purified_orphan.fastq'):
    os.remove(prefix+'_purified_orphan.fastq')

#for line in open(prefix+'_nontarget_read_alignment.sam','r').readlines():
for line in read_lines_from_large_file(prefix+'_nontarget_read_alignment.sam'):
    spl = line.split('\t')
    if len(spl) < 4:
        continue
    else:
        if spl[1] == '77':
            with open(prefix+'_purified_pe1.fastq','a') as out:
                out.write('@'+spl[0]+'\n'+spl[9]+'\n+\n'+spl[10]+'\n')
        elif spl[1] == '141':
            with open(prefix+'_purified_pe2.fastq','a') as out:
                out.write('@'+spl[0]+'\n'+spl[9]+'\n+\n'+spl[10]+'\n')
        elif spl[1] == '4':
            with open(prefix+'_purified_orphan.fastq','a') as out:
                out.write('@'+spl[0]+'\n'+spl[9]+'\n+\n'+spl[10]+'\n')
        else:
            continue

