#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 16:44:52 2017

@author: kevinamses
"""
import argparse
import subprocess
import os
import inspect
import settings
from lib import *
import shutil

bin_dir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
pkg_home = os.path.dirname(bin_dir)
path_to_this_file = inspect.getfile(inspect.currentframe())
this_file = os.path.split(path_to_this_file)[1]
this_module = "predict_and_blast"

### Argument parsing ###

parser = argparse.ArgumentParser()
parser.add_argument('-n','--nucl', metavar = "contig_fasta", action="store",required=True,help = "A FASTA file containing the genome assembly.")
parser.add_argument('-sp','--augustus_sp', metavar = "augustus_species", action="store",required=False, help = "Augustus species for gene predicition.")
parser.add_argument('-f','--prefix', metavar = 'prefix_for_output', required=False, default='scgid', help="The prefix that you would like to be used for all output files. DEFAULT = scgid")
parser.add_argument('-c','--coding', metavar = "what_to_do", action="store", required=False, default="y,y", help = "Which commands need to be run.")

# BLAST args

parser.add_argument('-p','--prot', metavar = "protein_fasta", action="store",required=True,help = "A FASTA file containing the proteins called from the genome.")
parser.add_argument('-e', '--evalue', metavar = 'e-value_cutoff', action = 'store', required = False, default = '1e-10', help = "The evalue cutoff for blast. Default: 1xe-10)")
parser.add_argument('-db', '--spdb', metavar = 'swissprot_db', action='store', required=False, help = "The path to and prefix of your swissprot blast database - can be compiled into blastdb OR in FASTA format.")
parser.add_argument('--cpus', metavar = 'cores', action = 'store', required = False, default = 1, help = "The number of cores available for blastp to use.")
args = parser.parse_args()

### Variable declaration ###
prefix = args.prefix
nucl = args.nucl
prot = args.prot
blastp_evalue = args.evalue
cores = args.cpus

logger = start_logging("predict_and_blast", sys.argv)

if args.coding not in ['1,1','1,0','0,1','0,0']:
    logger.critical("Invalid argument value for -c|--coding")
    raise ValueError
else:
    to_do = args.coding.split(",")


'''
## check and ammend file structure ##
if os.path.isdir(prefix+'_scgid_output') is False:
    os.mkdir(prefix+'_scgid_output')
if os.path.isdir(prefix+'_scgid_output/blob') is False:
    os.mkdir(prefix+'_scgid_output/blob')
'''

### Do we need to predict proteins from nucleotide fasta...? ###
if to_do[0] is '1':

    aug_sp = args.augustus_sp

# make sure that augustus is available in path
    if settings.AUGUSTUS_CONFIG_PATH is not None:
        AUGUSTUS_CONFIG_PATH = settings.AUGUSTUS_CONFIG_PATH
    else:
        try:
            AUGUSTUS_CONFIG_PATH = os.environ["AUGUSTUS_CONFIG_PATH"]
        except:
            logger.exception("Global variable $AUGUSTUS_CONFIG_PATH is not present in environment and is unset in settings.py. Unable to run augustus.")
            raise ValueError("Global variable $AUGUSTUS_CONFIG_PATH is not present in environment and is unset in settings.py. Unable to run augustus")

    try:
        subprocessT(['augustus','--AUGUSTUS_CONFIG_PATH='+AUGUSTUS_CONFIG_PATH,'--help'])
    except:
        logger.exception("Augustus is either not installed or unavailable.")
        raise IOError("Augustus is either not installed or unavailable.")

# Run augustus...
    logger.info("Predicting proteins with augustus...")

    run_aug_cmd = ['augustus','--AUGUSTUS_CONFIG_PATH='+AUGUSTUS_CONFIG_PATH,'--species='+aug_sp,'--gff3=on', nucl,'--uniqueGeneId=true']

    logger.info(' '.join(run_aug_cmd))
    out = subprocessP(run_aug_cmd, logger)

    logger.info("Writing augustus output to "+prefix+".aug.out.gff3")
    with open(prefix+'.aug.out.gff3','w') as f:
        f.write(out)


    gff3_to_fasta(prefix+'.aug.out.gff3', prefix+'.aug.out.fasta')
    logger.info("gff3 converted to fasta and writen to "+os.path.join(os.getcwd(),prefix+'.aug.out.fasta'))

    ##move gff3 to default locations if finished successfully in temp
    shutil.copyfile(prefix+'.aug.out.gff3','../'+prefix+'.aug.out.gff3')
    shutil.copyfile(prefix+'.aug.out.fasta','../'+prefix+'.aug.out.fasta')


### Do we need to BLAST...? ###
if to_do[1] is '1':

    path_to_spdb = args.spdb

    ### Check to make sure blast is available ###
    try:
        subprocessT(['blastp','-help'])
    except:
        logger.exception("BLAST is either not installed or unavailable.")
        raise IOError("BLAST is either not installed or unavailable.")

    ## Start BLASTing...
    logger.info("Blasting predicted proteins against the swissprot database located at: "+path_to_spdb)

    ## Make sure the blast database has been created, or create it
    fasta_ext = ['faa','fas','fasta']
    default = '.'.join([path_to_spdb,'phr'])
    print default
    if os.path.isfile(default):
        logger.info("Blast database detected at "+path_to_spdb)

    elif os.path.isfile(path_to_spdb):
        db_ext = path_to_spdb.split('.')[-1]
        if is_fasta(path_to_spdb):
            new_dbname = os.path.split(path_to_spdb)[1]
            logger.info("No blast database detected at "+path_to_spdb+". Building...")
            build_cmd = ['makeblastdb','-in',path_to_spdb,'-dbtype','prot','-title',new_dbname,'-out',path_to_spdb]
            logger.info(' '.join(build_cmd))
            subprocessP(build_cmd, logger)
        else:
            logger.critical("-db|--spdb denoted as FASTA format, but has been determined to not be in FASTA format")
            raise IOError
    else:
        for ext in fasta_ext:
            looking = '.'.join([path_to_spdb, ext])
            if os.path.isfile(looking):
                if is_fasta(looking):
                    path_to_spdb = looking
                    new_dbname = os.path.split(path_to_spdb)[1]
                    logger.info("No blast database detected at "+path_to_spdb+", but FASTA detercted at "+looking+". Building database from it...")
                    build_cmd = ['makeblastdb','-in',path_to_spdb,'-dbtype','prot','-title',new_dbname,'-out',path_to_spdb]
                    logger.info(' '.join(build_cmd))
                    subprocessP(build_cmd, logger)
                    break
            if fasta_ext.index(ext) == len(fasta_ext)-1:
                logger.critical("No blast database, nor FASTA with which to create one, found at "+path_to_spdb+".")
                raise IOError


## Move on to blastp once blast database has been found or created

    blastp_cmd = ['blastp','-query',prot,'-max_target_seqs','1','-evalue',blastp_evalue,'-db',path_to_spdb,'-outfmt','6','-out',prefix+'.spdb.blast.out','-num_threads',cores]
    logger.info(' '.join(blastp_cmd))
    subprocessP(blastp_cmd, logger)

    ##move blastout to default locations if finished successfully in temp
    shutil.copyfile(prefix+'.spdb.blast.out','../'+prefix+'.spdb.blast.out')

