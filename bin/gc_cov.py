#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 15:53:11 2017

@author: Kevin Amses (amsesk@umich.edu)
"""
#%% imports
import sys
import os
import subprocess
import operator
import inspect
import argparse
import re
import ast
import shutil
import logging
import pandas as pd
import numpy as np
import cPickle as pickle
from sequence import *
from lib import *
from infotable import infotable
import settings
#import matplotlib.pyplot as plt

#%%
## Argument Handling ##

# ALWAYS REQUIRED
parser = argparse.ArgumentParser()
parser.add_argument('-n','--nucl', metavar = "contig_fasta", action="store",required=True,help = "A FASTA file containing the genome assembly.")
parser.add_argument('-t','--taxdb', metavar = "taxonomy_db", action="store", required=False, default=settings.path_to_taxdb, help = "The location of the taxonomy database, likely provided by an earlier script.")
parser.add_argument('-s', '--stringency', metavar = "stringency_threshold", required=False, default=0.05, help = "The proportion of annotated non-target points that are allowed to be included in the final selection window. IMPORTANT NOTE: The non-target-annotated points can still be throw-out of the final genome fasta by specifying the --toss_nontarget option.")
parser.add_argument('-f','--prefix', metavar = 'prefix_for_output', required=False, default='scgid', help="The prefix that you would like to be used for all output files. DEFAULT = scgid")
parser.add_argument('-g', '--targets', metavar = 'target_taxa', action='store', required=True, help="A comma-separated list with NO spaces of the taxonomic levels that the gc-coverage window should be chosen with respect to including. EXAMPLE: '-g Fungi,Eukaryota,Homo'")
parser.add_argument('-x', '--exceptions', metavar = 'exceptions_to_target_taxa', action='store', required=False, default=None, help="A comma-separated list with NO spaces of any exlusions to the taxonomic levels specified in -g|--targets. For instance if you included Fungi in targets but want to exclude ascomycetes use: '-x Ascomycota'")

parser.add_argument('-sp','--augustus_sp', metavar = "augustus_species", action="store",required=False, default=None, help = "Augustus species for gene predicition.")
parser.add_argument('-e', '--evalue', metavar = 'e-value_cutoff', action = 'store', required = False, default = '1e-10', help = "The evalue cutoff for blast. Default: 1xe-10)")
parser.add_argument('-db', '--spdb', metavar = 'swissprot_fasta', action='store', required=False, default=settings.path_to_spdb,  help = "The path to your version of the swissprot database in FASTA format.")
parser.add_argument('--cpus', metavar = 'cores', action = 'store', required = False, default = '1', help = "The number of cores available for blastp to use.")

# MAYBE PROVIDED BY EARLIER PART OF SCRIPT
parser.add_argument('-b','--blastout', metavar = "blastout", action="store",required=False, help = "The blast output file from a search of the swissprot database with your proteins as query. Defaults to outfmt=6 and max_target_seqs=1. Provided by earlier script.")
parser.add_argument('-p','--prot', metavar = "protein_fasta", action="store",required=False, help = "A FASTA file containing the proteins called from the genome.")
args =  parser.parse_args()

#%%
bin_dir =os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
pkg_home = os.path.dirname(bin_dir)
path_to_this_file = inspect.getfile(inspect.currentframe())
this_module = "gc-cov"

nucl_path = os.path.abspath(args.nucl)
path_to_taxdb = os.path.abspath(args.taxdb)
path_to_spdb = os.path.abspath(args.spdb)
prefix = args.prefix
stringency_thresh = args.stringency

prot = args.prot
blastout = args.blastout

if prot is not None:
    prot = os.path.abspath(prot)
if blastout is not None:
    blastout = os.path.abspath(blastout)

blastout_check = True
augout_check = True

#%% navigate to head of working directory
try:
    os.chdir(args.prefix+'_scgid_output')
except:
    os.mkdir(args.prefix+'_scgid_output')
    os.chdir(args.prefix+'_scgid_output')

# Start logging now that we are in run working directory
logger = start_logging(this_module, sys.argv)

#%% in the absense of -p|--prot, see if a gff3 of predicted proteins was previously produced by rscu.py or there if a fasta is present here
if prot is None:
    prot = os.path.join(os.getcwd(),'blob',prefix+'.aug.out.fasta') #check for fasta here first
    if not os.path.isfile(prot):
        try_gff3 = os.path.join(os.getcwd(),'rscu',prefix+'.aug.out.gff3') #check for gff3 in rscu next, leave prot equal to default aug.out.fasta name - will work for either case if|else below
        if os.path.isfile(try_gff3):
            logger.info("Nothing given for -p|--prot, however found gff3 produced by rscu module at, "+prot+". Converting to fasta...")
            gff3_to_fasta(try_gff3, os.path.join(os.getcwd(),'blob',prefix+'.aug.out.fasta'))
        else:
            logger.info("Nothing given for -p|--prot and unable to locate augustus fasta or gff3 output files in default locations. Augustus will have to be run prior to blob-based genome selection.")
            augout_check = False #Will have to run Augustus below

else:
    if not os.path.isfile(prot):
        #have to raise IOError, user provided a nonexistent file, augout_check = False, in spirit but unecessary
        logger.critical("-p|--prot: No such file or directory, "+prot)
        raise IOError
    else:
        pass #augout_check set to True by default

#%% see if spdb.blast.out has been done here already
if blastout is None:
        #try to locate spdb blastout in default location
        blastout = os.path.join(os.getcwd(),'blob',prefix+'.spdb.blast.out') #check here first
        if not os.path.isfile(blastout):
            logger.info("Nothing given for -b|--blastout and unable to locate blast output file in default location. BLASTP search will have to be run prior to blob-derived genome selection.")
            blastout_check = False #Will have to run BLASTP below

else:
    if not os.path.isfile(blastout):
        #have to raise IOError, user provided a nonexistent file, blastout_check = False, in spirit but unecessary
        logger.critical("-b|--blastout: No such file or directory, "+blastout)
        raise IOError
    else:
        pass #blastout_check set to True by default
#%% navigate to the blob directory

try:
    os.chdir('blob')
except:
    os.mkdir('blob')
    os.chdir('blob')

#%% make a temp dir to deal with incomplete runs leaving incomplete files
if os.path.isdir('temp'):
    shutil.rmtree('temp')
    os.mkdir('temp')
    os.chdir('temp')
else:
    os.mkdir('temp')
    os.chdir('temp')

#%%

## Let's get started ##
logger.info('Starting GC-Coverage-based genome selection')

#%% Protein prediction and blasting
if not blastout_check or not augout_check: #this means we have to predict and/or blast first
    if not augout_check and args.augustus_sp is None:
        raise ValueError("ERROR: -sp|--augustus_sp required for gene prediction.")
    if not blastout_check and args.spdb is None:
        raise ValueError("ERROR: -db|--spdb required for BLASTP searching of predicted genes.")

    #tell predict_and_blast.py what it needs to do
    to_do = ['1','1']
    if augout_check:
        to_do[0] = '0'
        logger.info("Found protein fasta: "+prot+". Skipping augustus...")
    if blastout_check:
        to_do[1] = '0'
        logger.info("Found blast output file: "+blastout+". Skipping BLAST...")
    to_do_arg = ','.join(to_do)
    
    if args.augustus_sp is None:
        args.augustus_sp = "NA"

    arguments = ['-n',nucl_path,'-sp',args.augustus_sp,'-e',args.evalue,'-db',args.spdb,'--cpus',args.cpus,'-f',args.prefix,'-c',to_do_arg,'-p',prot]
    call = os.path.join(bin_dir,'predict_and_blast.py')
    arguments.insert(0,call)
    py = sys.executable
    arguments.insert(0,py)
    subprocessC(arguments)

    logger.info("Gene model prediction and blastp against swissprot-style database complete. Returning to GC-Coverage plot anaysis.")

else:
    logger.info("Found protein fasta: "+prot)
    logger.info("Found blast output file: "+blastout)

# get best blast hit
output = []
for key,best in best_blast_hit(blastout).iteritems():
    blastout_bestout = prefix+'.spdb.blast.out.best'
    output.append('\t'.join(best))

output='\n'.join(output)
with open(blastout_bestout,'w') as f:
    f.write(output)

logger.info("Best blast hits written to "+os.path.join(os.getcwd(),prefix+'.spdb.blast.out.best'))


# parse blastout
parsed_blastout_out = prefix+'.spdb.blast.out.best.parsed'
logger.info("Writing best, parsed blast output to "+parsed_blastout_out)
parsed_blastout = parse_spdb_blastout(path_to_spdb, blastout_bestout, logger)
with open(parsed_blastout_out,'w') as f:
    f.write('\n'.join(parsed_blastout))

#set blastout to the parsed,best blastout for further processing
blastout = os.path.abspath(parsed_blastout_out)

#%%

#Get targetting information from arguments
target_taxa = {}
target_taxa['target'] = args.targets.split(',')
if args.exceptions is None:
    target_taxa['exception'] = []
else:
    target_taxa['exception'] = args.exceptions.split(',')

nucl = pkl_fasta_in_out(nucl_path)
proteins = pkl_fasta_in_out(prot, seq_type='prot', contig_info=True)

logger.info("Nucleotide and protein fastas read-in successfully.")

# read-in taxdb or create one

if path_to_taxdb is None:
    pass
    ##### Download all-taxonomy and run build_spex-taxdb.py #####
    #arguments = ['-db',,'-t',,'-o',,'-q','-a',]
    call = os.path.join(bin_dir,'build_spex_taxdb.py')
    arguments.insert(0,call)
    py = sys.executable
    arguments.insert(0,py)
    subprocessC(arguments)

else:
    with open(path_to_taxdb,'r') as f:
        s = f.read()
        taxdb = ast.literal_eval(s)
    logger.info("Taxonomy database imported successfully.")

#### BUILD INFOTABLE ####
## (1st) Gather contig and protein information from sequence objects
ldict = [] 
gcs = {}
covs = {}
for i in nucl:
    node_num = "_".join(i.label.split("_")[0:2])+"_"
    gcs[node_num] = i.gc
    covs[node_num] = i.coverage
attrs = {}
for i in proteins:
    pid = reverse(reverse(i.label).split(".")[0])
    node_num = "_".join(i.label.split("_")[0:2])+"_"
    attrs[i.label] = [i.length,covs[node_num],gcs[node_num],pid]

## (2nd) Grab hit info from the parsed blast output file and use OS=<blah> to grab lineage information from the taxdb
pattern = re.compile(SPDB_OS_REGEXP_PATTERN)
for line in open(blastout).readlines():
    spl = line.strip().split("\t")
    query = spl[0]
    desc = spl[2]
    evalue = spl[3]
    s = re.search(pattern, desc)
    if s is None:
        logger.critical("Internal error. Incorrectly formatted parsed BLAST output.")
        raise IOError
    sp_os = s.group(1).strip()

    ## Raise useful error message if there is a mismatch between taxbd and blastout
    try:
        lineage = taxdb[sp_os]
    except KeyError:
        logger.warning("'"+sp_os+"' is missing from the taxdb. You probably specified a -t|--taxdb file built from a likely out-of-date swissprot database other than the one you blasted against. Rerun build_taxdb or specify a different taxdb to correct. Lineage set as 'Not_in_taxdb' for this run.")
        lineage = "Not_in_taxdb"
    
    # Can get rid of these once local taxdb is rebuilt
    sp_os = sp_os.replace("'","")
    sp_os = sp_os.replace("#","")
    
## (3rd) For each row, make a dictionary and append to ldict list
    newrow = {
            'contig': query,
            'prot_len': attrs[query][0],
            'coverage': attrs[query][1],
            'gc': attrs[query][2],
            'pid': attrs[query][3],
            'sp_os': sp_os,
            'lineage': lineage,
            'evalue': evalue,
            'parsed_lineage':None
            }
    ldict.append(newrow)

## (Finally) Put everything into an infotable object, df attribute and make some changes to columns with apply
info_table = infotable(target_taxa)
colnames = ["contig","prot_len","coverage","gc","pid","sp_os","lineage","evalue","parse_lineage"]
info_table.populate(ldict, colnames)

info_table.df.contig = info_table.df.contig.apply(node_num_only)
info_table.df.lineage = info_table.df.lineage.apply(str.replace,args=('; ','.'))
info_table.df.lineage = info_table.df.lineage.apply(str.split,args='.')
info_table.df.reset_index()

#Parse lineage info into target|nontarget|unclassified
info_table.parse_lineage()

''' replaced by infotable.parse_lineage()
parsed = []
for lin in info_table.lineage:
    if any(i in target_taxa['target'] for i in lin) is True and any(i in target_taxa['exception'] for i in lin) is False:
        parsed.append('target')
    elif lineage == "Not_in_taxdb":
        parsed.append('unclassified')
    else:
        parsed.append('nontarget')
info_table['parse_lineage'] = parsed
'''

## Output infotables (one for classifieds, one for unclassifieds

info_table.df.to_csv(prefix+'_info_table.tsv',sep='\t', index = False, header = False)

no_class_cols = ["contig","coverage","gc","taxonomy"]
no_class = pd.DataFrame(columns=no_class_cols)
for no_class_tig in [x for x in nucl if '_'.join(x.label.split('_')[0:2]) not in info_table.df.contig.values]:
    short_label = '_'.join(no_class_tig.label.split('_')[0:2])
    new_row = pd.DataFrame([[short_label, no_class_tig.coverage, no_class_tig.gc, "Unclassified"]],columns=no_class_cols)
    no_class = no_class.append(new_row)

no_class.to_csv(prefix+'_unclassified_info_table.tsv', sep = '\t', index = False, header = False)

logger.info("Protein information table generated and ready to go. Printed to "+os.path.join(os.getcwd(),prefix+".info_table.tsv"))
logger.info("GC and coverage information on protein-less contigs printed to "+os.path.join(os.getcwd(),prefix+'_unclassified_info_table.tsv'))

## Get target table (starting means for selection process
## Generate the selection windows 
target = info_table.df.loc[info_table.df['parse_lineage'] == 'target']
windows = generate_windows(info_table.df, target_taxa, target)

logger.info("GC-coverage windows successfully generated.")

## Pick the best window  
names = ["gc1cov1","gc1cov2",
         "gc2cov1","gc2cov2",
         "cov1gc1","cov1gc2",
         "cov2gc1","cov2gc2",
         "gc1noCov","gc2noCov",
         "cov1noGc","cov2noGc"]

window_stats = {
        'names': [],
        'tp': [],
        'ntp': [],
        'window': [],
        'gc_width': [],
        'cov_width': []
        }
i=0
#logger.info('\t'.join(['method','p(target)','p(nontarget)','gc_window','cov_window']))

with open(prefix+'.windows.all.out','w') as f:
    colnames = '\t'.join(['method','p_target', 'p_nontarget', 'gc_window', 'gc_window_width', 'cov_window', 'cov_window_width'])+'\n'
    f.write(colnames)

for win in windows:
    gc_window_table = get_window_table(info_table.df, win['gc'],"gc")
    window_table = get_window_table(gc_window_table, win['cov'],"coverage")

    t = float(window_table.loc[window_table['parse_lineage'] == 'target'].shape[0])
    all_t = float(info_table.df.loc[info_table.df['parse_lineage'] == 'target'].shape[0])
    nt = float(window_table.shape[0]-t)
    all_nt = float(info_table.df.shape[0]-all_t)
    disp_gc = str(round(win['gc'][0],4))+'-'+str(round(win['gc'][1],4))
    disp_cov = str(round(win['cov'][0],4))+'-'+str(round(win['cov'][1],4))
    gc_window_width = win['gc'][1] - win['gc'][0]
    cov_window_width = win['cov'][1] - win['cov'][0]
    #logger.info('\t'.join([names[i],str(round(t/all_t,4)), str(round(nt/all_nt,4)), disp_gc, disp_cov]))
    with open(prefix+'.windows.all.out','a') as f:
        #out = row+'\t'+'\t'.join(map(str,values))+'\n'
        out = '\t'.join([names[i],str(round(t/all_t,4)), str(round(nt/all_nt,4)), disp_gc, str(gc_window_width), disp_cov, str(cov_window_width)])+'\n'
        f.write(out)

    window_stats['names'].append(names[i])
    window_stats['tp'].append(t/all_t)
    window_stats['ntp'].append(nt/all_nt)
    window_stats['window'].append(win)
    window_stats['gc_width'].append(gc_window_width)
    window_stats['cov_width'].append(cov_window_width)
    i+=1

logger.info("Information on all windows printed to "+os.getcwd()+"/"+prefix+".windows.all.out")


## Note that this does not account for multiple maxima yet
keep_going = True
best_window_idx = None
best_window = None

while keep_going is True:
    if len(window_stats['window']) == 0:
        raise StopIteration("No best window at stringency threshold =",stringency_thresh," - raise it or add more options to target_taxa.")
    best_window_idx = window_stats['tp'].index(max(window_stats['tp']))
    if window_stats['ntp'][best_window_idx] <= stringency_thresh:
        keep_going = False

        # This is the best window based on the data and stringency threshold
        best_window = gc_cov_window(lower_gc=window_stats['window'][best_window_idx]['gc'][0],
                                    upper_gc=window_stats['window'][best_window_idx]['gc'][1],
                                    lower_cov=window_stats['window'][best_window_idx]['cov'][0],
                                    upper_cov=window_stats['window'][best_window_idx]['cov'][1],
                                    tailedness=window_stats['names'][best_window_idx],
                                    tp=window_stats['tp'][best_window_idx],
                                    ntp=window_stats['ntp'][best_window_idx])

        #print gc-cov plot and draw window lines
        print_plot_cmd = ['Rscript','--vanilla',os.path.join(bin_dir,'gc_cov.plot.R'),prefix+'_info_table.tsv',prefix+'_unclassified_info_table.tsv',str(window_stats['window'][best_window_idx]['gc']),str(window_stats['window'][best_window_idx]['cov']),prefix+"_gc_coverage_plots.pdf"]
        logger.info(" ".join(print_plot_cmd))
        out = subprocessP(print_plot_cmd, logger)

    else:
        window_stats['names'].pop(best_window_idx)
        window_stats['tp'].pop(best_window_idx)
        window_stats['ntp'].pop(best_window_idx)
        window_stats['window'].pop(best_window_idx)


logger.info("BEST WINDOW: "+str(best_window.__dict__))

## Decide on each contigs taxonomy based on all of its proteins
info_table.decide_taxonomy()

## Output list of contigs dropped by taxonomy alone
with open('excluded_by_sp_taxonomy','w') as ex:
    for nt in info_table.dump:
    #for nt in to_exclude
        ex.write(nt)
        ex.write('\n')

## Make gc-coverage decisions on unclassifieds for final draft genome
final_genome = []
non_target_bin = []
for tig in nucl:
    #if tig.name in to_add:
    if tig.name in info_table.keep:
        final_genome.append(tig)
    #elif tig.name in to_exclude:
    elif tig.name in info_table.dump:
        non_target_bin.append(tig)
    elif (tig.gc >= best_window.lower_gc and tig.gc <= best_window.upper_gc) and (tig.coverage >= best_window.lower_cov and tig.coverage <= best_window.upper_cov):
        final_genome.append(tig)
    else:
        non_target_bin.append(tig)

## Output final draft genome to fasta
logger.info("Printing final genome draft.")
genome_size = 0
with open(prefix+'_blob_final_genome.fasta','w') as f:
    for tig in final_genome:
        genome_size+=len(tig.sequence)
        #out= ">" + tig.label + "\n"+ tig.sequence + "\n" ## sequence.outFasta(), fix this
        out = tig.outFasta()
        f.write("{}\n".format(out))
with open(prefix+'_blob_final_nontarget_bin.fasta','w') as f:
    for tig in non_target_bin:
        #genome_size+=len(tig.sequence)
        #out= ">" + tig.label + "\n"+ tig.sequence + "\n" ## sequence.outFasta(), fix this
        out = tig.outFasta()
        f.write("{}\n".format(out))
logger.info("Done. Final draft genome FASTA has been written to "+prefix+'_blob_final_genome.fasta')
logger.info("Draft genome is comprised of "+str(genome_size)+" nucelotides on "+str(len(final_genome))+" contigs.")

## If you've gotten here, then everything should have worked, move everything from temp to ../ and change dir and delete temp
for f in os.listdir('.'):
    os.rename(f,"../"+f)
os.chdir("../")
os.rmdir("temp")



