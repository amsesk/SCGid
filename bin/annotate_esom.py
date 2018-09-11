#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 20 14:13:57 2017

@author: kevinamses
"""

import argparse
import os
import subprocess
#from ete3 import NCBITaxa
from lib import *
import settings


#%%
parser = argparse.ArgumentParser()

parser.add_argument('-n','--nucl', metavar = "contig_fasta", action="store",required=True, help = "A FASTA file containing the nucleotide assembly. (MANDATORY)")
parser.add_argument('-f','--prefix', metavar = 'prefix_for_output', required=False, default='scgid', help="The prefix that you would like to be used for all output files. DEFAULT = scgid")
parser.add_argument('-b','--blastout', metavar = "blastout", action="store",required=False, help = "The blast output file from a blastn search of the NCBI nt database with your contigs as query. If you have not done this yet, this script will do it for you.")

parser.add_argument('-m','--mintig', metavar = "minimum_contig_size", action="store",required=False, default="1000", help = "Contig size cutoff (in nucleotides) for inclusion in the ESOM training. Default = 1000 bp")
parser.add_argument('-w','--window', metavar = "window_size", action="store",required=False, default="1000", help = "Size of the window in which kmer frequencies are calculated. Default = 1000 bp")

parser.add_argument('--cpus', metavar = 'cores', action = 'store', required = False, default = "1", help = "The number of cores available for BLAST to use.")
parser.add_argument('-e', '--evalue', metavar = 'e-value_cutoff', action = 'store', required = False, default = '1e-5', help = "The evalue cutoff for blast. Default: 1xe-5)")
#parser.add_argument('-rm', '--rankmode', action = 'store_true', required = False, help = "Annotate contigs at the same single taxonomic rank across all contigs. (e.g. superkingdom)")
#parser.add_argument('-te', '--targetexcept', action = 'store_true', required = False, help = "Annotate contigs at a varietry of taxonomic ranks across all contigs. (e.g. Eukaryota, except Fungi). Must be used in combination with -s|--annotation_scheme.")
parser.add_argument('-s','--annotation_scheme', metavar = 'annotation_scheme', action = 'store', required = True, help = "The annotation scheme to use in target_except annotation mode (RECOMMENDED). Groups MUST be mutually exclusive to avoid overlap and unincluded groups will be arbitrarily marked as 'Unclassified'. You ABSOLUTELY MUST use this basic syntax: '-s|--annotation_schem target1^exception1,excetions2/target2^exception1/target3/etc...'. Example: '-s|--annotation_scheme Eukaryota^Fungi,Metazoa/Fungi/Metazoa/Bacteria^Proteobacteria/Proteobacteria'. See documention for detailed information and more examples.")
parser.add_argument('-i','--infotable', metavar = "infotable", action="store",required=False, help = "The scgid gc-cov-derived infotable generated by a blastp search of a swissprot-style protein database.")
#parser.add_argument('--mode', metavar = "mode", action="store",required=False, default ='blastn', help = "The type of blast results that you would like to use to topology ('blastp' or 'blastn'). This module will automatically do a blastn search of the NCBI nt database for you. At this time, a blastp search can not be run directly from this script. INSTEAD, if using mode 'blastp' (DEFAULT, recommended) you must specify a scgid blob-derived <prefix>_info_table.tsv file with -i|--infotable")

args =  parser.parse_args()

#%%
prefix = args.prefix
nucl_path = os.path.abspath(args.nucl)
file_prefix = os.path.split(nucl_path)[1]
bin_dir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
pkg_home = os.path.dirname(bin_dir)
this_module = "kmers"

blastout = args.blastout
info_table = args.infotable

## only target-except mode supported at this time, so set it to True
targetexcept = True

### mode, 'blastp' not yet supported for `scgid kmers annotate`
#mode = args.mode
mode = 'blastn'

if blastout is not None:
    blastout = os.path.abspath(blastout)
if info_table is not None:
    info_table = os.path.abspath(info_table)

blastout_check = True #use this to decide whether or not blastn has to be run or an IOError has to be raised 

#%% navigate to head of working directory
try:
    os.chdir(args.prefix+'_scgid_output')
except:
    os.mkdir(args.prefix+'_scgid_output')
    os.chdir(args.prefix+'_scgid_output')

#%%

logger = start_logging(this_module, sys.argv)

#%% see if nt.blast.out has been done for rscu or is present here in the esom folder under its default name
if mode == 'blastn':
    if blastout is None:
        #try to locate nt blastout in default locations
        blastout = os.path.join(os.getcwd(),'esom',prefix+".nt.blast.out") #check here first
        if not os.path.isfile(blastout):
            blastout = os.path.join(os.getcwd(),'rscu',prefix+".nt.blast.out") #check in rscu next
            if not os.path.isfile(blastout):
                logger.info("Nothing given for -b|--blastout and unable to locate nt blast output file in default locations. BLASTN search will have to be run prior to annotation.")
                blastout_check = False #Will have to run BLASTN below
                blastout = os.path.join(os.getcwd(),'esom',prefix+".nt.blast.out") #set blastout back to default
    else:
        if not os.path.isfile(blastout):
            #have to raise IOError since the user provided a nonexistent file, blastout_check = False, in spirit but unecessary
            logger.critical("-b|--blastout: No such file or directory, "+blastout)
            raise IOError
        else:
            pass #blastout_check set to True by default

#%% navigate to esom output directory
try:
    os.chdir('esom')
except:
    os.mkdir('esom')
    os.chdir('esom')

#%% BLASTN, if necessary ###
if mode == 'blastn': 
    if not blastout_check:

        #Make sure that BLAST is available
        try:
            subprocessT(['blastn','-help'])
        except:
            logger.exception("BLAST is either not installed or unavailable.")
            raise IOError("BLAST is either not installed or unavailable.")

        logger.info("No blastn output file detected. Running blastn on contig file.")

        outfmt = "6 qseqid sseqid pident qlen length mismatch gapope evalue bitscore staxids"
        blastn_cmd = ["blastn","-query",nucl_path,"-max_target_seqs","1","-num_threads",args.cpus,"-db","nt","-outfmt", outfmt, "-evalue",args.evalue,"-out",prefix+".nt.blast.out"]
        logger.info(' '.join(blastn_cmd))
        p = subprocessP(blastn_cmd, logger)

    for key,best in best_blast_hit(blastout).iteritems():
        with open(blastout+".best",'a') as f:
            f.write('\t'.join(best))
            f.write('\n')
        with open(blastout+".best.taxids",'a') as f:
            f.write(best[0]+"\t"+best[8]+"\n")
    
    logger.info("Pulling taxonomic information for nucleotide scaffolds from ncbi taxonomy database...")
    ncbi = NCBITaxa()
    ids = {}
    for i in open(blastout+".best.taxids").readlines():
        cols = i.split('\t')
        name = cols[0].strip()
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
else:
    info_table = pd.read_csv(info_table, sep="\t", header = None)
    info_table.columns = ["contig","prot_len","coverage","gc","pid","sp_os","lineage","evalue","parse_lineage"]
    

##### ANNOTATION #####
##using target-except mode
        
if targetexcept is True:
    annotation_file_towrite = ""
    scheme = args.annotation_scheme
    pair_spl = scheme.split('/')
    arch = {} #Architecture of rules for target-except mode
    pair_num = 0
    cls_legend = {}

    ## parse rule pair
    for pair in pair_spl:
        cls_legend[pair_num]=[pair]
        tar = None
        exceptions = None
        tar = pair.split('^')[0]
        try: #try to see if there are exceptions to target
            exceptions = pair.split('^')[1]
            exceptions = exceptions.split(',')
        except:
            pass #No exceptions to target, continue on with exceptions = None
        ## append rule pair to architecture
        arch['p_'+str(pair_num)] = {'t': tar, 'e': exceptions}
        pair_num += 1
    for key, val in ids.iteritems():
        #print key, val
        for pair_id,te_pair in arch.iteritems():
            pair_cls = int(pair_id.split('_')[1]) + 1 ## +1 to allow Unclassifieds to be cls 0
            if te_pair['e'] is None:
                if te_pair['t'] in val.values():
                    annotation_file_towrite += key+'\t'+str(pair_cls)+'\n'
            else:
                if te_pair['t'] in val.values() and not any(e in te_pair['e'] for e in val.values()):
                    annotation_file_towrite += key+'\t'+str(pair_cls)+'\n'

    #annot_output = os.path.join(os.getcwd(),file_prefix+'.annotation')
    annot_output = os.path.join(os.getcwd(),'annotation_temp')
    with open(annot_output,'w') as annot:
        annot.write(annotation_file_towrite)
    lines_in_annot = annotation_file_towrite.split('\n')
    lines_in_annot_count = len(lines_in_annot)-1
    print lines_in_annot_count
    print len(ids)
    in_annot = [x.split('\t')[0] for x in lines_in_annot]
    for key,val in ids.iteritems():
        if key not in in_annot:
            print key+":"+str(val)
    if lines_in_annot_count < len(ids):
        logger.warning("Some blastn-annotated contigs will be makred as 'Unclassified' because your target-exception groups are not completely inclusive. Did you mean to do this?" )
    elif lines_in_annot_count > len(ids):
        logger.critical("There is overlap between your target-exception groups. This is going to cause ambiguity based on execution order. You should fix this by making your target-exception pairs MUTUALLY EXCLUSIVE. See documentation for more detailed information.")

logger.info("ESOM target_exception-based annotation file written to "+os.path.join(os.getcwd(),nucl_path+'.annotation'))

## print tetramer freqs again to generate annotated .cls file for viewing #
if not os.path.isfile(file_prefix):
    arguments = ['cp',nucl_path,file_prefix]
    subprocessP(arguments)

if os.path.isfile(os.path.join(os.getcwd(),file_prefix+'.annotation')):
    os.remove(os.path.join(os.getcwd(),file_prefix+'.annotation'))

call = os.path.join(bin_dir,'print_tetramer_freqs_deg_filter_esom_VD.pl')
arguments = ['perl',call,'-s',file_prefix,'-m',args.mintig,'-w',args.window,'-k','4','-a',annot_output]
logger.info(' '.join(arguments))
subprocessC(arguments)
#os.remove(os.path.join(os.getcwd(),'annotation_temp'))

## Color-coding and naming for classes ##

to_add="%0\tUnclassified (0)\t255\t255\t255\n"
for pair_num, val in cls_legend.iteritems():
    print pair_num+1
    to_add+="%"+str(pair_num+1)+"\t"+val[0]+"("+str(pair_num+1)+")"+"\t"+"\t".join(map(str,random_color()))+"\n"

#read pertinent things from the file
with open(file_prefix+'.cls','r+') as cls:
    cls.seek(0,0)
    head = cls.readline()
    data = cls.read()

#reprint with legend
with open(file_prefix+'.cls','w') as cls:
    cls.write(head)
    cls.write(to_add)
    cls.write(data)

logger.info("Annotated class file written to "+file_prefix+'.cls')

'''
##### ANNOTATION #####
### If using rank mode
if args.rankmode is True:
    for key, val in ids.iteritems():
        #print key+"\t"+str(val)
        ranks = ncbi.get_rank(val.keys())
        rank_of_interest = [i for i in ranks if ranks[i] == 'superkingdom']
        #print key, rank_of_interest, ncbi.get_rank(rank_of_interest), ncbi.get_taxid_translator(rank_of_interest).values()
        with open(nucl_path+'.tax','a') as tax:
            l = key+'\t'+ncbi.get_taxid_translator(rank_of_interest).values()[0]+'\n'
            tax.write(l)
    #if len(rank_of_interest) == 0:
    #    print val
    #    print ranks
    ## Figure out how many classes will be needed. That is, how many unique taxons are there?
    with open(nucl_path+'.tax','r') as tax:
        uniq = []
        for line in tax.readlines():
            name = line.split('\t')[1].strip()
            if name not in uniq:
                uniq.append(name)
    ## Rewrite .tax file with class numbers instead of names
    with open(nucl_path+'.tax','r') as tax:
        for line in tax.readlines():
            node = line.split('\t')[0].strip()
            name = line.split('\t')[1].strip()
            with open(nucl_path+'.annotation','a') as annot:
                cls = uniq.index(name) + 1
                to_write = node+'\t'+str(cls)+'\n'
                annot.write(to_write)
    logger.info("ESOM rankmode-based annotation file written to "+os.path.join(os.getcwd(),nucl_path+'.annotation'))
'''
