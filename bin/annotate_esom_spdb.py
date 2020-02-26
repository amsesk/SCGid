import sys
import os
from lib import best_blast_hit, SPDB_OS_REGEXP_PATTERN, parse_spdb_blastout
import sys
import re
import ast
import os
import argparse


def generate_annotation_file(taxrpt, scheme):
    annotation_file_towrite = ""
    pair_spl = scheme.split('/')
    arch = {} #Architecture of rules for target-except mode
    pair_num = 1
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
        arch[pair_num] = {'t': tar, 'e': exceptions}
        pair_num += 1

    for key, val in taxrpt.iteritems():
        for pair_id,te_pair in arch.iteritems():
            pair_cls = pair_id 
            if te_pair['e'] is None:
                if te_pair['t'] in val:
                    annotation_file_towrite += key+'\t'+str(pair_cls)+'\n'
            else:
                if te_pair['t'] in val and not any(e in te_pair['e'] for e in val):
                    annotation_file_towrite += key+'\t'+str(pair_cls)+'\n'

    annot_output = os.path.join(os.getcwd(),'annotation_temp')
    with open(annot_output,'w') as annot:
        annot.write(annotation_file_towrite)

    lines_in_annot = annotation_file_towrite.split('\n')
    lines_in_annot_count = len(lines_in_annot)-1

    in_annot = [x.split('\t')[0] for x in lines_in_annot]
    for key,val in taxrpt.iteritems():
        if key not in in_annot:
            pass
            #print key+":"+str(val)
    if lines_in_annot_count < len(taxrpt):
        print "Some blastn-annotated contigs will be makred as 'Unclassified' because your target-exception groups are not completely inclusive. Did you mean to do this?" 
    elif lines_in_annot_count > len(taxrpt):
        print "There is overlap between your target-exception groups. This is going to cause ambiguity based on execution order. You should fix this by making your target-exception pairs MUTUALLY EXCLUSIVE. See documentation for more detailed information."

    print "ESOM target_exception-based annotation file written to "+annot_output

def generate_taxrpt_from_raw_blastout(blastout, spdb, taxdb):
    output = []
    blastout_bestout = "{}.best".format(sys.argv[1])
    for query, hit in best_blast_hit( blastout, bitcol = 11).iteritems():
        output.append('\t'.join(hit))

    output='\n'.join(output)
    with open(blastout_bestout,'w') as f:
        f.write(output)

    parsed_blastout_out = "{}.best.parsed".format(sys.argv[1])
    parsed_blastout = parse_spdb_blastout(sys.argv[2], blastout_bestout, None)
    with open(parsed_blastout_out, 'w') as f:
        f.write('\n'.join(parsed_blastout))

    #set blastout to the parsed,best blastout for further processing
    parsed_blastout = os.path.abspath(parsed_blastout_out)

    ## (2nd) Grab hit info from the parsed blast output file and use OS=<blah> to grab lineage information from the taxdb
    pattern = re.compile(SPDB_OS_REGEXP_PATTERN)

    ## This should be set somewhere else
    taxlvl_idx = 1

    taxrpt = {}
    for line in open(parsed_blastout).readlines():
        spl = [x.strip() for x in line.split('\t')]
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
        except:
            try:
                # Can get rid of these once local taxdb is rebuilt
                sp_os = sp_os.replace("'","")
                sp_os = sp_os.replace("#","")
                lineage = taxdb[sp_os]
            except KeyError:
                logger.warning("'"+sp_os+"' is missing from the taxdb. You probably specified a -t|--taxdb file built from a likely out-of-date swissprot database other than the one you blasted against. Rerun build_taxdb or specify a different taxdb to correct. Lineage set as 'Not_in_taxdb' for this run.")
                lineage = "Not_in_taxdb"

        taxrpt[query] = [x.strip() for x in lineage.split(";")]

    return taxrpt

blastout = sys.argv[1]
spdb = sys.argv[2]
taxdb = ast.literal_eval(open(sys.argv[3]).read())
scheme = sys.argv[4]

taxrpt = generate_taxrpt_from_raw_blastout(blastout, spdb, taxdb)
generate_annotation_file(taxrpt, scheme)