
#%% IMPORTS
import sys
import re
import os
import ast
from sequence import *
import subprocess
import threading
import argparse
import inspect
import settings
from lib import file_yield_lines, spdb_grab_os, alltaxtab2dict, do_long_wait, output_cols

#handle arguments
parser = argparse.ArgumentParser()
parser.add_argument('-db','--spdb', metavar="spdb",action="store",default=settings.path_to_spdb, required=True,help = "A file in FASTA format containing the swissprot database to be used. FASTA headers must be in standard swissprot format [... OS= ... GN= | PE=). See README for more information and examples.")
parser.add_argument('-t','--taxdb', metavar="taxonomy_db",action="store", required=True,help = "The uniprot taxonomy database downloaded in tsv format. See README for more information.")
parser.add_argument('-o','--output', metavar="output_path",action="store",required=True,help = "The destination and filename of the taxonomy database to be created.")
parser.add_argument('-q', '--quiet', action="store_true", required=False, default=False, help="In quiet mode, no status messages will be passed to stdout.")
parser.add_argument('-a', '--append', metavar="tax_to_add",action="store", required=False, help="Two-column tab-separated file containing the species names (column 1) and semi-colon-separated lineage information (column 2) for proteins you added to the input swissprot fasta. Note that that fasta labels must be in swissprot format ... see README.")
args =  parser.parse_args()

spdb = os.path.abspath(args.spdb)
output_file = os.path.abspath(args.output)
alltaxtab = os.path.abspath(args.taxdb)

for i in [spdb, alltaxtab]:
    if not os.path.isfile(i):
        raise IOError("File not found: %s" % (i))

bin_dir =os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
pkg_home = os.path.dirname(bin_dir)

if not args.quiet:
    print "PATH_TO_SPDB: %s" % spdb
    print "PATH_TO_TAXONOMY_ALL: %s" % alltaxtab
    print "> Starting build of taxonomy database"
    print "> Scanning database for USIDs (note, OS=<USID>)"
#%%
if not args.quiet:
    unique_sp_os = do_long_wait(lambda: spdb_grab_os(file_yield_lines(spdb)), 'list')
else:
    unique_sp_os = spdb_grab_os(file_yield_lines(spdb))

# flatten nested list that gets returned 
unique_sp_os = [item for sub in unique_sp_os for item in sub] 

''' both replaced by spdb_grab_os
all_sp_os = []
for line, sp_os in copy_spdb_grab_os(file_yield_lines(spdb)):
    if sp_os is not None:
        all_sp_os.append(sp_os)
unique_sp_os = list(set(all_sp_os))
'''
''' replaced by above code block
command = ["sh",os.path.join(bin_dir,"extract_os.sh"),spdb]
p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
out,err = p.communicate()
unique_sp_os = open("unique_os").readlines() #from extract_os.sh
unique_sp_os = map(str.strip, unique_sp_os)
os.remove("unique_os")
'''
#%%
if not args.quiet:
    print "> There are [%d] UISDs present in your version of the swissprot database." %(len(unique_sp_os))
    print "> Extracting lineage information for these USIDs from all-taxaonomy.tab file"

    db = do_long_wait(lambda: alltaxtab2dict(file_yield_lines(alltaxtab)), 'dict')
else:
    db = alltaxtab2dict(file_yield_lines(alltaxtab))

#%% Try to find links between unique_sp_os and all-taxonomy.tab, assign them
taxdb_new = {}
for line in unique_sp_os:
    sp_os = line.strip()
    sp_os = sp_os.replace("'","")
    orig_os = sp_os # save original sp_os for taxdb_new key later on
    try:
        lineage = db[sp_os]
        #print sp_os
    except:
        sp_os_spl = sp_os.split(" ") # split the damn thing once
        #print sp_os_spl
        if len(sp_os_spl) == 1: #if sp_os is one word, use that as database key
            #sp_os = sp_os
            pass
        elif len(sp_os_spl) == 2: #if sp_os is two words, use first word as database key
            sp_os = sp_os_spl[0]
        else:
            s = re.search("(^[^(]+)",sp_os) # capture everything up until a '(' to eliminate strain info, etc...
            if s is None: #If no match, get rid of last word???
                split = sp_os_spl(" ")
                sp_os = " ".join(split[0:len(split)-1]) ## Why did I want to get rid of the last word??? Typo? Leaving it for now...must have had to do with something? Maybe, maybe not. Only thing that would be None from re is an sp_os that started with a '('
            else:
                sp_os = s.group(1).strip()
        try:
            lineage = db[sp_os]
        except:
            try:
                lineage = next(lin for key,lin in db.iteritems() if sp_os.strip() in key) ## look for the sp_os in other database keys
            except:
                if len(sp_os.split(" ")) == 2: ## if sp_os has two words, look for the first one (hoping for genus) in other database keys, and if found, assign that lineage information to the key
                    try:
                        lineage = next(lin for key,lin in db.iteritems() if sp_os.split(" ")[0].strip() in key) 
                    except:
                        lineage = "Not in taxdb" ## if the first word still isn't found, apply 'Not in taxdb' to the key and give up
                else:
                    lineage = "Not in taxdb" ## if sp_os is only one word, that we already couldn't find in other keys, assign 'Not in taxdb' and give up
    taxdb_new[orig_os] = lineage


## Write taxonomy database in python dictionary format to output_file
with open(output_file, 'w') as t:
    i=1
    t.write("{\n")
    for sp_os, lineage in taxdb_new.iteritems():
        t.write("'%s': '%s'" % (sp_os, lineage))
        if i != len(taxdb_new):
            t.write(",\n")
        else:
            t.write("\n}")
        i+=1

#Check dictionary length and quality
new_dict = ast.literal_eval(open(output_file).read())
if args.quiet == False:
    if len(new_dict.keys()) == len(unique_sp_os):
        print "> Checking the size of the new taxonomy database... %s%s%s" % (output_cols['GREEN'],"[GOOD]",output_cols['RESET'])
    else:
        print "> Something went wrong. Length of new taxbd is not equal to that of list of unique taxon identifiers %d != %d %s%s%s" % (len(new_dict.keys()), len(unique_sp_os), output_cols['RED'],"[FAILURE]",output_cols['RESET'])

    nolink_count = new_dict.values().count("Not in taxdb")
    if nolink_count == 0:
        report = output_cols['GREEN']+"[GOOD]"+output_cols['RESET']
    else:
        output_cols['RED'],"[GOOD]",output_cols['RESET']
    print "> There are [%d] USIDs in your sequence database without a link in the new taxonomy database... %s" % (nolink_count, report)
    print "> You're taxonomy database has been written to: "+os.path.join(os.getcwd(),output_file)+"."
