"""
Created on Tue Nov 14 16:18:37 2017

@author: kevinamses
"""

help_msg = """
 ___   ___  ___  ____  ____  
/ __) / __)/ __)(_  _)(  _ \ 
\__ \( (__( (_-. _)(_  )(_) )
(___/ \___)\___/(____)(____/ 

Author: Kevin Amses (amsesk@umich.edu)
Version: 0.1b

Usage: scgid <module> [args...] --> try scgid <module> -h|--help to see module-specific arguments

Specify a module to continue:
    
    Core Modules
        gc-cov      Generate a draft genome based on gc-coverage-taxonomy plots
        kmers       Generate a draft genome based on an ESOM map
        codons      Generate a draft genome based on relative synonymous codon usage
        consensus   Draw majority rule between three draft genomes into final consensus-derived draft genome
        purify      Filter-out contigs identified by scgid as nontarget in order to reassemble  

    Utilities
        init        Install and configure scgid 
        update      Update scgid from GitHub repo
        buildtaxdb  Build a scgid-compatible taxonomy database from tab-seperated uniprot taxonomy file
        spdb-update    Check for new versions of the swissprot protein database and download if available
        spexpand    Expand your version of the swissprot database to include more proteins with lineage information 

"""


import sys
import subprocess
import os
import inspect
import logging
import logging.config
from scripts.loglib import LoggingEntity, logger_name_gen, ExitOnExceptionHandler
from scripts.gct import Gct 

bin_dir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
pkg_home = os.path.dirname(bin_dir)

logging.config.fileConfig(os.path.join(bin_dir, "logging_config.ini"))
log = logging.getLogger("SCGid")

class SCGid(LoggingEntity, object):
    def __init__(self, call):
        if call == "gct":
            Gct().run()

if len(sys.argv) == 1 or sys.argv[1] in ['-h','--help','-help']:
    print (help_msg)

elif sys.argv[1] == "init":
    arguments = sys.argv[2:]
    call = os.path.join(pkg_home,'bin','init_setup.py')
    arguments.insert(0,call)
    py = sys.executable
    arguments.insert(0,py)
    subprocess.call(arguments)

elif sys.argv[1] == "buildtaxdb":
    arguments = sys.argv[2:]
    call = os.path.join(pkg_home,'bin','build_taxdb.py')
    arguments.insert(0,call)
    py = sys.executable
    arguments.insert(0,py)
    subprocess.call(arguments)
    
elif sys.argv[1] == "spexpand":
    arguments = sys.argv[2:]
    call = os.path.join(pkg_home,'bin','expand_spdb.py')
    arguments.insert(0,call)
    py = sys.executable
    arguments.insert(0,py)
    subprocess.call(arguments)

elif sys.argv[1] == "spdb-update":
    arguments = sys.argv[2:]
    call = os.path.join(pkg_home,'bin','update_swissprot.py')
    arguments.insert(0,call)
    py = sys.executable
    arguments.insert(0,py)
    subprocess.call(arguments)

elif sys.argv[1] == "update":
    arguments = sys.argv[2:]
    call = os.path.join(pkg_home,'bin','update_scgid.py')
    arguments.insert(0,call)
    py = sys.executable
    arguments.insert(0,py)
    subprocess.call(arguments)

elif sys.argv[1] == "kmers":
    arguments = sys.argv[2:]
    call = os.path.join(pkg_home,'bin','esom.py')
    arguments.insert(0,call)
    py = sys.executable
    arguments.insert(0,py)
    subprocess.call(arguments)

elif sys.argv[1] == "gc-cov":
    arguments = sys.argv[2:]
    call = os.path.join(pkg_home,'bin','gc_cov.py')
    arguments.insert(0,call)
    py = sys.executable
    arguments.insert(0,py)
    subprocess.call(arguments)

elif sys.argv[1] == "codons":
    arguments = sys.argv[2:]
    call = os.path.join(bin_dir,'rscu.py')
    arguments.insert(0,call)
    py = sys.executable
    arguments.insert(0,py)
    subprocess.call(arguments)

elif sys.argv[1] == "consensus":
    arguments = sys.argv[2:]
    call = os.path.join(bin_dir,'consensus.py')
    arguments.insert(0,call)
    py = sys.executable
    arguments.insert(0,py)
    subprocess.call(arguments)

elif sys.argv[1] == "purify":
    arguments = sys.argv[2:]
    call = os.path.join(bin_dir,'purify.py')
    arguments.insert(0,call)
    py = sys.executable
    arguments.insert(0,py)
    subprocess.call(arguments)

elif sys.argv[1] == "qualcheck":
    arguments = sys.argv[2:]
    call = os.path.join(bin_dir,'qc.py')
    arguments.insert(0,call)
    py = sys.executable
    arguments.insert(0,py)
    subprocess.call(arguments)

elif sys.argv[1] == "gct":
    log.info("Calling Gct")
    result = SCGid(sys.argv[1])
    log.info("Final message from root.")

else:
    print("\nERROR: Bad module selection. Printing help screen...\n")
    print(help_msg)
