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
Version: 1.0.0

Usage: scgid <module> [args...] --> try scgid <module> -h|--help to see module-specific arguments

Specify a module to continue:
    
    Core Modules
        gct         Generate a draft genome based on gc-coverage-taxonomy plots
        kmers       Generate a draft genome based on an ESOM map
        codons      Generate a draft genome based on relative synonymous codon usage
        consensus   Draw majority rule between three draft genomes into final consensus-derived draft genome
        purify      Filter-out contigs identified by scgid as nontarget in order to reassemble  

    Utilities
        init        Install and configure scgid 
        update      Update SCGid from GitHub repo
        spdbup      Check for new versions of the swissprot protein database and download if available
        spexpand    Expand your version of the swissprot database to include more proteins with lineage information 

"""


import sys
import subprocess
import os
import inspect
import logging
import logging.config
import warnings
from scgid.config import InitialConfig
from scgid.modcomm import LoggingEntity, logger_name_gen, ExitOnExceptionHandler
from scgid.gct import Gct 
from scgid.codons import Codons
from scgid.kmers import Kmers
from scgid.consensus import Consensus
from scgid.update_swissprot import SPDBUpdater, SPDBExpander
from scgid.update_scgid import SCGIDUpdate

class SCGid(LoggingEntity, object):
    def __init__(self, call):

        # Show all warnings
        warnings.simplefilter("default")
        
        self.SCGID_SCRIPTS = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
        self.SCGID = os.path.dirname(self.SCGID_SCRIPTS)

        logging.config.fileConfig(os.path.join(self.SCGID_SCRIPTS, "logging_config.ini"))

        self.modcall = {
            'init': InitialConfig,
            'gct': Gct,
            'codons': Codons,
            'kmers': Kmers,
            'consensus': Consensus,
            'spdbup': SPDBUpdater,
            'spexpand': SPDBExpander,
            'update': SCGIDUpdate
        }

        self.logger = logging.getLogger("SCGid")
        
        # Try to update SCGid from repo in other module calls only if being run in interactive shell
        if self.modcall != "update" and sys.stdout.isatty():
            SCGIDUpdate(is_automated_update=True).run()
        
        self.logger.info(f"Calling {call}")

        if not os.path.isfile( os.path.join(self.SCGID, "config.yaml") ):
            
            self.modcall['init']().run()
        
        else:
            if call in self.modcall:
                self.modcall[call]().run()
            
            else:
                self.logger.critical(f"Bad module selection `{call}`")
                print(help_msg)


        self.logger.info("Final message from root.")

def main():
    if len(sys.argv) == 1 or sys.argv[1] in ['-h','--help','-help']:
        print (help_msg)
        sys.exit(0)

    else:
        SCGid(sys.argv[1])

'''

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

else:
    print("\nERROR: Bad module selection. Printing help screen...\n")
    print(help_msg)
'''
