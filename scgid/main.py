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
import itertools
import yaml
import scgid
from scgid.config import InitialConfig, FileConfig
from scgid.modcomm import LoggingEntity, ErrorHandler, logger_name_gen, ExitOnExceptionHandler, get_root, Root
from scgid.gct import Gct 
from scgid.codons import Codons
from scgid.kmers import Kmers
from scgid.consensus import Consensus
from scgid.update_swissprot import SPDBUpdater, SPDBExpander
from scgid.error import ModuleError
from scgid.library import flatten_dict, bcolors
#from scgid.update_scgid import SCGIDUpdate

class SCGidPipeline(object):
    def __init__(self, opts_path):
        self.opts_path = opts_path
        self.root = get_root()
        self.GLOBAL = {}
        self.GCT = {}
        self.KMERS = {}
        self.CODONS = {}

        self.simplelogger = logging.getLogger("SCGid.unfmt")

    def create_options_file(self) -> str:
        gct_args = Gct.generate_argparser()._actions
        kmersT_args = scgid.kmers.Train.generate_argparser()._actions
        kmersA_args = scgid.kmers.Annotate.generate_argparser()._actions
        codons_args = Codons.generate_argparser()._actions

        gct_vars = {o.metavar: o for o in gct_args if o.metavar is not None}
        kmers_vars = {o.metavar: o for o in kmersT_args if o.metavar is not None}
        kmers_vars.update({o.metavar: o for o in kmersA_args if o.metavar is not None})
        codons_vars = {o.metavar: o for o in codons_args if o.metavar is not None}

        #print(dir(gct_vars["assembly_fasta"]))

        self.GCT.update(gct_vars)
        self.KMERS.update(kmers_vars)
        self.CODONS.update(codons_vars)

        all_vars = dict(gct_vars, **kmers_vars)
        all_vars.update(codons_vars)

        for arg in iter( set(itertools.chain(gct_vars, kmers_vars, codons_vars)) ):
            if ( 
                (arg in gct_vars and arg in kmers_vars) or
                (arg in kmers_vars and arg in codons_vars) or
                (arg in gct_vars and arg in codons_vars)
                ):

                self.GLOBAL[arg] = all_vars[arg]

            else:
                pass

        cfg = FileConfig()
        cfg.load_yaml()
        
        for metavar, opt in all_vars.items():
            try:
                dest = opt.option_strings[1].replace("-","")
            except IndexError:
                continue

            if opt.default is None and f"default_{dest}" in cfg.settings:
                self.GLOBAL[metavar].default = cfg.settings[f"default_{dest}"]


        for s, _ in self.GLOBAL.items():
            for modopts in [self.GCT, self.KMERS, self.CODONS]:
                modopts.pop(s, None)


        opts_path = self.opts_path

        with open(opts_path, 'w') as opts_file:
            for mod in ["GLOBAL", "GCT", "KMERS", "CODONS"]:

                opts = getattr(self, mod)

                required = [ f"    {o.metavar}: {o.default}" for k,o in opts.items() if o.required ]
                optional = [ f"    {o.metavar}: {o.default}" for k,o in opts.items() if not o.required ]

                opts_file.write (f"{mod.upper()}:\n")
                if len(required) > 0:
                    opts_file.write ("  REQUIRED:\n")
                    opts_file.write ("\n".join(required))
                    opts_file.write ("\n")
                if len(optional) > 0:
                    opts_file.write ("  OPTIONAL:\n")
                    opts_file.write ("\n".join(optional))
                    opts_file.write ("\n")

        return opts_path


    def read_options_file(self):
        super_config = FileConfig(path=self.opts_path)
        super_config.load_yaml()

        self.GLOBAL = super_config.settings["GLOBAL"]
        self.GCT = super_config.settings["GCT"]
        self.KMERS = super_config.settings["KMERS"]
        self.CODONS = super_config.settings["CODONS"]

        return None

    def show (self):
        self.read_options_file()
        gct_opts = flatten_dict(self.GCT)
        global_opts = flatten_dict(self.GLOBAL)
        
        gct_opts.update(global_opts)

        gct_cmd = Gct(argdict = gct_opts, loglevel=logging.CRITICAL).cli_invocation()

        codons_opts = flatten_dict(self.CODONS)
        codons_opts.update(global_opts)

        codons_cmd = Codons(argdict = codons_opts, loglevel=logging.CRITICAL).cli_invocation()

        kmers_opts = flatten_dict(self.KMERS)
        kmers_opts.update(global_opts)

        kmersT_cmd = scgid.kmers.Train(argdict = kmers_opts, loglevel=logging.CRITICAL).cli_invocation()
        kmersT_cmd.insert(1, "kmers")

        kmersA_cmd = scgid.kmers.Annotate(argdict = kmers_opts, loglevel=logging.CRITICAL).cli_invocation()
        kmersA_cmd.insert(1, "kmers")

        self.simplelogger.info("")
        self.simplelogger.info(f"{' '.join(gct_cmd)}\n")
        self.simplelogger.info(f"{' '.join(codons_cmd)}\n")
        self.simplelogger.info(f"{' '.join(kmersT_cmd)}\n")
        self.simplelogger.info(f"{' '.join(kmersA_cmd)}\n")

        return None


    def run(self):
        self.root.logger.info ("Starting SCGid in pipeline mode (serial). SCGid will run automatically as far as it can.")

        self.root.logger.info (f"Order of module calls will be: {bcolors.MAGENTA}SCGID GCT{bcolors.ENDC} > {bcolors.OKBLUE}SCGID CODONS{bcolors.ENDC} > {bcolors.FAIL}SCGID KMERS TRAIN{bcolors.ENDC} > {bcolors.OKGREEN}SCGID KMERS ANNOTATE{bcolors.ENDC}")

        self.root.logger.info ("Following successful complettion, you need to visually evaluate and carve the ESOM topology. See README for more information.")

        self.root.logger.info ("Fatal errors will halt the entire pipeline, but we'll keep track of where this run dies so you can restart after correcting errors.")

        self.read_options_file()

        self.root.logger.info (f"Using options file at {self.opts_path}. Your options:")

        self.root.simplelogger.info(f"{self.GLOBAL}\n{self.GCT}\n{self.CODONS}\n{self.KMERS}")
        
        # Run GCT with GCT and global options
        gct_opts = flatten_dict(self.GCT)
        global_opts = flatten_dict(self.GLOBAL)
        
        gct_opts.update(global_opts)

        Gct(argdict = gct_opts).run()

        codons_opts = flatten_dict(self.CODONS)
        codons_opts.update(global_opts)

        Codons(argdict = codons_opts).run()

        kmers_opts = flatten_dict(self.KMERS)
        kmers_opts.update(global_opts)

        Kmers(argdict = kmers_opts).run()

        return 0

class SCGid(LoggingEntity, ErrorHandler, Root, object):
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
            #'update': SCGIDUpdate
        }

        self.logger = logging.getLogger("SCGid")
        self.simplelogger = logging.getLogger("SCGid.unfmt")
        
        ''' suspend autoupdate for now - installed SCGid HOME is not a git repo
        # Try to update SCGid from repo in other module calls only if being run in interactive shell
        if self.modcall != "update" and sys.stdout.isatty():
            SCGIDUpdate(is_automated_update=True).run()
        '''
        
        self.logger.info(f"Calling {call}")

        if not os.path.isfile( os.path.join(self.SCGID, "config.yaml") ):
            
            res, ret = self.modcall['init']().run()
        
        else:
            if call in self.modcall:
                self.modcall[call]().run()
            
            else:
                if call == "genopts":
                    try:
                        opts_path = sys.argv[2]
                    except IndexError:
                        opts_path = "scgid.opts"
                    except:
                        return ModuleError()

                    opts_path = SCGidPipeline(opts_path = opts_path).create_options_file()
                    self.logger.info(f"Wrote SCGid options file to `{opts_path}`. Fill it out and invoke by running `scgid run {opts_path}`.")

                elif call == "pipeline":
                    SCGidPipeline(opts_path = sys.argv[2]).run()

                elif call == "show_pipeline":
                    SCGidPipeline(opts_path = sys.argv[2]).show()

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


