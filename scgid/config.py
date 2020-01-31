import os
import yaml
import inspect
import subprocess
import logging
import sys
import time
import threading
import urllib
import urllib3
from importlib_metadata import version
from scgid.modcomm import logger_name_gen, LoggingEntity, ErrorHandler, pkgloc
from scgid.dependencies import Dependencies
from scgid.reuse import ReusableOutputManager
from scgid.error import ConfigError
from scgid.library import report_outcome, is_fasta, file_grep, output_cols, CURSOR_UP_ONE, ERASE_LINE, ow_last_stdout
from scgid.db import UniprotFTP

scgid_banner = """
 ___   ___  ___  ____  ____  
/ __) / __)/ __)(_  _)(  _ \ 
\__ \( (__( (_-. _)(_  )(_) )
(___/ \___)\___/(____)(____/ 
        """
class FileConfig(object):
    def __init__(self, path = None):
        self.SCGID_SCRIPTS = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
        self.SCGID = os.path.dirname(self.SCGID_SCRIPTS)
        self.settings = None
        
        if path is None:
            self.path = os.path.join(self.SCGID, "config.yaml")
        else:
            self.path = path

    def load_yaml(self) -> None:
        if os.path.isfile( f"{self.path}.local" ):
            loc = f"{self.path}.local"

        else:
            loc = self.path

        try:
            with open(loc, 'r') as cfg:

                self.settings = yaml.load(cfg, Loader=yaml.BaseLoader)
                return None

        except IOError:
            return ConfigError("Unable to locate configuration file")
        except:
            return ConfigError("Malformed configuration file")

    def write_yaml(self, path = None) -> None:
        if path is None:
            path = self.path
        with open(path, 'w') as cfg:
            yaml.dump(self.settings, cfg)
        return None


class Config(LoggingEntity, ErrorHandler):
    def __init__(self):
        self.logger = logging.getLogger( logger_name_gen() )
        self.SCGID_SCRIPTS = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
        self.SCGID = os.path.dirname(self.SCGID_SCRIPTS)
        self.OUTPUTSUFFIX = "_scgid_output"
        self.dependencies = Dependencies()
        self.reusable = ReusableOutputManager()
        self.rundir = None

    def __repr__(self):
        return "\n".join(["{}: {}".format(key, setting) for key,setting in self.__dict__.items()])

    def get(self, value):
        try: 
            attr = getattr(self, value)
            return attr

        except KeyError:
            raise KeyError(f"Member `{value}` not in Config")
        

    def load_cmdline(self, parsed_args, case_args={}):
        for unset in [x for x,p in parsed_args.__dict__.items() if p is None]:
            try:
                setattr(parsed_args, unset, getattr(self, f"default_{unset}"))
            except:
                self.logger.debug(f"No loaded value from config.yaml to set None-type argument `{unset}`")
        
        self.__dict__.update( vars(parsed_args) )
    
    def load_argdict(self, argdict):
        self.__dict__.update(argdict)

    def load_yaml(self):
        #loc = os.path.join(pkg_resources.resource_string(__name__, "config.yaml"))
        loc = os.path.join(self.SCGID, "config.yaml")
        if os.path.isfile( "{}.local".format(loc) ):
            loc = "{}.local".format(loc)
        try:
            with open(loc, 'r') as cfg:
                self.logger.info (f"Using package configuration file located at `{loc}`")
                self.__dict__.update( yaml.load(cfg, Loader=yaml.BaseLoader) )
        except IOError:
            return ConfigError("Unable to locate configuration file")
        except:
            return ConfigError("Malformed configuration file")
    
    def check_case_args (self, mapper_dict):
        for option_arg, options in mapper_dict.items():
            choice = options[self.get(option_arg)]
            if not choice["bool"]:
                self.logger.warning(choice["warning"])

    def check_esom_path (self):
        # Check to see if esomtrn exists in configured `esom_path`
        if not os.path.isfile(
            os.path.join(
                self.get("esom_path"),
                "bin",
                "esomtrn"
            )):
            self.logger.warning(f"Invalid path to ESOM specified in config.yaml: `{self.get('esom_path')}`")

            # If not, try to find it in environment and update module config
            p = subprocess.Popen( ["which", "esomtrn"], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
            path, _ = p.communicate()
            if p.returncode:
                return ConfigError("ESOM path either unset or invalide. Checked in $PATH and config.yaml")

            else:
                path = os.path.dirname(path).decode("utf-8")
            
                # Updating moduleconfig here
                setattr(self, "esom_path", path)

class InitialConfig(object):
    def __init__(self):
        self.HOME, self.SCRIPTS = pkgloc()
        self.VERSION = "1.0.0"
        self.DB = os.path.join(self.HOME, 'database')
        self.initial_config = {}
        self.initial_config["mpicmd"] = "mpirun"

    def welcome(self):
        print (scgid_banner)
        print (f"----- Welcome to SCGid, version {version('SCGid')} -----")
        print ("> SCGid is provided under the GNU General Public License v3.0\n")
        print ("> This script is going to help you configure SCGid by defining some package-wide variables and downloading databases.\n")
        if os.getcwd() != self.HOME:
            os.chdir(self.HOME)
            print (f"> Navigating to scgid bin directory at {self.HOME} ...\n")

        if sys.version[0] != '3':
            print (f"{output_cols['RED']}ERROR: SCGid v1.0.0 is only compatible with python 3.x.x. Exitting...{output_cols['RESET']}")
            sys.exit(1)

    def config_third_party_deps(self):
        tpdeps = {
            "esom_path": {
                'target': lambda path: os.path.isfile(os.path.join(path,"bin","esomstart")),
                'question': "Enter full path to ESOM install: "
                },
            "clams_path": {
                'target': lambda path: os.path.isfile(path),
                'question': "Enter full path to ClaMS-CLI.jar: "
                },
            "path_to_Rscript": {
                'target': lambda path: os.path.isfile(path),
                'question': "Enter full path to Rscript: "
                }
            }
        for var, stuff in tpdeps.items():
            while True:
                entry = input(stuff['question'])
                if entry == "DEBUG_SKIP":
                    break
                else:
                    if not stuff["target"](entry):
                        report_outcome (stuff['question'], "RED", "[NOT FOUND, TRY AGAIN]")
                    else:
                        report_outcome (stuff['question'], "GREEN", "[GOOD]")
                        tpdeps[var] = entry
                        break

        self.initial_config.update(tpdeps)
    
    def config_spdb (self):

        self.initial_config['default_spdb'] = None
        self.initial_config['spdb_version'] = None

        print ("\n> SCGid uses the UniProt swissprot database as the baseline database for blastp annotations.")
        while True:
            entry = input("Download the current version of the swissprot database? [y/n] ").strip().lower()
            if entry not in ['y','n']:
                continue
            else:
                if entry == 'y': 
                    ## download database and update settings with version and path (both handled by update_swissprot.py
                    entry = input(f"Specify protein database download directory: [default = '{self.DB}'] ")
                    
                    if len(entry) != 0:
                        self.DB = entry
                    else:
                        pass
                    
                    try:
                        os.makedirs(self.DB)
                    
                    except FileExistsError:
                        pass
                    
                    with UniprotFTP() as uniprot:
                        path = uniprot.retr_spdb(self.DB)
                        path_to_spdb = uniprot.decompress(path)
                        
                        print(f"SPDB @ {path_to_spdb}")

                    self.initial_config['default_spdb'] = path_to_spdb
                    self.initial_config['spdb_version'] = uniprot.remote_reldate

                    break
                else: #entry.lower() is 'n'
                    while True:
                        question = "Path to your database (FASTA): "
                        u_db = input(question)
                        
                        ## Verify the existence and swissprot-like formatting of user-supplied database ##
                        if not os.path.isfile(u_db):
                            report_outcome(question, "RED", "[NOT FOUND, TRY AGAIN]")
                        else:
                            if not is_fasta(u_db):
                                report_outcome(question, "RED", "[SPECIFIED DB NOT IN FASTA FORMAT]")
                            else:
                                all_with_OS = True
                                for header in file_grep ("^>", u_db, 'multiple'):
                                    if "OS=" not in header:
                                        report_outcome(question, "RED", "[HEADERS NOT SWISSPROT-STYLE, see README]")
                                        all_with_OS = False
                                        break
                                if not all_with_OS:
                                    continue
                                else:
                                    report_outcome(question, "GREEN", "[GOOD]")
                                    break
                    ### add location of user-defined swissprot-style database to settings
                    self.initial_config['default_spdb'] = u_db
                    break
    
    def config_taxdb (self):
        #%% Figure out the swissprot taxonomy database ##
        self.initial_config['taxonomy_all_tab'] = None
        self.initial_config['default_taxdb'] = None
        "/Users/aimzez/dev/SCGid/database/uniprot-taxonomy-all-120919.tab"

        dwnld_link = "https://www.uniprot.org/taxonomy/?query=*&format=tab"
        date = time.strftime("%m%d%y")

        print ("\n> SCGid requires a copy of the swissprot taxonomy database.")
        entry = input("Download the current version of the swissprot taxonomy database? This can take awhile. [y/n] ").strip().lower()
        while True:    
            if entry not in ['y','n']:
                continue
            else:
                if entry == 'y':
                    self.initial_config['taxonomy_all_tab'] = os.path.join(self.DB, f"uniprot-taxonomy-all-{date}.tab")
                    dwnld_to = self.initial_config['taxonomy_all_tab']
                    if os.path.isfile(dwnld_to):
                        os.remove(dwnld_to)
                    print(f"\n> Downloading uniprot taxonomy database to: {dwnld_to}...\n")
                    
                    http = urllib3.PoolManager()
                    r = http.request('GET', dwnld_link, preload_content=False)
                    with open(dwnld_to, 'wb') as taxall:
                        while True:
                            chunk = r.read(100000)
                            if not chunk:
                                break
                            taxall.write(chunk)
                            csize = os.path.getsize(dwnld_to)
                            print(CURSOR_UP_ONE,ERASE_LINE,csize)
                    break 
                else:
                    self.initial_config['taxonomy_all_tab'] = None
                    print ("> Okay, but you'll have to download and build the taxonomy database manually before using SCGid. See README for more information.")
                    break

        #%% Build taxonomy database
        if self.initial_config['taxonomy_all_tab'] is not None:
            try_build = subprocess.call([
                sys.executable, os.path.join(self.SCRIPTS,'build_taxdb.py'), 
                '-db', self.initial_config['default_spdb'], 
                '-t', self.initial_config['taxonomy_all_tab'], 
                '-o', f"{self.initial_config['default_spdb']}.taxdb"
                ])
            if try_build: ## STOP if buildtaxdb exit code != 0
                print (f"> {output_cols['RED']}[ERROR]{output_cols['RESET']} Fatal error encountered while building taxonomy database.")
                sys.exit(1)

            self.initial_config['default_taxdb'] = f"{self.initial_config['default_spdb']}.taxdb"

        else:
            self.initial_config['default_taxdb'] = None

    def write_config(self):
        cfg_path = os.path.join(self.HOME, "config.yaml")
        print (f"\n> Writing your settings to {cfg_path}")
        print (f"\nSCGid Config\n{'-'*24}")
        with open(cfg_path, 'w') as cfg:
            yaml.dump(self.initial_config, cfg)
            print('\n'.join([f"{k}: {v}" for k,v in self.initial_config.items()]))
        print ("-"*24+"\n")

    def run(self):
        
        self.welcome()

        self.config_third_party_deps()

        self.config_spdb()

        self.config_taxdb()

        self.write_config()

