import os
import time
import yaml
import sys
import subprocess
import argparse
import ast
from scgid.modcomm import pkgloc, Head, ErrorHandler, LoggingEntity, logger_name_gen
from scgid.parsers import PathAction, OutputPathStore
from scgid.module import Module
from scgid.db import UniprotFTP
from scgid.library import output_cols
from scgid.error import ModuleError
from scgid.parsers import SPDBTaxonomy
from scgid.config import FileConfig

class SPDBUpdater:
    def __init__(self):
        self.HOME, self.SCRIPTS = pkgloc()
        self.config_path = os.path.join(self.HOME, "config.yaml")

    def update_config(self, yamlcfg, new_spdb, new_taxdb, new_version):

        yamlcfg["default_spdb"] = new_spdb
        yamlcfg["default_taxdb"] = new_taxdb
        yamlcfg["spdb_version"] = new_version

        with open(self.config_path, 'w') as cfg:
            yaml.dump(yamlcfg, cfg)
        
        return None
        
    def run(self):
        
        with UniprotFTP() as uniprot:
            if not uniprot.needs_retr (
                os.path.join(self.HOME, "config.yaml"),
                uniprot.remote_reldate
            ):
                
                print (f"Nothing to do. SPDB alrady up to date.")

            else:
                with open(self.config_path, 'r') as cfg:
                    yamlcfg = yaml.load(cfg, Loader=yaml.BaseLoader)

                dest_dir = os.path.dirname(yamlcfg["default_spdb"])

                print(f"Remote database newer (released {uniprot.remote_reldate})")

                
                new_spdb = uniprot.retr_spdb(dest_dir)
                new_spdb = uniprot.decompress(new_spdb)

                new_taxdb = f"{new_spdb}.taxdb"
                
                #%% Build taxonomy database
                try_build = subprocess.call([
                    sys.executable, os.path.join(self.SCRIPTS,'build_taxdb.py'), 
                    '-db', new_spdb, 
                    '-t', yamlcfg['taxonomy_all_tab'], 
                    '-o', new_taxdb
                    ])
                if try_build: ## STOP if buildtaxdb exit code != 0
                    print (f"> {output_cols['RED']}[ERROR]{output_cols['RESET']} Fatal error encountered while building taxonomy database.")
                    sys.exit(1)

                self.update_config(yamlcfg, new_spdb, new_taxdb, uniprot.remote_reldate)
                

class SPDBExpander(Module, Head, LoggingEntity, ErrorHandler):
    def __init__(self):
        super().__init__(self.__class__)
        self.HOME, self.SCRIPTS = pkgloc()

        self.argparser = self.generate_argparser()
        self.parsed_args = self.argparser.parse_args()
        self.config.load_cmdline( self.parsed_args )

        self.taxdb = SPDBTaxonomy(self.config.get("taxdb"))
    
    def generate_argparser(self):
        parser = argparse.ArgumentParser()
        parser.add_argument("mod", nargs="*")
        parser.add_argument('-db','--spdb', metavar="swissprot_fasta",action=PathAction, required=False, help = "A file in FASTA format containing the swissprot database to be used. FASTA headers must be in standard swissprot format [... OS= ... GN= | PE= | OX= | EOL ]. See README for more information and examples.")
        parser.add_argument('-t','--taxdb', metavar='taxonomy_db', action=PathAction, required=False, help='The path to the taxonomy database associated with -db|--spdb')
        parser.add_argument('-p','--proteins', metavar='proteins_to_add', action=PathAction, required=True, help='A file in FASTA format that contains the protein sequences that you would like to add to the swissprot-style database.')
        parser.add_argument('-l','--lineages', metavar='lineages_to_add', action=PathAction, required=True, help='A two-column, tab-separated list of all unique OSs present in -p|--proteins (column 1) and semicolon-separated lineage information (column2).')
        parser.add_argument('-o','--output', metavar='output', action=OutputPathStore, required=True, help='The output path for your expanded swissprot-style database.')
        parser.add_argument('-d','--defaults', action='store_true', required=False, help='Include this flag if you would like scgid to update package settings to reflect the newly created databases as your defaults.')
        
        return parser

    def run(self):

        old_spdb = self.config.get('spdb')
        old_taxdb  = self.config.get('taxdb')

        new_spdb = self.config.get('output')
        new_taxdb = f"{self.config.get('output')}.taxdb"

        # Write expanded SPDB fasta
        self.logger.info(f"Loading and expanding SPDB [{old_spdb}]")
        with open(new_spdb, 'w') as f_out:
            for f_in in [old_spdb, self.config.get("proteins")]:
                for line in open(f_in, 'r'):
                    f_out.write(line)

        # Write expanded SPDB taxonomy database
        self.logger.info(f"Loading and expanding SP taxdb [{old_taxdb}]")
        with open(new_taxdb, 'w') as f_out:
            self.taxdb.expand( self.config.get("lineages") )
            self.taxdb.write(f_out)

        # Update config.yaml if -d|--defaults is True
        if self.config.defaults:
            config = FileConfig()
            config.load_yaml()
            config.settings["default_spdb"] = new_spdb
            config.settings["default_taxdb"] = new_taxdb
            
            # Overwrite existing config.yaml (setting path == None overwrites same path read from)
            config.write_yaml(path = None)

            self.simplelogger.info(f"\nOverwrote defaults in `{config.path}`:\n{'-'*40}")
            self.simplelogger.info(f"default_spdb\t{new_spdb}")
            self.simplelogger.info(f"default_taxdb\t{new_taxdb}\n")