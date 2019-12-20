import os
import time
import yaml
import sys
import subprocess
import argparse
from scgid.modcomm import pkgloc, Head, ErrorHandler, LoggingEntity, logger_name_gen
from scgid.parsers import PathAction, OutputPathStore
from scgid.module import Module
from scgid.db import UniprotFTP
from scgid.library import output_cols
from scgid.error import ModuleError

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
        self.start_logging()
        self.check_path_args()

        # Write expanded SPDB
        with open(self.config.get("output"), 'w') as f_out:
            for f_in in [self.config.get("spdb"), self.config.get("proteins")]:
                for line in open(f_in, 'r'):
                    f_out.write(line)
        
        with open(self.config.get("lineages", 'r')) as supplied_lineages:
            lines = supplied_lineages.readlines()
            
            # Split file into key and lineage for taxdb (hopefully)
            spl = [l.split("\t") for l in lines]

            # Check length of splits to ensure that there are TWO columns, one key with one tab-separated lineage
            right_ncols = [len(s) != 2 for s in spl]
            if not any(right_ncols):
                return ModuleError(f"Supplied lineage file is not correctly formatted. Two-column tab-separated list expected, but found rows with ncols {','.join([i for i in right_ncols if i != 2])}")
            
            else:
                


'''
bin_dir =os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))

if len(settings.path_to_spdb) == 0:
    db_loc = os.path.abspath(sys.argv[1])
else:
    db_loc = os.path.split(settings.path_to_spdb)[0]

try:
    os.chdir(db_loc)
except:
    os.mkdir(db_loc)
    os.chdir(db_loc)

# Open FTP connection
server = "ftp.uniprot.org"
file_to_get = "uniprot_sprot.fasta.gz"

ftp = FTP(server)
ftp.login()
ftp.cwd("pub/databases/uniprot/current_release/knowledgebase/complete/")

# Get swissprot version information from reltdate.txt
ftp.retrbinary("RETR %s" % "reldate.txt", open("reldate.txt","wb").write)
with open("reldate.txt") as f:
    dateline = f.readlines()[1]
    s = re.search("[0-9]{2}-[a-zA-Z]{3}-[0-9]{4}", dateline)
    date = s.group(0)
    fname = "uniprot_sprot_"+date+".fasta.gz"
ts = datetime.datetime.strptime(date, "%d-%b-%Y").strftime("%s")
os.remove("reldate.txt")

decomp_cmd = ["gzip", "-fd", fname]

# Download swissprot if necessary ... necessary if seconds from epoch of current local version < FTP version

## account for <scgid init> case where there is no current spdb, set version to 0 seconds after epoch
try:
    current_spdb = settings.spdb_version
except:
    current_spdb = "01-Jan-1970"

changes_needed = False

if os.path.isfile("uniprot_sprot_"+current_spdb+".fasta"):
    current_ts = datetime.datetime.strptime(current_spdb, "%d-%b-%Y").strftime("%s")
    if ts == current_ts:
        print "> Nothing to do... swissprot database up to date."
    else:
        changes_needed = True
        os.remove("uniprot_sprot_"+current_spdb+".fasta")
        ftp_retr_and_report(ftp, server, file_to_get)
        print "> Decompressing %s" % (file_to_get)
        subprocess.call(decomp_cmd)
        db_fasta = '.'.join(fname.split('.')[0:-1])
        db_path = os.path.join(os.getcwd(), db_fasta)
        print "> Updated to new version of swissprot database at "+db_path
        
else:
    changes_needed = True
    ftp_retr_and_report(ftp, server, file_to_get, fname)
    print "> Decompressing %s" % (file_to_get)
    subprocess.call(decomp_cmd)
    db_fasta = '.'.join(fname.split('.')[0:-1])
    db_path = os.path.join(os.getcwd(), db_fasta)
    print "> Downloaded current version of swissprot database to "+db_path
ftp.quit()

## update or write values to settings.py
if changes_needed:
    try:
        replace_line_by_pattern(os.path.join(bin_dir,"settings.py"), "spdb_version=", "spdb_version=\"%s\"" % (date))
    except:
        with open(os.path.join(bin_dir,"settings.py"),'a') as s:
            s.write("spdb_version=\""+date+"\"")
            s.write("\n")

    try:
        replace_line_by_pattern(os.path.join(bin_dir,"settings.py"), "path_to_spdb=", "path_to_spdb=\"%s\"" % (db_path))
    except:
        with open(os.path.join(bin_dir, "settings.py"),'a') as s:
            s.write("path_to_spdb=\""+db_path+"\"")
            s.write("\n")

'''
