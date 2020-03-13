import os
import sys
import numpy as np
import time
import io
import re
import yaml
import gzip
import shutil
import logging
import ast
from ftplib import FTP
from scgid.modcomm import pkgloc, LoggingEntity, ErrorHandler, logger_name_gen
from scgid.library import CURSOR_UP_ONE, ERASE_LINE
from scgid.error import ModuleError, ErrorClassNotImplemented, Ok, check_result

class MalformedTaxonomyDatabaseError(ModuleError):
    def __init__(self, taxdbpath):
        super().__init__()
        self.msg = f"Taxonomy database at `{taxdbpath}` improperly formatted. Run `build_taxdb.py` to generate a new one."
        self.catch()

class MalformedLineageFileError(ModuleError):
    def __init__(self):
        super().__init__()
        self.msg = "Malformed lineage file supplied as argument to `-l|--lineages`."
        self.catch()

def ftp_retr_progress (block, dest, tsize):
    with open(dest, 'ab') as d:
        ret = d.write(block)
    csize = os.path.getsize(dest)
    prog = int(np.floor(((csize)/(tsize))*50))
    sys.stdout.write(CURSOR_UP_ONE)
    sys.stdout.write(ERASE_LINE)
    #sys.stdout.write("["+"#"*(prog)+" "*(50-prog)+"] "+str(prog*2)+"%\n")
    sys.stdout.write(f"[{'#'*prog}{' '*(50-prog)}] {prog*2}%\n")
    return(ret)


def ftp_retr_and_report (ftp_inst, src, dest):
    tsize = ftp_inst.size(src)
    print (f"> Retrieving {src} from {ftp_inst.host} ...\n[ {' '*50} ] 0%")
    ftp_inst.retrbinary("RETR %s" % (src), lambda block: ftp_retr_progress(block, dest, tsize))

class UniprotFTP(object):
    def __init__(self):
        self.server = "ftp.uniprot.org"
        self.ftp_object = FTP(self.server)
        self.cd_to = "pub/databases/uniprot/current_release/knowledgebase/complete/"
        self.spdb_fasta = "uniprot_sprot.fasta.gz"
        self.path_to_reldate_txt = "reldate.txt"
        self.pull_reldate_str = lambda string: re.search("[0-9]{2}-[a-zA-Z]{3}-[0-9]{4}", string)

    def __enter__(self):
        self.ftp_object.login()
        self.ftp_object.cwd(self.cd_to)
        self.remote_reldate = self.get_remote_reldate()

        return self

    def __exit__(self, *kwargs):
        self.ftp_object.__exit__()
    
    def get_remote_reldate(self):
        with io.StringIO() as buffer:
            self.ftp_object.retrlines(f"RETR {self.path_to_reldate_txt}", buffer.write)
            text = buffer.getvalue()
        
        s = self.pull_reldate_str(text)
        
        return s.group(0)

    def needs_retr (self, config_path, remote_reldate):
        with open(config_path, 'r') as cfg:
            config_dict = yaml.load(cfg, Loader=yaml.BaseLoader)
        
        s = self.pull_reldate_str(config_dict["default_spdb"])
        if s is None:
            return True

        else:
            if not os.path.isfile(config_dict["default_spdb"]):
                raise FileNotFoundError(f"SPDB listed in comfig.yaml does not exist")

            current_reldate = time.mktime( time.strptime( s.group(0), "%d-%b-%Y" ) )
            
            if current_reldate < time.mktime( time.strptime(self.remote_reldate, "%d-%b-%Y") ):
                return True
            
            else:
                return False

    def decompress (self, path):
        d_path = path.replace(".gz", "")
        with gzip.open(path, 'rb') as f_in:
            with open(path.replace(".gz", ""), 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        return d_path       

    def retr_spdb (self, dest_dir):
        dest = os.path.join(
            dest_dir,
            self.spdb_fasta.replace(".fasta.gz", f"_{self.remote_reldate}.fasta.gz")
        )
        if os.path.isfile(dest):
            os.remove(dest)
        ftp_retr_and_report(self.ftp_object, self.spdb_fasta, dest)
        return dest

class SPDBTaxonomy(LoggingEntity, ErrorHandler):
    def __init__ (self, path_to_taxdb):
        self.logger = logging.getLogger ( logger_name_gen() )
        self.taxdb = {}
        self.load(path_to_taxdb)

    @check_result
    def load(self, path_to_taxdb):
        with open(path_to_taxdb,'r') as f:
            try:
                self.taxdb = ast.literal_eval( f.read() )
                return Ok(None)

            except SyntaxError:
                raise MalformedTaxonomyDatabaseError(taxdbpath = self.head.config.get("taxdb"))

            except:
                raise ErrorClassNotImplemented

    def __repr__(self):
        return '\n'.join([f"{k}: {v}" for k,v in self.taxdb.items()])

    def expand(self, path_to_lineage) -> None:
        with open(path_to_lineage) as f_in:
            lines = f_in.readlines()
        
        # Split file into key,lineage for taxdb (hopefully)
        spl = [ [i.strip() for i in l.split("\t")] for l in lines]

        # Check length of splits to ensure that there are TWO columns, one key with one tab-separated lineage
        right_ncols = [len(s) == 2 for s in spl]
        
        if not all(right_ncols):
            raise MalformedLineageFileError
        
        else:
            lines_parsed = {i[0]: i[1] for i in spl} 
            self.taxdb.update(lines_parsed)

        return None

    # Takes open file object - use within with statment
    # This should be done with most functions that read/write to files I think - makes tests with StringIO easy
    def write(self, f_out) -> None:
        f_out.write("{\n")
        f_out.write(",\n".join([f"'{key}': '{lin}'" for key,lin in self.taxdb.items()]))
        f_out.write("\n}")
        return None
