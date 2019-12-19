import os
import sys
import numpy as np
import time
import io
import re
import yaml
import gzip
import shutil
from ftplib import FTP
from scgid.modcomm import pkgloc
from scgid.library import CURSOR_UP_ONE, ERASE_LINE

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
