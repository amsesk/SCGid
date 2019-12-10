import os
from scgid.modcomm import pkgloc
from scgid.db import UniprotFTP

HOME, _ = pkgloc()
uniprot = UniprotFTP()

uniprot.login()

if uniprot.needs_retr (
    os.path.join(HOME, "config.yaml"),
    uniprot.get_remote_reldate()
    ):

    print ("Needs to be retrieved")
else:
    print ("Up to date")

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
