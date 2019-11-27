import subprocess
import re
import sys
import os
import pickle
from io import StringIO
from collections import OrderedDict
from scripts.sequence import DNASequence, DNASequenceCollection, AASequence, AASequenceCollection

def is_fasta(f, strict=False, verbose = False):
    lines = [l.strip() for l in open(f).readlines()]

    fail = False
    has_headers = True
    made_up_of_seqs = True
    just_had_header = False
    idx = 0

    #lines = go_through_list(lines)
    for line in lines:
        if not just_had_header:
            if re.match('^>.+',line) is not None:
                has_headers = True
                just_had_header = True
            else:
                has_headers = False
        elif just_had_header:
            if strict:
                if re.match('^[-atcguATCGUgalmfwkqespvicyhrndtbzxuGALMFWKQESPVICYHRNDTBZXU*]+$',line) is not None:
                    made_up_of_seqs = True
                else:
                    made_up_of_seqs = False
            else:
                if re.match('^[-A-Za-z*]+$',line) is not None:
                    made_up_of_seqs = True
                else:
                    made_up_of_seqs = False
            try:
                if re.match('^>.+',lines[idx+1]) is not None:
                    just_had_header = False
            except:
                break
        if not has_headers or not made_up_of_seqs:
            fail = True
            if verbose:
                print( f"Deviation from fasta format on line #{idx+1}" )
        idx+=1
    if not fail:
        return True
    else:
        return False

def gff3_to_fasta(gff3, outname):
    prots = []
    recording = False
    for line in open(gff3).readlines():
        if recording is False:
            s = re.search("[#][ ]protein[ ]sequence[ ][=][ ]",line)
            if s is not None:
                recording = True
                s = re.search("[#][ ]protein[ ]sequence[ ][=][ ][[]([A-Za-z]+)",line)
                protein_sequence = s.group(1)
        elif recording is True:
            s = re.search("[#][ ]end[ ]gene[ ]([A-Za-z0-9_.]+)",line)
            if s is not None:
                recording = False
                prots.append(AASequence(s.group(1), protein_sequence))
                protein_sequence = ""
            else:
                s = re.search("[#][ ]([a-zA-Z]+)",line)
                protein_sequence += s.group(1)
    with open(outname,'w') as fasta:
        for i in prots:
            fasta.write(i.to_fasta()+"\n")

def log_command_line_errors (err, log_inst):
    if err != '':
        print(err)
        return 0
        err = err.split('\n')
        for i in err:
            if len(i) > 0:
                log_inst.critical(i)

def log_command_line_stdout (out, log_inst):
    for l in [x.strip() for x in StringIO(out.decode("utf-8")).readlines()] + [""]: ### Add extra line to deal with weird ESOM OOM stdout that cuts off
        log_inst.info(l)
    return 0

def subprocessC (args):
    p = subprocess.call(args)
    if p != 0:
        sys.exit(1)

def subprocessP (args, log_inst, log_stdout=False):
    p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out,err = p.communicate()
    if log_stdout:
        assert isinstance(log_inst,list), "Internal error. In order to log stdout, a list of logs must be passed to subprocessP"
        log_command_line_stdout(out,log_inst[1])
        if p.returncode != 0:
            log_command_line_errors(err, log_inst[0])
            sys.exit(1)
    else:
        if p.returncode != 0:
            log_command_line_errors(err, log_inst)
            sys.exit(1)
    return out

def subprocessT (args):
    p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()

def pickle_loader (pklFile):
    try:
        while True:
            yield pickle.load(pklFile)
    except EOFError:
            pass
#%%
def pkl_fasta_in_out (org_fname, seq_type = "nucl", spades = False):
    #pkl_fname = org_fname.split("/")[-1]+'.pkl'
    pkl_fname = "{}.pkl".format(org_fname)
    #if os.path.isfile(os.getcwd()+"/"+pkl_fname):
    if os.path.isfile(pkl_fname) and os.stat(org_fname).st_ctime < os.stat(pkl_fname).st_ctime:
        #print "LOADING FROM
        objlist = OrderedDict()
        with open(pkl_fname,'rb') as infile:
            headers = []
            obj = []
            for header, sequence in pickle_loader(infile):
                headers.append(header)
                obj.append(sequence)
            objlist = dict(zip(headers,obj))
            if seq_type == "nucl":
                return DNASequenceCollection().from_dict(objlist)
            else:
                return AASequenceCollection().from_dict(objlist)

    else:
        #print "MAKING NEW PKL"
        if seq_type == "nucl":
            objlist = DNASequenceCollection().from_fasta(org_fname, spades)
        else:
            objlist = AASequenceCollection().from_fasta(org_fname)

        with open(pkl_fname,'wb') as output:
            for seq in objlist.odict.items():
                pickle.dump(seq, output, pickle.HIGHEST_PROTOCOL)
    return objlist


def file_grep (pattern, file, mode='first'):
    if mode not in ['first','multiple']:
        raise NameError("Invalid option, "+mode)
    result_list = []
    with open(file,'r') as f:
        for line in f.readlines():
            if re.search(pattern,line) is not None:
                if mode == 'first':
                    return line.strip()
                elif mode == 'multiple':
                    result_list.append(line.strip())
    return result_list
