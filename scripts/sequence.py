import sys
import numpy
import re
from collections import OrderedDict

def reverse (string):
    chars = ["x"]*len(string)
    pos = len(string) - 1
    for l in string:
        chars[pos] = l
        pos -= 1
    return ''.join(chars)

def complement(string):
    convert = {
            'A': 'T',
            'T': 'A',
            'A': 'U',
            'U': 'A',
            'G': 'C',
            'C': 'G',
            'N': 'N'
        }
    comp = [convert[l] for l in string]
    return ''.join(comp)

def revcomp(string):
    convert = {
            'A': 'T',
            'T': 'A',
            'G': 'C',
            'C': 'G',
            'N': 'N'
        }
    comp = [convert[l] for l in string]
    revcomp = comp[::-1]

    return ''.join(revcomp)
    

class DNASequenceCollection(object):

    def __init__(self):
        self.index = {}
        self.seqtype = DNASequence

    def get(self, header):
        return self.index[header]
    
    def pop(self, header):
        return self.index.pop(header)

    def seqs(self):
        return self.index.values()

    def from_dict(self, d):
        self.index.clear()
        self.index.update(d)
        return self
    
    def remove_small_sequences (self, minlen):
        self.index = {header: seqobj for header, seqobj in self.index.items() if seqobj.length >= minlen}

    def gc_cov_filter(self, gc_range, coverage_range):
        gc_min, gc_max = gc_range
        coverage_min, coverage_max = coverage_range
        
        sort = {
            "keep": {},
            "dump": {}
        }

        for key, seqobj in self.index.items():
            if seqobj.gc >= gc_min and seqobj.gc <= gc_max and seqobj.coverage >= coverage_min and seqobj.coverage <= coverage_max:
                sort["keep"][key] = self.get(key)
            else:
                sort["dump"][key] = self.get(key)
        
        return sort

    def header_list_filter(self, header_list):
        seq_dict = dict()
        for header in header_list:
            seq_dict[header] = self.index.get(header)
        return DNASequenceCollection().from_dict(
            { h: seq_dict[h] for h in sorted(seq_dict) }
        )

    def from_fasta(self, fasta, spades = False):
        header_pattern = re.compile("^>(.+)")
        header = str()
        sequence = str()

        with open(fasta, 'r') as f:
            for line in f.readlines():
                line = line.strip()
                if len(line) == 0:
                    continue

                m = re.match(header_pattern, line)
                if m is not None:
                    if len(header) > 0:
                        
                        if header in self.index:
                            raise KeyError("FASTA has duplicated headers")

                        self.index[header] = self.seqtype(header, sequence, spades)
                        header = m.group(1)
                        sequence = str()
                    else:
                        header = m.group(1)
                else:
                    sequence += line
            # Add the last sequence to the OrderedDict
            self.index[header] = DNASequence(header, sequence, spades)

        return self

    def rekey_by_shortname (self):
        self.index = {"_".join(k.split("_")[0:2]): v for k,v in self.index.items()}

    def write_fasta (self, outpath):
        with open(outpath, 'w') as o:
            for _,s in self.index.items():
                o.write(f"{s.to_fasta()}\n")

class DNASequence(object):

    def __init__(self, header, string, spades = False):
        self.header = header
        self.string = string
        self.length = len(self.string)

        self.coverage = None
        self.shortname = None
        if spades:
            spl = header.split('_')
            self.coverage = float(spl[5])
            self.shortname = "_".join(spl[0:2])

        self.gc = float()
        self.gcCount = 0
        for letter in self.string:
            if letter == "G" or letter == "C":
                self.gcCount += 1
        self.gc = float(self.gcCount)/float(self.length)

    def to_fasta(self):
        return f">{self.header}\n{self.string}"
    
    def revcomp(self, inplace = False):
        if inplace:
            self.string = revcomp(self.string)
        else:
            return revcomp(self.string)
    
    def transcribe(self):
        trans = {
                'A':'U',
                'T':'A',
                'G':'C',
                'C':'G',
                'N':'N'
                }
        transcript = ''.join([trans[l] for l in self.string])
        return transcript

class AASequenceCollection(object):

    def __init__(self):
        self.index = {}

    def get(self, header):
        return self.index[header]

    def seqs(self):
        return self.index.values()

    def from_dict (self, odict):
        self.index.update(odict)
        return self

    def from_fasta(self, fasta):
        header_pattern = re.compile("^>(.+)")
        header = str()
        sequence = str()

        with open(fasta, 'r') as f:
            for line in f.readlines():
                line = line.strip()
                if len(line) == 0:
                    continue

                m = re.match(header_pattern, line)
                if m is not None:
                    if len(header) > 0:
                        self.index[header] = AASequence(header, sequence)
                        header = m.group(1)
                        sequence = str()
                    else:
                        header = m.group(1)
                else:
                    sequence += line
            ## add the last sequence to the list
            self.index[header] = AASequence(header, sequence)

        return self

class AASequence(object):

    def __init__(self, header, string):
        self.header = header
        self.string = string
        self.length = len(self.string)
    
    def to_fasta(self):
        return f">{self.header}\n{self.string}"

# Old shit
class Sequence(object):
    def __init__(self, label, seq, seq_type = "nucl", contig_info = False):
        self.label = label
        self.id = label

        if seq_type == "prot" and contig_info == True:
            try:
                self.coverage = float(".".join(label.split("_")[5].split(".")[0:len(label.split("_")[5].split("."))-1]))
            except:
                self.coverage = float(".".join(label.split("_")[5].split(".")[0:1]))
        elif seq_type == "nucl" and contig_info == True:
            self.coverage = float(label.split("_")[5])
        elif seq_type == "nucl" or seq_type == "prot" and contig_info == False:
            self.coverage = None
        if contig_info == True:
            self.name = "_".join(self.label.split("_")[0:2])
            self.shortname = self.name
        self.sequence = seq
        self.length = len(self.sequence.strip())

        ## Check format to shorten Node names if possible
        m = re.match("NODE[_][0-9]+[_]length[_][0-9]+[_]cov[_][0-9]+",self.label)
        if m is not None:
            spl = self.label.split("_")
            sid = '_'.join(spl[0:2])
            sid = sid.replace("NODE_","N")
            self.id = sid

        self.nr_domain = ""
        self.nr_phylum = ""
        if seq_type is "nucl":
            self.gc = float()
            self.gcCount = 0
            for letter in self.sequence:
                if letter == "G" or letter == "C":
                    self.gcCount += 1
                self.gc = float(self.gcCount)/float(self.length)
    def outFasta(self):
        out= ">" + self.label + "\n"+ self.sequence
        return out
    def revcomp(self):
        comp=""
        for i in self.sequence:
            if i == "A":
                comp+="T"
            elif i == "T":
                comp+="A"
            elif i == "C":
                comp+="G"
            elif i == "G":
                comp+="C"
        out = reverse(comp)
        return out

def revcomp_str(dna_string):
    comp=""
    for i in dna_string:
        if i == "A":
            comp+="T"
        elif i == "T":
            comp+="A"
        elif i == "C":
            comp+="G"
        elif i == "G":
            comp+="C"
    out = reverse(comp)
    return out
#read sequences in a fasta file into individual sequence() objects
def readFastq(infile):
    allSeqs=[]
    lengths=[]
    p=1
    m=0
    for i in open(infile).readlines():
        if m == 100:
            break
        if p == 2:
            lengths.append(len(i.strip()))
        if p == 4:
            p = 0
        p+=1
        m+=1
        #print i.split("\n")[1]
        #print "XXXXXXX"
        #p+=1
    return numpy.mean(lengths)

def getFastaHeaders (fasta):
    out = list()
    with open(fasta,'r') as f:
        for line in f:
            if line.startswith(">"):
               line = line.replace(">","",1)
               out.append(line.strip())
    return out

def readFasta (fasta, seq_type = "nucl", contig_info = False):
    all_seqs = []
    header_pattern = re.compile("^>(.+)")
    label = str()
    seq = str()

    with open(fasta, 'r') as f:
        for line in f.readlines():
            line = line.strip()
            if len(line) == 0:
                continue

            m = re.match(header_pattern, line)
            if m is not None:
                if len(label) > 0:
                    all_seqs.append( Sequence(label, seq, seq_type, contig_info) )
                    label = m.group(1)
                    seq = str()
                else:
                    label = m.group(1)
            else:
                seq += line
        ## add the last sequence to the list
        all_seqs.append( Sequence(label, seq, seq_type, contig_info) )
    return all_seqs

def readFasta_dumber(file, seq_type = "nucl", contig_info = False):
    allSeqs=[]
    for i in open(file).read().split(">"):
        if len(i) == 0:
            continue
        l=0
        seq=""
        #print i
        for x in i.split("\n"):
            if l==0:
                lab=x
            elif l!=0:
                seq+=x
            l+=1
        if seq_type == "prot" and contig_info == True:
            allSeqs.append(Sequence(lab,seq,"prot",True))
        elif seq_type == "nucl" and contig_info == True:
            allSeqs.append(Sequence(lab,seq,"nucl",True))
        elif seq_type == "prot" and contig_info == False:
            allSeqs.append(Sequence(lab,seq,"prot",False))
        elif seq_type == "nucl" and contig_info == False:
            allSeqs.append(Sequence(lab,seq,"nucl",False))
            #allSeqs.append(contig(lab,seq))
    return allSeqs

def get_attribute_list (instanceList, attribute):
    listToOutput = []
    for i in instanceList:
        listToOutput.append(getattr(i,attribute))
    return listToOutput

