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

class DNASequenceCollection(object):
    
    def __init__(self):
        self.odict = OrderedDict()
    
    def get(self, header):
        return self.odict[header]
    
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
                        self.odict[header] = DNASequence(header, sequence, spades) 
                        header = m.group(1)
                        sequence = str()
                    else:
                        header = m.group(1)
                else:
                    sequence += line
            ## add the last sequence to the list
            self.odict[header] = DNASequence(header, sequence, spades)
        
        return self

    def rekey_by_shortname (self):
        self.odict = {"_".join(k.split("_")[0:2]): v for k,v in self.odict.items()}

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

class AASequenceCollection(object):
    
    def __init__(self):
        self.odict = OrderedDict()
    
    def get(self, header):
        return self.odict[header]
    
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
                        self.odict[header] = AASequence(header, sequence) 
                        header = m.group(1)
                        sequence = str()
                    else:
                        header = m.group(1)
                else:
                    sequence += line
            ## add the last sequence to the list
            self.odict[header] = AASequence(header, sequence)
        
        return self

class AASequence(object):

    def __init__(self, header, string):
        self.header = header
        self.string = string
        self.length = len(self.string)

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

