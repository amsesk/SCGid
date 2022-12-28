import sys
import numpy
import re
import numpy as np
from collections import OrderedDict
from scgid.error import ErrorClassNotImplemented

SPADES_HEADER_PATTERN = re.compile("NODE_[0-9]+_length_[0-9]+_cov_[0-9.]+")
IUPAC_DEGENERATES = ['R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V']

def reverse(string):
    chars = ["x"] * len(string)
    pos = len(string) - 1
    for l in string:
        chars[pos] = l
        pos -= 1
    return ''.join(chars)


def complement_dna(string):
    convert = {
        'A': 'T',
        'T': 'A',
        'G': 'C',
        'C': 'G',
        'N': 'N'
    }
    comp = [convert[l] if l in convert.keys() else 'N' for l in string]
    return ''.join(comp)


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
    comp = []
    for l in string:
        l = l.upper()
        try:
            comp.append(convert[l])
        except KeyError:
            if l in IUPAC_DEGENERATES:
                print(f"WARNING: Degenerate base designation '{l}' converted to 'N'")
                comp.append("N")
            else:
                raise IOError(f"Sequence string contains illegal character: '{l}'")
    revcomp = comp[::-1]

    return ''.join(revcomp)

'''Deprecated I think
class DnaSequence(object):

    def __init__(self, accession, description, string):
        self.accession = accession
        self.description = description
        self.string = "".join([c.upper() for c in string])
        self.length = len(string)

    def gc_content(self):
        gc_count = 0
        for letter in self.string:
            if letter in ["G", "C"]:
                gc_count += 1
        gc_content = float(gc_count) / float(self.length)
        return gc_content

    def to_fasta(self):
        header = f">{self.accession} {self.description}".strip()
        return f"{header}\n{self.string}"

    def revcomp(self, inplace=False):
        if inplace:
            self.string = revcomp(self.string)
        else:
            return revcomp(self.string)

    def transcribe(self):
        trans = {
            'A': 'U',
            'T': 'A',
            'G': 'C',
            'C': 'G',
            'N': 'N'
        }
        transcript = []
        for l in self.string:
            try:
                transcript.append(trans[l])
            except KeyError:
                if l in IUPAC_DEGENERATES:
                    print(f"WARNING: Degenerate base designation '{l}' converted to 'N'")
                    transcript.append("N")
                else:
                    raise IOError(f"Sequence string contains illegal character: '{l}'")


        return ''.join(transcript)
'''

class Sequences(object):

    def __init__(self, seqtype):
        self.index = OrderedDict()
        self.seqtype = seqtype

    def from_fasta(fasta_handle, seqtype):
        header_pattern = re.compile("^>(.+)")

        sequences = Sequences(seqtype)

        accession = ""
        description = ""
        sequence = ""

        for line in fasta_handle:
            line = line.strip()
            if len(line) == 0:
                continue

            m = re.match(header_pattern, line)

            ##############################
            ## Are we on a header line? ##
            ##############################
            if m is not None:
                ################
                ## Yes we are ##
                ################

                #############################################
                ## Are we currently working on a sequence? ##
                ## If so, add the current sequence to the  ##
                ## index before moving on to the new one   ##
                #############################################
                if len(accession) != 0:
                    sequences.index[accession] = sequences.seqtype(
                        accession, description, sequence)
                    accession = ""
                    description = ""
                    sequence = ""

                ###############################
                ## Now move on to the header ##
                ## on this line              ##
                ###############################
                header = m.group(1).split(" ")
                if len(header) > 1:
                    accession = header[0]
                    description = " ".join(header[1:])
                else:
                    accession = header[0]
                    description = ""

                # Make sure this accession is not yet in the index
                if accession in sequences.index:
                    raise KeyError("FASTA has duplicated headers")

            else:
                ################################
                ## This is not a header line, ##
                ## so add to the sequence    ###
                ################################
                sequence += line

        ########################################
        ## Add the last sequence to the index ##
        ########################################
        if accession in sequences.index:
            raise KeyError("FASTA has duplicated headers")
        sequences.index[accession] = sequences.seqtype(
            accession, description, sequence)
        print(sequences.index[accession].string)

        return sequences

    def has_spades_headers(self):
        spades_header_pattern = SPADES_HEADER_PATTERN
        for k in self.index:
            if not re.match(spades_header_pattern, k):
                return False

        return True

    def seqs(self):
        return self.index.values()

    def __str__(self):
        return "\n".join([s.to_fasta() for _k, s in self.index.items()])


class DNASequenceCollection(object):

    def __init__(self):
        self.index = {}
        self.seqtype = DNASequence

    def get(self, header):
        return self.index[header]

    def pop(self, header):
        return self.index.pop(header)

    def check_size(self, cutoff):
        if len(self.seqs()) < cutoff:
            return False
        else:
            return True

    def seqs(self):
        return self.index.values()

    def headers(self):
        return self.index.keys()

    def from_dict(self, d):
        self.index.clear()
        self.index.update(d)
        return self

    def remove_small_sequences(self, minlen):
        self.index = {header: seqobj for header,
                      seqobj in self.index.items() if seqobj.length >= minlen}

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
            {h: seq_dict[h] for h in sorted(seq_dict)}
        )

    def from_fasta(self, fasta, spades=False, coverage_dict=None):
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

                        if coverage_dict is not None:
                            try:
                                coverage = coverage_dict.index[header]
                            except KeyError:
                                raise scgid.gct.MalformedCoverageTableError(f"Contig `{header}` is not present in the supplied coverage table.")
                        else:
                            coverage = np.nan

                        if header in self.index:
                            raise KeyError("FASTA has duplicated headers")

                        self.index[header] = self.seqtype(
                            header, sequence, spades, coverage)
                        header = m.group(1)
                        sequence = str()
                    else:
                        header = m.group(1)
                else:
                    sequence += line

            # Add the last sequence to the OrderedDict
            if coverage_dict is not None:
                try:
                    coverage = coverage_dict.index[header]
                except KeyError:
                    raise scgid.gct.MalformedCoverageTableError(f"Contig `{header}` is not present in the supplied coverage table.")
            else:
                coverage = np.nan

            if header in self.index:
                raise KeyError("FASTA has duplicated headers")
            self.index[header] = DNASequence(
                header, sequence, spades, coverage)

        return self

    def rekey_by_shortname(self):
        self.index = {
            "_".join(k.split("_")[0:2]): v for k, v in self.index.items()}

    def write_fasta(self, outpath):
        with open(outpath, 'w') as o:
            for _, s in self.index.items():
                o.write(f"{s.to_fasta()}\n")


class DNASequence(object):

    def __init__(self, header, string, spades=False, coverage=np.nan):
        self.header = header
        self.string = string
        self.length = len(self.string)

        self.coverage = None
        self.shortname = None

        try:
            spl = header.split('_')
            self.shortname = "_".join(spl[0:2])
        except IndexError:
            self.shortname = self.header
        except:
            raise ErrorClassNotImplemented

        if np.isnan(coverage):
            try:
                spl = header.split('_')
                self.coverage = float(spl[5])
            except IndexError:
                pass
            except ValueError:
                pass
            except:
                raise ErrorClassNotImplemented
        else:
            self.coverage = coverage

        self.gc = float()
        self.gcCount = 0
        for letter in self.string:
            if letter.upper() == "G" or letter.upper() == "C":
                self.gcCount += 1
        self.gc = np.divide(self.gcCount, self.length)

    def to_fasta(self):
        return f">{self.header}\n{self.string}"

    def revcomp(self, inplace=False):
        if inplace:
            self.string = revcomp(self.string)
        else:
            return revcomp(self.string)

    def transcribe(self):
        trans = {
            'A': 'U',
            'T': 'A',
            'G': 'C',
            'C': 'G',
            'N': 'N'
        }
        transcript = []
        for l in self.string:
            l = l.upper()
            try:
                transcript.append(trans[l])
            except KeyError:
                if l in IUPAC_DEGENERATES:
                    print(f"WARNING: Degenerate base designation '{l}' converted to 'N'")
                    transcript.append("N")
                else:
                    raise IOError(f"Sequence string contains illegal character: '{l}'")


        return ''.join(transcript)

class AASequenceCollection(object):

    def __init__(self):
        self.index = {}

    def get(self, header):
        return self.index[header]

    def seqs(self):
        return self.index.values()

    def from_dict(self, odict):
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
            # add the last sequence to the list
            self.index[header] = AASequence(header, sequence)

        return self


class AASequence(object):

    def __init__(self, header, string):
        self.header = header
        self.string = string
        self.length = len(self.string)

    def to_fasta(self):
        return f">{self.header}\n{self.string}"

'''
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

        # Check format to shorten Node names if possible
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
# read sequences in a fasta file into individual sequence() objects
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
        # print i.split("\n")[1]
        # print "XXXXXXX"
        # p+=1
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
        # add the last sequence to the list
        all_seqs.append( Sequence(label, seq, seq_type, contig_info) )
    return all_seqs

def readFasta_dumber(file, seq_type = "nucl", contig_info = False):
    allSeqs=[]
    for i in open(file).read().split(">"):
        if len(i) == 0:
            continue
        l=0
        seq=""
        # print i
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
            # allSeqs.append(contig(lab,seq))
    return allSeqs

def get_attribute_list (instanceList, attribute):
    listToOutput = []
    for i in instanceList:
        listToOutput.append(getattr(i,attribute))
    return listToOutput

'''
