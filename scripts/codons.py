import argparse
import sys
import re
import pandas as pd
from collections import namedtuple
from scripts.module import Module
from scripts.modcomm import LoggingEntity, Head
from scripts.reuse import ReusableOutput, ReusableOutputManager, augustus_predict, nucleotide_blast, protein_blast
from scripts.dependencies import CaseDependency
from scripts.parsers import PathAction
from scripts.sequence import DNASequenceCollection, DNASequence, revcomp, complement

SYNONYMOUS_CODONS = {
    'Phe': ['UUU','UUC'],
    'Leu': ['UUA','UUG','CUU','CUC','CUA','CUG'],
    'Ile': ['AUU','AUC','AUA'],
    'Met': ['AUG'],
    'Val': ['GUU','GUC','GUA','GUG'],
    'Ser': ['UCU','UCC','UCA','UCG','AGU','AGC'],
    'Pro': ['CCU','CCC','CCA','CCG'],
    'Thr': ['ACU','ACC','ACA','ACG'],
    'Ala': ['GCU','GCC','GCA','GCG'],
    'Tyr': ['UAU','UAC'],
    'STOP': ['UAA','UAG','UGA'],
    'His': ['CAU','CAC'],
    'Gln': ['CAA','CAG'],
    'Asn': ['AAU','AAC'],
    'Lys': ['AAA','AAG'],
    'Asp': ['GAU','GAC'],
    'Glu': ['GAA','GAG'],
    'Cys': ['UGU','UGC'],
    'Trp': ['UGG'],
    'Arg': ['CGU','CGC','CGA','CGG','AGA','AGG'],
    'Gly': ['GGU','GGC','GGA','GGG']
    }

class CDSConcatenate(DNASequence):
    def __init__(self, header, string):
        super().__init__(header, string, spades = False)
        self.codon_counts = pd.Series(
            {
            'UUU': 0, 'UUC': 0, 'UUA': 0, 'UUG': 0,
            'CUU': 0, 'CUC': 0, 'CUA': 0, 'CUG': 0,
            'AUU': 0, 'AUC': 0, 'AUA': 0, 'AUG': 0,
            'GUU': 0, 'GUC': 0, 'GUA': 0, 'GUG': 0,
            'UCU': 0, 'UCC': 0, 'UCA': 0, 'UCG': 0,
            'AGU': 0, 'AGC': 0, 'CCU': 0, 'CCC': 0,
            'CCA': 0, 'CCG': 0, 'ACU': 0, 'ACC': 0,
            'ACA': 0, 'ACG': 0, 'GCU': 0, 'GCC': 0,
            'GCA': 0, 'GCG': 0, 'UAU': 0, 'UAC': 0,
            'UAA': 0, 'UAG': 0, 'UGA': 0, 'CAU': 0,
            'CAC': 0, 'CAA': 0, 'CAG': 0, 'AAU': 0,
            'AAC': 0, 'AAA': 0, 'AAG': 0, 'GAU': 0,
            'GAC': 0, 'GAA': 0, 'GAG': 0, 'UGU': 0,
            'UGC': 0, 'UGG': 0, 'CGU': 0, 'CGC': 0,
            'CGA': 0, 'CGG': 0, 'AGA': 0, 'AGG': 0,
            'GGU': 0, 'GGC': 0, 'GGA': 0, 'GGG': 0
            }
        )
        self.rscu_table = {}
    def split_codons(self):
        return [ complement(self.transcribe()[i:i+3]) for i in range(0, self.length, 3) ]
    
    def count_codons(self):
        for codon in self.split_codons():
            if 'N' in codon:
                continue
            self.codon_counts[codon] += 1
    
    def calculate_rscu(self):
        for amino_acid, codons in SYNONYMOUS_CODONS.items():
            synonymous_codon_counts = self.codon_counts[codons]
            amino_acid_occurences = sum(synonymous_codon_counts)
            for c in codons:
                if amino_acid_occurences == 0:
                    rscu = 0.00
                else:
                    pass
                print (self.header, amino_acid, c, self.codon_counts[c])

class Codons(Module, LoggingEntity, Head):
    def __init__(self, argdict = None):
        super().__init__(self.__class__)
        if argdict is not None:
            self.config.load_argdict(argdict)
        else:
            self.argparser = self.generate_argparser()
            self.parsed_args = self.argparser.parse_args()
            self.config.load_cmdline( self.parsed_args ) # Copy command line args defined by self.argparser to self.config
        
        self.config.reusable.populate(
                ReusableOutput(
                    arg = "gff3",
                    pattern = ".*[.]aug[.]out[.]gff3$",
                    genfunc = augustus_predict,
                    genfunc_args = {
                        "prefix": self.config.get("prefix"),
                        "nucl": self.config.get("nucl"),
                        "augustus_sp": self.config.get("augustus_sp"),
                        "outpath": f"{self.config.get('prefix')}.aug.out.gff3"
                        }
                    )
                )
        self.config.dependencies.populate(
                CaseDependency("augustus", "gff3", None),
            )
        
        if self.config.get("mode") == "blastp":
            self.config.reusable.add(
                ReusableOutput(
                    arg = "blastout",
                    pattern = ".*[.]spdb[.]blast[.]out$",
                    genfunc = protein_blast,
                    genfunc_args = {
                        "prot_path": self.config.get("prot"),
                        "db": self.config.get("spdb"),
                        "evalue": self.config.get("evalue"),
                        "cpus": self.config.get("cpus"),
                        "outpath": f"{self.config.get('prefix')}.spdb.blast.out"
                        }
                    )
            )
            self.config.dependencies.add(
                CaseDependency("blastp", "blastout", None)
            )

        elif self.config.get("mode") == "blastn":
            self.config.reusable.add(
                ReusableOutput(
                    arg = "blastout",
                    pattern = ".*[.]nt[.]blast[.]out$",
                    genfunc = nucleotide_blast,
                    genfunc_args = {
                        "nucl_path": self.config.get("nucl"),
                        "db": "nt",
                        "evalue": self.config.get("evalue"),
                        "cpus": self.config.get("cpus"),
                        "outpath": f"{self.config.get('prefix')}.nt.blast.out"
                        }
                    )
            )
            self.config.dependencies.add(
                CaseDependency("blastn", "blastout", None)
            )
        
        else:
            print("Bad mode selection.")
            sys.exit(1)
            
    
    def generate_argparser(self):
        parser = argparse.ArgumentParser()
        parser.add_argument("mod", nargs="*")
        parser.add_argument('-m','--gff3', metavar='gene_models', action = PathAction, required = False, default = None, help ="A gff3 file from Augustus (one is generated by scgid blob) that contains the gene models for your metagenome.")
        parser.add_argument('-n','--nucl', metavar='contig_fasta', action = PathAction, required = True, help ="The contig fasta associated with your metagenome.")
        parser.add_argument('-p','--prot', metavar = "protein_fasta", action=PathAction, required=False, help = "A FASTA file containing the proteins called from the genome.")

        parser.add_argument('-g', '--targets', metavar = 'target_taxa', action='store', required=True, help="A comma-separated list with NO spaces of the taxonomic levels that the gc-coverage window should be chosen with respect to including. EXAMPLE: '-g Fungi,Eukaryota,Homo'")
        parser.add_argument('-x', '--exceptions', metavar = 'exceptions_to_target_taxa', action='store', required=False, default=None, help="A comma-separated list with NO spaces of any exlusions to the taxonomic levels specified in -g|--targets. For instance if you included Fungi in targets but want to exclude ascomycetes use: '-x Ascomycota'")
        parser.add_argument('-f','--prefix', metavar = 'prefix_for_output', required=False, default='scgid', help="The prefix that you would like to be used for all output files. DEFAULT = scgid")
        parser.add_argument('--cpus', metavar = 'cores', action = 'store', required = False, default = "1", help = "The number of cores available for BLAST to use.")
        parser.add_argument('--mode', metavar = "mode", action="store",required=False, default ='blastp', help = "The type of blast results that you would like to use to annotate the tips of the RSCU tree ('blastp' or 'blastn'). This module will automatically do a blastn search of the NCBI nt database for you. At this time, a blastp search can not be run directly from this script. INSTEAD, if using mode 'blastp' (DEFAULT, recommended) you must specify a scgid blob-derived _info_table.tsv file with -i|--infotable")
        parser.add_argument('--minlen', metavar = 'minlen', action = 'store', required = False, default = '3000', help = 'Minimum length of CDS concatenate to be kept and used to build RSCU tree. Highly fragmented assemblies will need this to be reduced. Reduce in response to `Tree too small.` error.')
        parser.add_argument('-sp','--augustus_sp', metavar = "augustus_species", action="store",required=False, default=None, help = "Augustus species for gene predicition. Type `augustus --species=help` for list of available species designations.")
        parser.add_argument('-e', '--evalue', metavar = 'e-value_cutoff', action = 'store', required = False, default = '1e-5', help = "The evalue cutoff for blast. Default: 1xe-5)")
        parser.add_argument('-b','--blastout', metavar = "blastout", action=PathAction, required=False, help = "The blast output file from a blastn search of the NCBI nt database with your contigs as query. If you have not done this yet, this script will do it for you.")
        parser.add_argument('-i','--infotable', metavar = "infotable", action=PathAction, required=False, help = "The scgid gc-cov-derived infotable generated by a blastp search of a swissprot-style protein database.")
        parser.add_argument("--noplot", action="store_true", default=False, required=False, help="Turns of plotting of annotated trees to PDF.")
        parser.add_argument('--Xmx', metavar = "available_memory", action="store",required=False, default = "2g", help = "Set memoray available to run ClaMs. Specicy as such: X megabytes = Xm, X gigabytes = Xg")

        return parser

    def locate_cds_gff3 (self, gff3_path):

        contig_chunks = {}

        with open (gff3_path, 'r') as gff3:

            CDSChunk = namedtuple("CDSChunk", ["start", "end", "strand"])
            for line in gff3:

                # Skip comment lines
                if line[0] == "#":
                    continue
                
                spl = line.split('\t')

                # Ignore all but CDS lines in gff3
                if spl[2] == "CDS":

                    shortname = '_'.join( spl[0].split('_')[0:2] )

                    # Capture pid     
                    s = re.search("[.](g[0-9]+)[.]",spl[8])
                    pid = s.group(1)

                    # Group CDS lines in gff3 by parent contig (by shortname) and protein (by pid)
                    if shortname in contig_chunks:
                        if pid in contig_chunks[shortname]:
                            contig_chunks[shortname][pid].append( CDSChunk(spl[3], spl[4], spl[6]) )
                        else:
                            contig_chunks[shortname][pid] = [ CDSChunk(spl[3], spl[4], spl[6]) ]
                    else:
                        contig_chunks[shortname] = {
                                pid: [ CDSChunk(spl[3], spl[4], spl[6]) ]
                                }
        return contig_chunks

    def concatenate_cds (self, contig_chunks, nucl):
        cds_concatenates = {}

        # Iterate through CDS chunks of predicted proteins on each contig and pull CDS sequences from nucleotide fasta
        for shortname, pids in contig_chunks.items():

            contig_cds_cat = str()
            
            for _, chunks in pids.items():

                gene_cds = str()

                for chunk in chunks:
                    
                    # Ignore zero-length CDS chunks
                    if chunk.start == chunk.end:
                        continue
                    
                    # Fetch CDS sequence from contig by start/stop indices listed in gff3 and revcomp if on reverse strand
                    if chunk.strand == '-':
                        chunk_string = revcomp(nucl.index[shortname].string[
                            int(chunk.start)-1: int(chunk.end)
                            ])
                    else: # i.e., chunk.strand == '+'
                        chunk_string = nucl.index[shortname].string[
                            int(chunk.start)-1: int(chunk.end)
                            ]

                    # Combine chunk sequences for each gene
                    gene_cds += chunk_string

                # Toss out predicted CDS if they aren't divisible by 3 to avoid introducing frameshifts into CDS concatenate
                # Combine gene_cds into contig_cds_cat beacuse they occur on the same contig
                if len(gene_cds) % 3 == 0:
                    contig_cds_cat += gene_cds
            
            # Store contig_cds_cat in DNASequence object and add to dict
            if len(contig_cds_cat) != 0:
                cds_concatenates[shortname] = CDSConcatenate(shortname, contig_cds_cat)
        
        # Return all contig-level CDS concatenates as a DNASequenceCollection object
        return DNASequenceCollection().from_dict(cds_concatenates)

    def run(self):
        self.start_logging()
        self.setwd( __name__, self.config.get("prefix") )
        self.config.reusable.check()
        self.config.dependencies.check(self.config)

        self.logger.info(f"Running in {self.config.get('mode')} mode.")

        # Read in nucleotide FASTA
        nucl = DNASequenceCollection().from_fasta(self.config.get("nucl"))

        # Rekey nucl by shortname
        nucl.rekey_by_shortname()

        # Concatenate all CDS sequences on each contig
        cds_coords = self.locate_cds_gff3(
            self.config.get("gff3")
        )
        for k,v in cds_coords.items():
            print(k,v)
        
        #sys.exit()

        cds_concatenates = self.concatenate_cds(
            cds_coords,
            nucl
        )
        
        # Remove CDS concatenates shorter than supplied minlin
        cds_concatenates.remove_small_sequences( int(self.config.get("minlen")) )

        cds_concatenates.write_fasta("large_concat.fasta")

        for c in cds_concatenates.seqs():
            c.count_codons()
            c.calculate_rscu()