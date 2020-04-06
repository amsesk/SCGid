import argparse
import sys
import os
import re
import logging
import pandas as pd
import numpy as np
from collections import namedtuple
from ete3 import Tree, NCBITaxa
from scgid.module import Module
from scgid.modcomm import LoggingEntity, Head, get_head, ErrorHandler
from scgid.reuse import ReusableOutput, ReusableOutputManager, augustus_predict, nucleotide_blast, protein_blast
from scgid.dependencies import CaseDependency
from scgid.parsers import PathStore, StoreInt
from scgid.sequence import DNASequenceCollection, DNASequence, revcomp, complement
from scgid.library import subprocessP
from scgid.infotable import InfoTable, get_by_idx, count_unique
from scgid.error import ModuleError, ArgumentError, Ok, check_result

class SmallTreeError(ModuleError):
    def __init__(self, ntips, mincladesize, minlen, error_catch = True):
        super().__init__()
        self.msg = f"The RSCU tree is smaller ({ntips} tips) than -c|--mincladesize so a best clade cannot be chosen. You can try reducing -c|--mincladesize ({mincladesize} tips) or --minlen ({minlen} bp) to address."
        if error_catch:
            self.catch()

class NoGoodCladesError(ModuleError):
    def __init__(self, error_catch = True):
        super().__init__()
        self.msg = f"All clades of sufficient size contain more nontarget-annotated tips than target-annotated tips. Consider expanding the SPDB with `scgid spexpand` or changing the target designation supplied to `-g|--target`."
        if error_catch:
            self.catch()

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
        self.rscu_table = pd.Series(
            {
            'UUU': 0.00, 'UUC': 0.00, 'UUA': 0.00, 'UUG': 0.00,
            'CUU': 0.00, 'CUC': 0.00, 'CUA': 0.00, 'CUG': 0.00,
            'AUU': 0.00, 'AUC': 0.00, 'AUA': 0.00, 'AUG': 0.00,
            'GUU': 0.00, 'GUC': 0.00, 'GUA': 0.00, 'GUG': 0.00,
            'UCU': 0.00, 'UCC': 0.00, 'UCA': 0.00, 'UCG': 0.00,
            'AGU': 0.00, 'AGC': 0.00, 'CCU': 0.00, 'CCC': 0.00,
            'CCA': 0.00, 'CCG': 0.00, 'ACU': 0.00, 'ACC': 0.00,
            'ACA': 0.00, 'ACG': 0.00, 'GCU': 0.00, 'GCC': 0.00,
            'GCA': 0.00, 'GCG': 0.00, 'UAU': 0.00, 'UAC': 0.00,
            'UAA': 0.00, 'UAG': 0.00, 'UGA': 0.00, 'CAU': 0.00,
            'CAC': 0.00, 'CAA': 0.00, 'CAG': 0.00, 'AAU': 0.00,
            'AAC': 0.00, 'AAA': 0.00, 'AAG': 0.00, 'GAU': 0.00,
            'GAC': 0.00, 'GAA': 0.00, 'GAG': 0.00, 'UGU': 0.00,
            'UGC': 0.00, 'UGG': 0.00, 'CGU': 0.00, 'CGC': 0.00,
            'CGA': 0.00, 'CGG': 0.00, 'AGA': 0.00, 'AGG': 0.00,
            'GGU': 0.00, 'GGC': 0.00, 'GGA': 0.00, 'GGG': 0.00
            }
        )
    def split_codons(self):
        self.transcomp = complement(self.transcribe())
        return [ self.transcomp[i:i+3] for i in range(0, self.length, 3) ]

    def count_codons(self):
        for codon in self.split_codons():
            if 'N' in codon:
                continue
            self.codon_counts[codon] += 1
        logrow = "[{:<10}] ".format(self.header) + " ".join(["{:<5}".format(x) for x in self.codon_counts.to_list()])
        return logrow

    def calculate_rscu(self):
        for _, codons in SYNONYMOUS_CODONS.items():
            synonymous_codon_counts = self.codon_counts[codons]
            amino_acid_occurences = sum(synonymous_codon_counts)
            for c in codons:
                if amino_acid_occurences == 0:
                    rscu = 0.00
                else:
                    rscu = self.codon_counts[c] / (amino_acid_occurences * (1/len(self.codon_counts[codons])))

                self.rscu_table[c] = rscu

class RSCUTree(object):
    def __init__(self, treepath, head = None):
        self.dendrogram = Tree(treepath)
        
        if head is None:
            self.head = get_head()
        else:
            self.head = head

        self.mode = self.head.config.get("mode")

        self.ntips = len(self.dendrogram.get_leaves())
        self.nnodes = len(self.dendrogram.get_descendants())
        self.ninternalnodes = self.nnodes - self.ntips

        if self.mode == "blastp":
            self.infotable = InfoTable()
            self.infotable.set_target(
                self.head.config.get("targets"),
                self.head.config.get("exceptions")
            )
            self.infotable.load(self.head.config.get("infotable"))
            self.infotable.parse_lineage()
            
            self.best_clade = None

    def annotate(self):
        self.infotable.decide_inclusion()
        if self.mode == "blastp":
            self._protannot()

        elif self.mode == "blastn":
            raise NotImplementedError("Annotation with NCBI nt blast results not yet implemented.")

        else:
            raise NotImplementedError(f"Bad mode selection: `{self.mode}`.")

        return 0

    def _protannot(self):
       for n in self.dendrogram.traverse():
            if n.is_leaf():
                n.add_feature("annotation","unclassified")
                if n.name in self.infotable.keep:
                    n.annotation = "target"
                elif n.name in self.infotable.dump:
                    n.annotation = "nontarget"
                else:
                    pass # All leaves set to unclassified above, so no need to do anything

    def _nuclannot(self):
        pass

    def pick_clade(self, error_catch = True):
        best = [] #"best" here means worth looking into, although it is going to include some real shitty trees

        # Include self.dendrogram in list of clades to encompass whole tree, but warn if the entire tree is selected as the best clade (warned at base of function)

        clades_of_sufficient_size = [c for c in [self.dendrogram]+self.dendrogram.get_descendants() if len(c.get_leaves()) >= int(self.head.config.get("mincladesize")) and not c.is_leaf()]
        if len(clades_of_sufficient_size) == 0:

            return SmallTreeError(
                ntips = self.ntips,
                mincladesize = self.head.config.get("mincladesize"),
                minlen = self.head.config.get("minlen"),
                error_catch = error_catch)


        for n in clades_of_sufficient_size:
            count_t = float(len([l for l in n if l.annotation == "target"]))
            count_nt = float(len([l for l in n if l.annotation == "nontarget"]))
            count_unclass = float(len([l for l in n if l.annotation == "unclassified"]))

            total_classified = float(count_t + count_nt)

            if total_classified == 0:
                continue

            measure = float( (count_t - count_nt) / total_classified )
            n.add_feature("measure",measure)
            n.add_feature("leaves", len(n.get_leaves()))
            best.append(n)

        #Remove clades with more nontarget than target
        best = [n for n in best if n.measure > 0.00]

        if len(best) == 0:
            return NoGoodCladesError(error_catch = error_catch)

        # Order clades by measure, descending
        best = sorted(best, key=lambda x: x.measure, reverse=True)

        #Bin trees based on common edges, or in other words look for sets of trees that aren't just nested within eachother from the list of best trees
        #However, based on how it's written now, if a monophyletic clade is nested within another, some (likely) poorer quality clades will include both and lead to a node ending up in multiple bins, BUT hopefully selection of the best clade from each bin sorts this out
        bins = {}
        i=0 #index iterator to add additional bins
        for p in best:
            if p == self.dendrogram:
                continue
            n_bins = len(bins) #number of bins present on each iteration
            found_match = False
            for n in range(0,n_bins):
                for t in bins[n]:
                    if len(p.compare(t)['common_edges']) != 0:
                        bins[n].append(p)
                        found_match = True
                        break
            if not found_match:
                bins[i] = [p]
                i+=1
        bins[i] = [self.dendrogram]
        #Now go through each bin and find the best tree from that bin
        best_trees_by_bin = []
        for b in range(0,len(bins)):
            measures = [i.measure for i in bins[b]]
            leaves = [i.leaves for i in bins[b]]
            max_measure_indices = [i for i,v in enumerate(measures) if v == max(measures)]
            if len(max_measure_indices) == 1:
                best_idx = measures.index(max(measures))
                best_trees_by_bin.append(bins[b][best_idx])
            else:
                leaves_of_maxes = [v for i,v in enumerate(leaves) if i in max_measure_indices]
                max_leaves = max(leaves_of_maxes)
                best_idx = leaves.index(max_leaves)
                best_trees_by_bin.append(bins[b][best_idx])

        #for i in best_trees_by_bin:
        #    print i.leaves,"\t",i.measure
        #Now compare the best trees from each bin to decide on a best tree (or set of best trees) for subsequent actions

        max_measure = max([i.measure for i in best_trees_by_bin])
        best_trees = [i for i in best_trees_by_bin if i.measure == max_measure]
        if len(best_trees) > 1:
            max_leaves = max([i.leaves for i in best_trees])
            best_trees = [i for i in best_trees if i.leaves == max_leaves]
            if len(best_trees) > 1:
                self.head.logger.warning("More than one best tree detemined from codon analysis.")
            else:
                final_tree = best_trees[0]
            #for tree in best_trees:
                #tree.show(tree_style=circ)
        else:
            final_tree = best_trees[0]

        self.best_clade = final_tree

        if self.best_clade == self.dendrogram:
            self.head.logger.warning("The entire RSCU tree was selected as the best clade. This either means: (i) there are no nontarget tips in your tree or (ii) your RSCU tree is small with few or no clades larger than -c|--mincladesize. Try raising -c|--mincladesize or --minlen.")

        return Ok(None)

    def write_tree_annotation(self, trainset, outpath):
        taxlvl = self.infotable.taxon_level(level=1)

        taxlvl = taxlvl.groupby("contig").agg({ "lineage": lambda x: ','.join(x).split(','),
            "evalue": lambda x: [e for e in x]
            })
        ### get all evalues equal to best evalues, and their indices
        taxlvl['maxes'] = taxlvl.evalue.apply(lambda x: [i for i,e in enumerate(x) if e == min(x)])

        ### match hits with the evalue for that hit
        taxlvl['maxtax'] = taxlvl.apply(get_by_idx, axis=1)

        ### Get best taxonomy based on counts of best evalue taxa
        taxlvl['decide'] = taxlvl['maxtax'].apply(count_unique)
        taxlvl_decisions = taxlvl[["decide"]]

        with open(outpath, 'w') as f:
            for leaf in self.dendrogram.iter_leaves():

                try:
                    taxon = taxlvl_decisions.loc[ leaf.name ].values[0]
                except KeyError:
                    taxon = "unclassified"

                if leaf.name in trainset.index.keys():
                    in_trainset = "trainset"
                else:
                    in_trainset = "not_selected"


                f.write(f"{leaf.name},{taxon},{leaf.annotation},{in_trainset}\n")


class Codons(Module, LoggingEntity, Head, ErrorHandler):
    def __init__(self, argdict = None, loglevel = logging.INFO):
        super().__init__(self.__class__, loglevel=loglevel)
        if argdict is not None:
            self.translated_args = self.translate_argdict(argdict, Codons.generate_argparser())
            self.config.load_argdict(self.translated_args)
            #print(self.config)
            self.parsed_args = self.config
        else:
            self.argparser = Codons.generate_argparser()
            self.parsed_args = self.argparser.parse_args()
            self.config.load_cmdline( self.parsed_args ) # Copy command line args defined by self.argparser to self.config

        self.set_rundir(self.config.get("prefix"))

        ''' This could be useful later for converting config back to argdict
        for l in self.config.__repr__().split("\n"):
            elements = [x.strip() for x in l.split(":")]
            print (f"\"{elements[0]}\": {elements[1]},")
        sys.exit()
        '''

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
                        "prot": self.config.get("prot"),
                        "db": self.config.get("spdb"),
                        "evalue": self.config.get("evalue"),
                        "cpus": self.config.get("cpus"),
                        "outpath": f"{self.config.get('prefix')}.spdb.blast.out"
                        }
                    )
                )

            self.config.reusable.add(
                ReusableOutput(
                    arg = "infotable",
                    pattern = f"{self.config.get('prefix')}[.]infotable[.]tsv$",
                    genfunc = ArgumentError,
                    genfunc_args = {
                        "msg": "Option -i|--infotable required to run in blastp mode. Run `scgid gct` to generate one or set `--mode blastn`."
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

        self.keep = {}


    def generate_argparser():
        parser = argparse.ArgumentParser()
        parser.add_argument("mod", nargs="*")
        parser.add_argument('-m','--gff3', dest="gff3", metavar='gene_models', action = PathStore, required = False, default = None, help ="A gff3 file from Augustus (one is generated by scgid blob) that contains the gene models for your metagenome.")
        parser.add_argument('-n','--nucl', metavar='assembly_fasta', action = PathStore, required = True, help ="The contig fasta associated with your metagenome.")
        parser.add_argument('-p','--prot', metavar = "protein_fasta", action=PathStore, required=False, help = "A FASTA file containing the proteins called from the genome.")

        parser.add_argument('-g', '--targets', metavar = 'target_taxa', action='store', required=True, help="A comma-separated list with NO spaces of the taxonomic levels that the gc-coverage window should be chosen with respect to including. EXAMPLE: '-g Fungi,Eukaryota,Homo'")
        parser.add_argument('-x', '--exceptions', metavar = 'exception_taxa', action='store', required=False, default=None, help="A comma-separated list with NO spaces of any exlusions to the taxonomic levels specified in -g|--targets. For instance if you included Fungi in targets but want to exclude ascomycetes use: '-x Ascomycota'")
        parser.add_argument('-f','--prefix', metavar = 'output_prefix', required=False, default='scgid', help="The prefix that you would like to be used for all output files. DEFAULT = scgid")
        parser.add_argument('--cpus', metavar = 'cpus', action = 'store', required = False, default = "1", help = "The number of cores available for BLAST to use.")
        parser.add_argument("-am", '--mode', metavar = "annotation_mode", action="store",required=False, choices=["blastp", "blastn"], default ='blastp', help = "The type of blast results that you would like to use to annotate the tips of the RSCU tree ('blastp' or 'blastn'). This module will automatically do a blastn search of the NCBI nt database for you. At this time, a blastp search can not be run directly from this script. INSTEAD, if using mode 'blastp' (DEFAULT, recommended) you must specify a scgid blob-derived _info_table.tsv file with -i|--infotable")
        parser.add_argument("-l", '--minlen', metavar = 'minlen', action = StoreInt, required = False, default = 3000, help = 'Minimum length of CDS concatenate to be kept and used to build RSCU tree. Highly fragmented assemblies will need this to be reduced. Reduce in response to `Tree too small.` error.')
        parser.add_argument('-c', '--mincladesize', metavar = 'mincladesize', action = StoreInt, required = False, default = 30, help = 'Minimum size of clade to serve as training set for ClaMS.')
        parser.add_argument('-sp','--augustus_sp', metavar = "augustus_species", action="store",required=False, default=None, help = "Augustus species for gene predicition. Type `augustus --species=help` for list of available species designations.")
        parser.add_argument('-e', '--evalue', metavar = 'blast_evalue_cutoff', action = 'store', required = False, default = '1e-5', help = "The evalue cutoff for blast. Default: 1xe-5)")
        parser.add_argument('-b','--blastout', metavar = "spdb_blast_output", action=PathStore, required=False, help = "The blast output file from a blastn search of the NCBI nt database with your contigs as query. If you have not done this yet, this script will do it for you.")
        parser.add_argument('-i','--infotable', metavar = "infotable", action=PathStore, required=False, help = "The scgid gc-cov-derived infotable generated by a blastp search of a swissprot-style protein database.")
        parser.add_argument("-np", "--noplot", action="store_true", default=False, required=False, help="Turns of plotting of annotated trees to PDF.")
        parser.add_argument('--Xmx', metavar = "available_memory", action="store",required=False, default = "2g", help = "Set memoray available to run ClaMs. Specicy as such: X megabytes = Xm, X gigabytes = Xg")

        parser.add_argument('-db', '--spdb', metavar = 'swissprot_fasta', action=PathStore, required=False, default=None,  help = "The path to your version of the swissprot database in FASTA format.")
        parser.add_argument('-t','--taxdb', metavar = "swissprot_taxdb", action=PathStore, required=False, default=None, help = "The location of the taxonomy database, likely provided by an earlier script.")

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

    # Investigate inclusion of STOP codons in CDS here later - see how the trees look first excluding them?
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

    # Function that confirms size of CDSConcatenate DNASequenceCollection is larger than minlen arg
    @check_result
    def check_n_concatenates (self, seq_collection, error_catch = True):
        if not seq_collection.check_size(int(self.config.get("mincladesize"))):
            
            return SmallTreeError (
                ntips = len(seq_collection.seqs()), 
                mincladesize = self.config.get("mincladesize"),
                minlen = self.config.get("minlen"),
                error_catch = error_catch
                )
        
        else:
            return Ok(None)

    def rscu_distances(self, cds_concatenates):
        headers = [s.header for s in cds_concatenates.seqs()]
        matrix = pd.DataFrame(columns = headers, index = headers)
        for p in cds_concatenates.seqs():
            row = {}
            for r in cds_concatenates.seqs():
                row[r.header] = (1/59)*sum( abs( p.rscu_table.to_numpy() - r.rscu_table.to_numpy() ) )
            matrix.loc[p.header] = row

        matrix = np.tril(matrix)

        matrix_frame = pd.DataFrame(matrix, columns = headers, index = headers)
        #matrix_frame.iloc[np.triu_indices(n = len(headers), k=0)] = np.nan

        return matrix_frame

    def write_nexus(self, distance_matrix, outpath):

        header = f"#NEXUS\nbegin distances;\ndimensions\nntax={distance_matrix.shape[0]};\nmatrix\n"
        trailer = ";\nEND"

        with open(outpath,'w') as f:
            f.write(header)

        distance_matrix.to_csv(outpath, sep=' ',index = True, header = False, mode = 'a')

        with open(outpath,'a') as f:
            f.write(trailer)

    def nj_tree(self, distmat_csv, outpath):
        cmd = [
            self.config.get("path_to_Rscript"),
            "--vanilla",
            os.path.join(self.config.SCGID_SCRIPTS,"ape_nj.R"),
            distmat_csv,
            outpath
            ]
        self.logger.info(' '.join(cmd))
        subprocessP(cmd, self.logger)

    def draw_annotated_tree(self, treefile, annotfile):
        cmd = [
            self.config.get("path_to_Rscript"),
            "--vanilla",
            os.path.join(self.config.SCGID_SCRIPTS, "codons_phytools.R"),
            treefile,
            annotfile,
            f"{self.config.get('prefix')}.tree.annotated.pdf"
            ]
        self.logger.info(' '.join(cmd))
        subprocessP(cmd, self.logger)

    def run_clams(self):
        prefix = self.config.get('prefix')
        with open(f"{prefix}.clams.trainset.tsv",'w') as f:
            f.write("rscu_derived_ts1\t")
            f.write(f"{prefix}.trainset.fasta")
            f.write("\n")

        self.logger.info("Running ClaMs...")
        cmd = [
            "java",
            f"-Xmx{self.config.get('Xmx')}",
            "-jar", os.path.join(self.config.get("clams_path")),
            os.path.join(f"{prefix}.clams.trainset.tsv"),
            self.config.get("nucl"),
            f"{prefix}.clams.out",
            "DBC", "2",
            "0.01217181"
            ]
        self.logger.info(' '.join(cmd))
        subprocessP(cmd, [self.logger, self.simplelogger], log_stdout=True)

        return None

    def parse_clams_out(self) -> list:
        to_keep = []
        count = 0
        with open(f"{self.config.get('prefix')}.clams.out") as clams_out:
            for line in clams_out.readlines():
                spl = [x.strip() for x in line.split("\t")]
                if spl[1] == "rscu_derived_ts1":
                    count += 1
                    to_keep.append("_".join(spl[0].split("_")[0:2]))
        self.logger.info(f"ClaMs finished, {len(to_keep)} matches to trainset detected.")
        return to_keep

    @check_result
    def run(self) -> DNASequenceCollection:
        #self.start_logging()
        self.setwd( __name__ )

        ##############################################
        ######## Skip this if called directly ########
        ######## (i.e., in tests)             ########
        ##############################################
        if self.root is not None:
            self.log_to_rundir(type(self).__name__)
            self.log_config()
        
        self.config.reusable.check()
        self.config.dependencies.check(self.config)
        self.config.reusable.generate_outputs()

        self.logger.info(f"Running in {self.config.get('mode')} mode.")

        # Read in nucleotide FASTA
        nucl = DNASequenceCollection().from_fasta(self.config.get("nucl"))
        self.logger.info(f"Read nucleotide FASTA at `{self.config.get('nucl')}` into memory.")

        # Rekey nucl by shortname
        nucl.rekey_by_shortname()

        ####################################
        #'''
        # Get coordinates and strand of CDS chunks
        self.logger.info(f"Pulling locations of CDS chunks from gff3 at `{self.config.get('gff3')}`")
        cds_coords = self.locate_cds_gff3(
            self.config.get("gff3")
        )

        # Concatenate all CDS sequences on each contig
        self.logger.info(f"Concatenating CDS chunks from each contig")
        cds_concatenates = self.concatenate_cds(
            cds_coords,
            nucl
        )

        # Remove CDS concatenates shorter than supplied minlin
        self.logger.info(f"Removing CDS concatenates <{self.config.get('minlen')} bp in length")
        cds_concatenates.remove_small_sequences( int(self.config.get("minlen")) )

        # Check that the DNASequenceCollection of CDSConcatenate is large enough compared to self.config.minlen
        self.check_n_concatenates(cds_concatenates)
        
        # Write CDS concatenates of sufficient length to fasta
        cds_concatenates.write_fasta(f"cds_concatenates_larger_than{self.config.get('minlen')}bp.fasta")

        # Count codons on each CDS concatenate
        self.logger.info(f"Counting codon occurences on each CDS concatenate")
        for c in cds_concatenates.seqs():
            self.simplelogger.info( c.count_codons() )

        # Remove CDS concatenates that contain STOP codons, becaause they shouldn't... I don't think
        self.logger.info(f"Removing CDS concatenates that contain STOP codons, because they shouldn\'t")
        cds_concatenates = DNASequenceCollection().from_dict(
            { header: seqobj for header, seqobj in cds_concatenates.index.items() if sum( seqobj.codon_counts[ SYNONYMOUS_CODONS["STOP"] ] ) == 0 }
        )

        # Calculate RSCU for each codon on each concatenate
        self.logger.info(f"Calculating contig RSCU profiles")
        for c in cds_concatenates.seqs():
            c.calculate_rscu()

        # Compute RSCU distance matrix
        self.logger.info(f"Constructing RSCU distance matrix")
        distance_matrix = self.rscu_distances(cds_concatenates)

        distance_matrix.to_csv("../distmat.csv", sep=",")
        #'''
        ###########################
        distance_matrix = pd.read_csv("../distmat.csv", sep=",", index_col=0)
        #print (distance_matrix.head)

        # Write RSCU distance matrix to NEXUS format and CSV format (for R)
        self.write_nexus(distance_matrix, f"{self.config.get('prefix')}.rscu.distance.matrix.nex")
        distance_matrix.to_csv(f"{self.config.get('prefix')}.rscu.distance.matrix.csv", sep=',', index = True, header = False, mode = 'w')
        self.logger.info(f"Printed RSCU distance matrix in NEXUS format to `{self.config.get('prefix')}.rscu.distance.matrix.csv`")

        # Compute NJ tree and write to .tre file (use R because biopython reads big trees really slow)
        self.logger.info(f"Building NJ tree from RSCU distance matrix in R")
        treefile = f"{self.config.get('prefix')}_rscu_nj.tre"
        self.nj_tree(
            f"{self.config.get('prefix')}.rscu.distance.matrix.csv",
            treefile
        )

        # Instantiate RSCUTree object and pass it NJ tree generated by R
        nj_tree = RSCUTree(treefile)

        # Annotate tips of the dendrogram. Mode is handled inside function body.
        nj_tree.annotate()

        # Pick the best clade for ClaMS training by iteratively binning and comparing clades' target:nontarget ratios
        nj_tree.pick_clade()

        # Return filtered DNASequenceCollection with non-trainset contigs excluded
        trainset = nucl.header_list_filter(nj_tree.best_clade.get_leaf_names())

        # Write trainset to FASTA
        trainset.write_fasta(f"{self.config.get('prefix')}.trainset.fasta")

        # Write tree annotation file for plotting with phytools in R
        nj_tree.write_tree_annotation(
            trainset,
            f"{self.config.get('prefix')}_njtree_annotations.csv"
            )

        # Draw annotated tree with phytools in R
        '''
        self.draw_annotated_tree(
            treefile,
            f"{self.config.get('prefix')}_njtree_annotations.csv"
        )
        '''

        self.run_clams()
        to_keep = self.parse_clams_out()

        # Construct DNASequenceCollection from final filtered assembly,
        # Resorted by DNASequenceCollection.header_list_filter()
        final_assembly = nucl.header_list_filter(to_keep)

         # Compute final filtered assembly stats
        filtered_size = sum([len(s.string) for s in final_assembly.seqs()])
        filtered_ncontigs = len(final_assembly.seqs())

        self.logger.info(f"Filtered assembly contains {filtered_ncontigs:,} contigs with a cumulative size of {filtered_size:,} bp ({filtered_size/1e6:.2f} Mbp).")

        # Print final filtered assembly to FASTA
        final_fname = f"{self.config.get('prefix')}.codons.filtered.assembly.fasta"
        final_assembly.write_fasta( final_fname )

        self.logger.info(f"Final filtered assembly written in FASTA format to `{final_fname}`")

        self.logger.info("RSCU-based filtering complete. Returning to SCGid.")

        # Migrate and then remove temp dir, cd back to starting dir
        self.migrate_temp_dir()
        self.resetwd()

        # Return final filtered assembly to SCGid root
        return Ok(final_assembly)