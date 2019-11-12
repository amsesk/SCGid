import logging
import inspect
import sys
import re
import os
from Depends import CaseDependency, ConstDependency
from Module import Module, Config
import argparse

class NewModule(Module):
    def __init__(self):
        super(NewModule, self).__init__(self.__class__)
        frm = inspect.stack()[1]
        self.parent = inspect.getmodule(frm[0])
        self.log = logging.getLogger("{}.{}".format(self.parent.__name__, __name__))
        self.argparser = self.generate_argparser()
        self.parsed_args = self.argparser.parse_args()
        self.config.load_cmdline( self.parsed_args )
        self.locate_input_files(
            [".*[.]spdb[.]",
            ".*[.]aug[.]out[.]fasta",
            ".*[.]aug[.]out[.]gff3",
            ".*[.]spdb[.]blast[.]out"]
            )
        '''
        self.config.dependencies.populate(
            CaseDependency("blastp", "blastout", None),
            CaseDependency("augustus", "prot", None),
        )
        '''
        self.config.dependencies.check( self.parsed_args )
    
    def generate_argparser (self):

        parser = argparse.ArgumentParser()
        parser.add_argument("mod", nargs="*")
        parser.add_argument('-n','--nucl', metavar = "contig_fasta", action="store",required=True,help = "A FASTA file containing the genome assembly.")
        parser.add_argument('-t','--taxdb', metavar = "taxonomy_db", action="store", required=False, default=None, help = "The location of the taxonomy database, likely provided by an earlier script.")
        parser.add_argument('-s', '--stringency', metavar = "stringency_threshold", required=False, default=0.05, help = "The proportion of annotated non-target points that are allowed to be included in the final selection window. IMPORTANT NOTE: The non-target-annotated points can still be throw-out of the final genome fasta by specifying the --toss_nontarget option.")
        parser.add_argument('-f','--prefix', metavar = 'prefix_for_output', required=False, default='scgid', help="The prefix that you would like to be used for all output files. DEFAULT = scgid")
        parser.add_argument('-g', '--targets', metavar = 'target_taxa', action='store', required=True, help="A comma-separated list with NO spaces of the taxonomic levels that the gc-coverage window should be chosen with respect to including. EXAMPLE: '-g Fungi,Eukaryota,Homo'")
        parser.add_argument('-x', '--exceptions', metavar = 'exceptions_to_target_taxa', action='store', required=False, default=None, help="A comma-separated list with NO spaces of any exlusions to the taxonomic levels specified in -g|--targets. For instance if you included Fungi in targets but want to exclude ascomycetes use: '-x Ascomycota'")

        parser.add_argument('-sp','--augustus_sp', metavar = "augustus_species", action="store",required=False, default=None, help = "Augustus species for gene predicition.")
        parser.add_argument('-e', '--evalue', metavar = 'e-value_cutoff', action = 'store', required = False, default = '1e-10', help = "The evalue cutoff for blast. Default: 1xe-10)")
        parser.add_argument('-db', '--spdb', metavar = 'swissprot_fasta', action='store', required=False, default=None,  help = "The path to your version of the swissprot database in FASTA format.")
        parser.add_argument('--cpus', metavar = 'cores', action = 'store', required = False, default = '1', help = "The number of cores available for blastp to use.")

        # MAYBE PROVIDED BY EARLIER PART OF SCRIPT
        parser.add_argument('-b','--blastout', metavar = "blastout", action="store",required=False, help = "The blast output file from a search of the swissprot database with your proteins as query. Defaults to outfmt=6 and max_target_seqs=1. Provided by earlier script.")
        parser.add_argument('-p','--prot', metavar = "protein_fasta", action="store",required=False, help = "A FASTA file containing the proteins called from the genome.")
        
        return parser
        
    def locate_input_files(self, patterns):
        relpath = os.path.join(
            "{}_scgid_output".format(self.parsed_args.prefix),
            __name__)
        #patterns = [os.path.join(relpath, p) for p in patterns]
        print patterns
        compiled = [re.compile(p) for p in patterns]
        listdir = os.listdir(relpath)
        for c in compiled:
            print [x for x in listdir if re.match(c, x) ]
            
            
        



    def run(self):
        for k,i in self.config.__dict__.iteritems():
            self.log.info("{}: {}".format(k,i))

if __name__ == "__main__":
    pass
