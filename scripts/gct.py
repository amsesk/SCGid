import logging
import inspect
import sys
import re
import os
from Depends import CaseDependency, ConstDependency
from Module import Module, Config
import argparse

class Gct(Module):
    def __init__(self):
        super(Gct, self).__init__(self.__class__)
        frm = inspect.stack()[1]
        self.parent = inspect.getmodule(frm[0])
        self.argparser = self.generate_argparser()
        self.parsed_args = self.argparser.parse_args()
        self.config.load_cmdline( self.parsed_args )
        self.start_logging()
        setattr(self.config, "files_to_generate", self.locate_input_files(
            {
            "prot": ".*[.]aug[.]out[.]fasta$",
            "blastout": ".*[.]spdb[.]blast[.]out$"
            }
        ))

        self.config.dependencies.populate(
            CaseDependency("blastp", "blastout", None),
            CaseDependency("augustus", "prot", None),
        )

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
    
    def locate_slow_step_outputs(self, slowsteps):
        for i,step in enumerate(slowsteps):
            if self.config.get(step.arg) is not None:
                slowsteps.pop(i)
        to_generate = {}
        update_config = {}
        relpath = os.path.join(
            "{}_scgid_output".format(self.parsed_args.prefix),
            __name__)
        compiled = { arg: re.compile(p) for (arg,p) in patterns.iteritems() }
        listdir = os.listdir(relpath)
        for arg, c in compiled.iteritems():
            matches = [os.path.join(relpath, fname) for fname in listdir if re.match(c,fname)]
            if len(matches) == 0:
                self.log.info("No match found for missing argument `%s`, planning to run additional step to generate", arg)
                to_generate[arg] = True
            elif len(matches) == 1:
                updated_arg = os.path.abspath(matches[0])
                update_config[arg] = { arg: updated_arg }
                to_generate[arg] = False
                self.log.info("Found matching file for missing argument `%s` at `%s`", arg, updated_arg)
            else:
                self.log.critical("Found multiple files matching pattern for argument `%s`, doing nothing")

    def local_input_files(self, patterns):
        for arg in patterns.keys():
            if self.config.get(arg) is not None:
                patterns.pop(arg)
        to_generate = {}
        update_config = {}
        relpath = os.path.join(
            "{}_scgid_output".format(self.parsed_args.prefix),
            __name__)
        compiled = { arg: re.compile(p) for (arg,p) in patterns.iteritems() }
        listdir = os.listdir(relpath)
        for arg, c in compiled.iteritems():
            matches = [os.path.join(relpath, fname) for fname in listdir if re.match(c,fname)]
            if len(matches) == 0:
                self.log.info("No match found for missing argument `%s`, planning to run additional step to generate", arg)
                to_generate[arg] = True
            elif len(matches) == 1:
                updated_arg = os.path.abspath(matches[0])
                update_config[arg] = { arg: updated_arg }
                to_generate[arg] = False
                self.log.info("Found matching file for missing argument `%s` at `%s`", arg, updated_arg)
            else:
                self.log.critical("Found multiple files matching pattern for argument `%s`, doing nothing")

        for arg, fname in update_config.iteritems():
            setattr(self.config, arg, fname)
        return to_generate

    def run(self):
        self.setwd( __name__, self.config.get("prefix") )
        

if __name__ == "__main__":
    pass
