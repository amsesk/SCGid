#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 19:58:39 2017

@author: kevinamses
"""

esom_help = """

Usage: scgid kmers <task> [args...] --> try scgid kmers <task> -h|--help to see task-specific arguments

Please specify a task to continue:

    train       Train the ESOM topology
    annotate    Annotate the trained ESOM toplogy with taxonomic information
    extract     Pull-out the target region by class number and generate ESOM-filtered draft
"""

if __name__ == "__main__":
    raise NotImplementedError(f"Calling this script directly is not implemented.")

else:
    import sys
    import os
    import re
    import argparse
    import subprocess
    import numpy as np
    from scripts.module import Module
    from scripts.library import file_grep, subprocessP, subprocessC
    from scripts.modcomm import LoggingEntity, Head, logger_name_gen
    from scripts.parsers import PathAction
    from scripts.dependencies import CaseDependency
    from scripts.reuse import ReusableOutput, nucleotide_blast

    class Kmers(Module, LoggingEntity, Head):
        def __init__(self,  argdict = None):
            super().__init__(self.__class__)

        def update_ESOM_HOME(self, path_to_esomstart, new_home):
            with open(path_to_esomstart, 'r') as f:
                match = 'ESOM_HOME[=]["]'
                reprint = []
                for line in f.readlines():
                    if re.search(match,line) is not None:
                        current_path = re.search('[\"]([^\"]+)[\"]',line).group(1)
                        new_home_line = line.replace(current_path,new_home)
                        reprint.append(new_home_line)
                    else:
                        reprint.append(line)
            with open(path_to_esomstart,'w') as f:
                f.write(''.join(reprint))

        def run(self):

            # Print ESOM module help message if a task has not been selected - malformed argument string
            if len(sys.argv) < 3:
                print (esom_help)
                return 1
            
            elif sys.argv[2].startswith("-"):
                print (esom_help)
                return 1
            
            else:
                if sys.argv[2] == "train":
                    Train().run()
                
                elif sys.argv[2] == "annotate":
                    Annotate().run()

                elif sys.argv[2] == "extract":
                    Extract().run()
                
                else:
                    print("Bad Esom mode selection.")
                    sys.exit(2)


    class Train(Kmers):
        def __init__(self, argdict = None):
            super().__init__(self.__class__)
            if argdict is not None:
                self.config.load_argdict(argdict)
            else:
                self.argparser = self.generate_argparser()
                self.parsed_args = self.argparser.parse_args()
                self.config.load_cmdline( self.parsed_args) # Copy command line args defined by self.argparser to self.config
                
                # Probably better to implement this as a class, but let's wait and see how much it comes up before pouring time into it
                self.case_args = {
                    "mode": {
                        "det": {
                            "bool": self.config.get("cpus") == 1,
                            "warning": "Since `--mode det` cannot run on multiple processors, argument passed to `--cpus` is being ignored. Running on one CPU..."
                        },
                        "somoclu": {
                            "bool": self.config.get("cpus") != 1,
                            "warning": "Since `--mode somoclu` runs MPI, argument `--cpus` should be set. Running on one CPU...` "
                        }
                    }
                }
                self.config.check_case_args(self.case_args)
            
            self.config.dependencies.populate(
                CaseDependency("somoclu", "mode", "s"),
                CaseDependency(self.config.get("mpicmd"), "mode", "s"),
            )

            # Set self.config.esom_path with path to esomtrn in the case that esom_path is not set in config.yaml, although this would be abnormal
            if not os.path.isdir(self.config.get("esom_path")):
                p = subprocess.Popen( ["which", "esomtrn"], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
                path, _ = p.communicate()
                if p.returncode != 0:
                    self.logger.critical("Unable to find ESOM bin directory in PATH or config.yaml. Correct by adding it to PATH or config.yaml.")
                path = os.path.dirname(path).decode("utf-8")
                setattr(self.config, "esom_path", path)
        
        def generate_argparser(self):
            parser = argparse.ArgumentParser()
            parser.add_argument("mod", nargs="*")
            parser.add_argument('-n','--nucl', metavar = "contig_fasta", action=PathAction,required=True, help = "A FASTA file containing the nucleotide assembly. (MANDATORY)")
            parser.add_argument('-f','--prefix', metavar = 'prefix_for_output', required=False, default='scgid', help="The prefix that you would like to be used for all output files. DEFAULT = scgid")

            # print_tetramer_freqs options
            parser.add_argument('-m','--mintig', metavar = "minimum_contig_size", action="store",required=False, default="1000", help = "Contig size cutoff (in nucleotides) for inclusion in the ESOM training. Default = 1000 bp")
            parser.add_argument('-w','--window', metavar = "window_size", action="store",required=False, default="1000", help = "Size of the window in which kmer frequencies are calculated. Default = 1000 bp")
            parser.add_argument('-k','--kmer', metavar = "kmer_size", action="store",required=False, default="4", help = "Kmer size for which frequencies will be calculated. Default = 4")

            # esomtrn options
            parser.add_argument('--mode', metavar = "mode", action="store", required=False, choices=["det","somoclu"], default = "det", help = "Mode to train the ESOM. [det|s]")
            parser.add_argument('--cpus', metavar = "mode", action="store", required=False, default=1, help = "Number of CPUs to use for training (Somoclu only)")
            parser.add_argument('-r','--rows', metavar = "rows_in_map", action="store",required=False, help = "The number of rows to be present in the output ESOM. Default = 5.5x the number of neurons")
            parser.add_argument('-c','--cols', metavar = "columns_in_map", action="store",required=False, help = "The number of columns to be present in the output ESOM. Default = 5.5x the number of neurons")
            parser.add_argument('-sr','--start_radius', metavar = "start_radius", action="store",required=False, default = '50', help = "Start radius for the ESOM. (MANDATORY)")
            parser.add_argument('-e','--epochs', metavar = "training_epochs", action="store",required=False, default = "20", help = "Number of epochs to train over. Default = 20")
            parser.add_argument('--Xmx', metavar = "available_memory", action="store",required=False, default = "512m", help = "Set memoray available to train the ESOM. Specicy as such: X megabytes = Xm, X gigabytes = Xg")

            return parser

        def specify_Xmx_esom(self, path_to_esomstart, mem):

            with open(path_to_esomstart,'r') as f:
                pos = 0
                trunc_point = 0
                for line in f.readlines():
                    if line.split(' ')[0] == 'java':
                        trunc_point = pos
                        break
                    pos += len(line)

            with open(path_to_esomstart,'a') as f:
                f.seek(trunc_point)
                f.truncate()
                new_line = "java -cp \"$CP\" -Xmx"+mem+" $MAIN \"$@\""
                f.write(new_line)

        def call_print_tetramer_freqs(self):
            self.logger.info("Calculating 4-mer frequencies across scaffolds and generating ESOM input files.")
            
            nucl_path = self.config.get("nucl")
            
            # Copy assembly to ESOM working directory since perl script will place outputs whereever it is, and we want them here
            cmd = ['cp', nucl_path, os.path.basename(nucl_path)]
            subprocessC(cmd)

            # Call `print_tetramer_freqs_deg_filter_esom_VD.pl`
            script_path = os.path.join(self.config.SCGID_SCRIPTS, 'print_tetramer_freqs_deg_filter_esom_VD.pl')
            cmd = ['perl', script_path,
                '-s', os.path.basename(nucl_path),
                '-m', self.config.get("mintig"),
                '-w', self.config.get("window"),
                '-k', self.config.get("kmer")
                ]

            self.logger.info(' '.join(cmd))

            subprocessP(cmd, self.logger)

        def determine_map_size(self):
            if self.config.get("rows") is None or self.config.get("cols") is None:
                nodes = len(open(f"{os.path.basename(self.config.get('nucl'))}.lrn").readlines()) - 4 # number of lines in .lrn file minus 4 header lines
                neurons = float(nodes*5.5)
                dim = float(np.sqrt(neurons))
                setattr( self.config, "rows", str(int(np.ceil(dim))) )
                setattr( self.config, "cols", str(int(np.ceil(dim))) )
                self.logger.info(f"Nothing specified for -r|--rows or -c|--cols. There are {nodes:,} nodes in the lrn file, using defult dimensions: {self.config.get('rows')}x{self.config.get('cols')}" )
            else:
                self.logger.info("Using user-specified dimensions for map: %sx%s",self.config.get("rows"), self.config.get("cols"))

        def train_det (self):
            rows = self.config.get("rows")
            cols = self.config.get("cols")
            start_radius = self.config.get("start_radius")
            epochs = self.config.get("epochs")
            
            nucl_basename = os.path.basename(self.config.get("nucl"))
            wts = f"{nucl_basename}.{rows}x{cols}e{epochs}.wts"
            bm = f"{nucl_basename}.{rows}x{cols}e{epochs}.bm"
            clsfile = f"{nucl_basename}.cls"
            lrn = f"{nucl_basename}.lrn"



            #Alter hard-coded memory-availability in esomstart to match user-specified resources
            self.specify_Xmx_esom(
                os.path.join( self.config.get("esom_path"), "esomstart"), self.config.get("Xmx")
                )

            self.logger.info(f"Training ESOM in mode `{self.config.get('mode')}`")
            cmd = [
                os.path.join( self.config.get("esom_path"), "esomtrn"),
                "--permute",
                "--out", wts,
                "-b", bm,
                "--cls", clsfile,
                "--lrn", lrn,
                "--algorithm", "kbatch",
                "--rows", rows,
                "--columns", cols,
                "-bmc", "6",
                "--start-radius", start_radius,
                "--epochs", epochs,
                "-k", "0.15",
                "--bmsearch", "standard",
                "--dist", "euc"
                ]

            self.logger.info(' '.join(cmd))

            subprocessP(cmd, [self.logger, self.simplelogger], log_stdout=True)

            self.logger.info("ESOM train call to Databionic ESOM Tools returned.")
            self.logger.warning("SCGid is unable to catch Java OutOutMemory (OOM) handlers. Check logfile to verify successful training completion.")

            if not os.path.isfile(wts) or not os.path.isfile(bm):
                self.logger.critical("Unable to locate *.wts or *.bm files associated with successful training. Check logfile for error output - Java OutOfMemory error the likely culprit.")
            
            return 0
        
        def train_somoclu(self):
            cpus = self.config.get("cpus")
            rows = self.config.get("rows")
            cols = self.config.get("cols")
            start_radius = self.config.get("start_radius")
            epochs = self.config.get("epochs")

            nucl_basename = os.path.basename(self.config.get("nucl"))
            lrn = f"{nucl_basename}.lrn"

            self.logger.info(f"Training ESOM in mode `{self.config.get('mode')}`")

            cmd = [
                self.config.get("mpicmd"),
                "-np", cpus,
                "somoclu",
                "-e", epochs,
                "-l", "0.5",
                "-L", "0.1",
                "-m", "toroid",
                "-r", start_radius,
                "-x", rows,
                "-y", cols,
                "-v", "2",
                lrn,
                nucl_basename
                ]

            self.logger.info(" ".join(cmd))

            subprocessP(cmd, [self.logger, self.simplelogger], log_stdout=True)

            self.logger.info("ESOM train call to somoclu returned.")

            return 0

        def run(self):
            self.start_logging()
            self.setwd( __name__, self.config.get("prefix") )
            self.config.reusable.check()
            self.config.dependencies.check(self.config)
            self.config.reusable.generate_outputs()

            self.call_print_tetramer_freqs()
            self.determine_map_size()

            ##Train the ESOM
            if self.config.get("mode") == "det":

                self.train_det()

            elif self.config.get("mode") == "somoclu":
                
                self.train_somoclu()
            
            # Literally cannot happen because of argparser choices, but just for the principle of the matter
            else:
    
                pass
    
    class Annotate(Kmers):
        def __init__(self, argdict = None):
            super().__init__(self.__class__)
            if argdict is not None:
                self.config.load_argdict(argdict)
            else:
                self.argparser = self.generate_argparser()
                self.parsed_args = self.argparser.parse_args()
                self.config.load_cmdline( self.parsed_args) # Copy command line args defined by self.argparser to self.config

            self.config.reusable.populate(
                ReusableOutput(
                    arg = "blastout",
                    pattern = ".*[.]nt[.]blast[.]out$",
                    genfunc = nucleotide_blast,
                    genfunc_args = {
                        "nucl": self.config.get("nucl"),
                        "db": "nt",
                        "evalue": self.config.get("evalue"),
                        "cpus": self.config.get("cpus"),
                        "outpath": f"{self.config.get('prefix')}.nt.blast.out"
                    }
                )
            )

            self.config.dependencies.populate(
                CaseDependency("blastn", "blastout", None)
            )
        
        def generate_argparser(self):
            parser = argparse.ArgumentParser()

            parser.add_argument("mod", nargs="*")
            parser.add_argument('-n','--nucl', metavar = "contig_fasta", action=PathAction,required=True, help = "A FASTA file containing the nucleotide assembly. (MANDATORY)")
            parser.add_argument('-f','--prefix', metavar = 'prefix_for_output', required=False, default='scgid', help="The prefix that you would like to be used for all output files. DEFAULT = scgid")
            parser.add_argument('-b','--blastout', metavar = "blastout", action=PathAction,required=False, help = "The blast output file from a blastn search of the NCBI nt database with your contigs as query. If you have not done this yet, this script will do it for you.")

            parser.add_argument('-m','--mintig', metavar = "minimum_contig_size", action="store",required=False, default="1000", help = "Contig size cutoff (in nucleotides) for inclusion in the ESOM training. Default = 1000 bp")
            parser.add_argument('-w','--window', metavar = "window_size", action="store",required=False, default="1000", help = "Size of the window in which kmer frequencies are calculated. Default = 1000 bp")

            parser.add_argument('--cpus', metavar = 'cores', action = 'store', required = False, default = "1", help = "The number of cores available for BLAST to use.")
            parser.add_argument('-e', '--evalue', metavar = 'e-value_cutoff', action = 'store', required = False, default = '1e-5', help = "The evalue cutoff for blast. Default: 1xe-5)")
            #parser.add_argument('-rm', '--rankmode', action = 'store_true', required = False, help = "Annotate contigs at the same single taxonomic rank across all contigs. (e.g. superkingdom)")
            #parser.add_argument('-te', '--targetexcept', action = 'store_true', required = False, help = "Annotate contigs at a varietry of taxonomic ranks across all contigs. (e.g. Eukaryota, except Fungi). Must be used in combination with -s|--annotation_scheme.")
            parser.add_argument('-s','--annotation_scheme', metavar = 'annotation_scheme', action = 'store', required = True, help = "The annotation scheme to use in target_except annotation mode (RECOMMENDED). Groups MUST be mutually exclusive to avoid overlap and unincluded groups will be arbitrarily marked as 'Unclassified'. You ABSOLUTELY MUST use this basic syntax: '-s|--annotation_schem target1^exception1,excetions2/target2^exception1/target3/etc...'. Example: '-s|--annotation_scheme Eukaryota^Fungi,Metazoa/Fungi/Metazoa/Bacteria^Proteobacteria/Proteobacteria'. See documention for detailed information and more examples.")
            #parser.add_argument('-i','--infotable', metavar = "infotable", action="store",required=False, help = "The scgid gc-cov-derived infotable generated by a blastp search of a swissprot-style protein database.")
            #parser.add_argument('--mode', metavar = "mode", action="store",required=False, default ='blastn', help = "The type of blast results that you would like to use to topology ('blastp' or 'blastn'). This module will automatically do a blastn search of the NCBI nt database for you. At this time, a blastp search can not be run directly from this script. INSTEAD, if using mode 'blastp' (DEFAULT, recommended) you must specify a scgid blob-derived <prefix>_info_table.tsv file with -i|--infotable")

            return parser

        def run(self):
            self.start_logging()
            self.setwd( __name__, self.config.get("prefix") )
            self.config.reusable.check()
            self.config.dependencies.check(self.config)
            self.config.reusable.generate_outputs()

    class Extract(Kmers):
        def __init__(self):
            pass
        
        def run(self):
            pass

    '''
    ## Make sure that esom path is set correctly in the esomstart and esomtrn ##
    esomstart = os.path.join(settings.esom_path,"bin","esomstart")
    esomstart_declare = file_grep("ESOM_HOME=",esomstart)
    current_path = re.search('[\"]([^\"]+)[\"]',esomstart_declare).group(1)
    if current_path is not settings.esom_path:
        update_ESOM_HOME (esomstart,settings.esom_path)

    esomtrn = os.path.join(settings.esom_path,"bin","esomtrn")
    esomtrn_declare = file_grep("ESOM_HOME=",esomtrn)
    current_path = re.search('[\"]([^\"]+)[\"]',esomtrn_declare).group(1)
    if current_path is not settings.esom_path:
        update_ESOM_HOME (esomtrn,settings.esom_path)


    if sys.argv[1] == "train":
        arguments = sys.argv[2:]
        call = os.path.join(pkg_home,'bin','train_esom.py')
        arguments.insert(0,call)
        py = sys.executable
        arguments.insert(0,py)
        subprocess.call(arguments)

    elif sys.argv[1] == "annotate":
        arguments = sys.argv[2:]
        call = os.path.join(pkg_home,'bin','annotate_esom.py')
        arguments.insert(0,call)
        py = sys.executable
        arguments.insert(0,py)
        subprocess.call(arguments)

    elif sys.argv[1] == "extract":
        arguments = sys.argv[2:]
        call = os.path.join(pkg_home,'bin','extract_esom.py')
        arguments.insert(0,call)
        py = sys.executable
        arguments.insert(0,py)
        subprocess.call(arguments)
    '''
