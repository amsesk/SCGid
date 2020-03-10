if __name__ == "__main__":
    raise NotImplementedError(f"Calling this script directly is not implemented.")

else:
    import logging
    import sys
    import re
    import os
    import argparse
    import shutil
    import io
    from scgid.dependencies import CaseDependency, ConstDependency
    from scgid.module import Module
    from scgid.reuse import ReusableOutputManager, ReusableOutput, augustus_predict, protein_blast
    from scgid.modcomm import LoggingEntity, Head
    from scgid.parsers import BlastoutParser, PathStore, SPDBTaxonomy
    from scgid.infotable import InfoTable, it_get_taxonomy_level
    from scgid.sequence import DNASequenceCollection, AASequenceCollection
    from scgid.flexwindow import generate_windows
    from scgid.error import Ok
    from scgid.plotter import PlotlyPlotter

    class Gct (Module, LoggingEntity, Head):
        def __init__(self, argdict = None, loglevel=logging.INFO):
            super().__init__(self.__class__, loglevel=loglevel)
            if argdict is not None:
                self.translated_args = self.translate_argdict(argdict, Gct.generate_argparser())
                self.config.load_argdict(self.translated_args)
                self.parsed_args = self.config
            else:
                self.argparser = Gct.generate_argparser()
                self.parsed_args = self.argparser.parse_args()
                self.config.load_cmdline( self.parsed_args ) # Copy command line args defined by self.argparser to self.config

            self.set_rundir(self.config.get("prefix"))

            self.log_to_rundir(type(self).__name__)
            
            self.config.reusable.populate(
                ReusableOutput(
                    arg = "prot",
                    pattern = ".*[.]aug[.]out[.]fasta$",
                    genfunc = augustus_predict,
                    genfunc_args = {
                        "prefix": self.config.get("prefix"),
                        "nucl": self.config.get("nucl"),
                        "augustus_sp": self.config.get("augustus_sp"),
                        "outpath": f"{self.config.get('prefix')}.aug.out.fasta"
                    }
                ),
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
                    })
            )
            self.config.dependencies.populate(
                CaseDependency("blastp", "blastout", None),
                CaseDependency("augustus", "prot", None),
            )

            self.infotable = InfoTable()
            self.unclassified_infotable = InfoTable()

            self.keep = {}
            self.dump = {}

        def generate_argparser ():

            parser = argparse.ArgumentParser()
            parser.add_argument("mod", nargs="*")
            parser.add_argument('-n','--nucl', metavar = "assembly_fasta", action=PathStore,required=True,help = "A FASTA file containing the genome assembly.")
            parser.add_argument('-s', '--stringency', metavar = "stringency_threshold", required=False, default=0.05, help = "The proportion of annotated non-target points that are allowed to be included in the final selection window. IMPORTANT NOTE: The non-target-annotated points can still be throw-out of the final genome fasta by specifying the --toss_nontarget option.")
            parser.add_argument('-f','--prefix', metavar = 'output_prefix', required=False, default='scgid', help="The prefix that you would like to be used for all output files. DEFAULT = scgid")
            parser.add_argument('-g', '--targets', metavar = 'target_taxa', action='store', required=True, help="A comma-separated list with NO spaces of the taxonomic levels that the gc-coverage window should be chosen with respect to including. EXAMPLE: '-g Fungi,Eukaryota,Homo'")
            parser.add_argument('-x', '--exceptions', metavar = 'exception_taxa', action='store', required=False, default=None, help="A comma-separated list with NO spaces of any exlusions to the taxonomic levels specified in -g|--targets. For instance if you included Fungi in targets but want to exclude ascomycetes use: '-x Ascomycota'")

            parser.add_argument('-sp','--augustus_sp', metavar = "augustus_species", action="store",required=False, default=None, help = "Augustus species for gene predicition.")
            parser.add_argument('-e', '--evalue', metavar = 'blast_evalue_cutoff', action = 'store', required = False, default = '1e-10', help = "The evalue cutoff for blast. Default: 1xe-10)")
            parser.add_argument('-db', '--spdb', metavar = 'swissprot_fasta', action=PathStore, required=False, default=None,  help = "The path to your version of the swissprot database in FASTA format.")
            parser.add_argument('-t','--taxdb', metavar = "swissprot_taxdb", action=PathStore, required=False, default=None, help = "The location of the taxonomy database, likely provided by an earlier script.")
            parser.add_argument('--cpus', metavar = 'cpus', action = 'store', required = False, default = '1', help = "The number of cores available for blastp to use.")

            # MAYBE PROVIDED BY EARLIER PART OF SCRIPT
            parser.add_argument('-b','--blastout', metavar = "spdb_blast_output", action=PathStore,required=False, help = "The blast output file from a search of the swissprot database with your proteins as query. Defaults to outfmt=6 and max_target_seqs=1. Provided by earlier script.")
            parser.add_argument('-p','--prot', metavar = "protein_fasta", action=PathStore, required=False, help = "A FASTA file containing the proteins called from the genome.")

            return parser

        def run(self):
            self.setwd( __name__ )

            self.print_opts()

            self.config.reusable.check()
            self.config.dependencies.check(self.config)
            self.config.reusable.generate_outputs()

            nucl = DNASequenceCollection().from_fasta(self.config.get("nucl"), spades = True)
            nucl.rekey_by_shortname()
            self.logger.info(f"Read nucleotide fasta at `{self.config.get('nucl')}`")

            prot = AASequenceCollection().from_fasta(self.config.get("prot"))
            self.logger.info(f"Read protein fasta at `{self.config.get('prot')}`")

            p = BlastoutParser()
            p.load_from_file(self.config.get("blastout"))
            p.get_best_hits()
            p.crossref_spdb(nucl, prot) # This function needs to be split-up and its functionality spread out
            spdb_tax = SPDBTaxonomy(self.config.get("taxdb"))
            p = spdb_tax.add_lineage_info(p)


            # Fill and do some reformatting of infotable - then parse the lineage information down to target/nontarget
            taxlvl_idx = 1
            self.infotable.set_target(
                self.config.get("targets"),
                self.config.get("exceptions")
                )
            colnames = p.parsed_hits[0].keys()
            self.infotable.populate(p.parsed_hits, colnames)

            # Various clean-up actions (e.g., remove problem characters and split lineage into list)
            self.infotable.tidy(taxlvl_idx)

            #Parse lineage info into target|nontarget|unclassified
            self.infotable.parse_lineage()

            # Create infotable object for unclassified contigs as well (i.e., contigs with no protein hits)
            self.unclassified_infotable = self.infotable.collect_unclassifieds(nucl)

            # Write infotables for classified and unclassified contigs to csv
            self.infotable.df.drop("sseqid", axis=1).to_csv(f"{self.config.get('prefix')}.infotable.tsv", sep='\t', index = False, header = False)
            self.unclassified_infotable.df.to_csv(f"{self.config.get('prefix')}.unclassified.infotable.tsv", sep = '\t', index = False, header = False)

            # Generate all 13 windows
            windows = generate_windows(self.infotable)
            
            # Print PDFs of windows to pdf in directory `windows` and stats on all windows to tsv
            '''
            if os.path.isdir('../windows'):
                shutil.rmtree('../windows')
            os.mkdir('windows')
            #windows.print_all_pdf("windows")
            '''
            plotout = f"{self.config.get('prefix')}.gctplt.html"
            plotter = PlotlyPlotter(infotable = self.infotable, n = 10)
            plotter.plot(outpath = plotout)

            windows.print_all_tsv(f"{ self.config.get('prefix') }.windows.all.out")

            # Pick the best window and store in `best`
            best_window = windows.pick( self.config.get("stringency") )
            self.logger.info( f"Best window at stringency level `s = {self.config.get('stringency')}`:\n\n{best_window.show()} ")

            # Decide which CLASSIFIED contigs to keep based on taxonomy of their proteins
            # Populates infotable.keep and infotable.dump
            self.infotable.decide_inclusion()

            # Pop DNASequences from nucl based on infotable decisions
            # MUST POP because of subsequent filtering by GC-Coverage
            for shortname in self.infotable.keep:
                self.keep[shortname] = nucl.pop(shortname)
            for shortname in self.infotable.dump:
                self.dump[shortname] = nucl.pop(shortname)

            # Decide which UNCLASSIFIED contigs to keep based on GC-Coverage cut-offs defined by the best window
            sort = nucl.gc_cov_filter(best_window.gc_range, best_window.coverage_range)
            self.keep.update(sort["keep"])
            self.dump.update(sort["dump"])

            # Construct DNASequenceCollection from final filtered assembly,
            # Resort so contigs are in original order
            final_assembly = DNASequenceCollection().from_dict(
                { k: self.keep[k] for k in sorted(self.keep) }
            )

            # Compute final filtered assembly stats
            filtered_size = sum([len(s.string) for s in final_assembly.seqs()])
            filtered_ncontigs = len(final_assembly.seqs())

            self.logger.info(f"Filtered assembly contains {filtered_ncontigs:,} contigs with a cumulative size of {filtered_size:,} bp ({filtered_size/1e6:.2f} Mbp).")
            
            # Print final filtered assembly to FASTA
            final_fname = f"{self.config.get('prefix')}.gct.filtered.assembly.fasta"
            final_assembly.write_fasta( final_fname )

            self.logger.info(f"Final filtered assembly written in FASTA format to `{final_fname}`")

            self.logger.info("GCT-based filtering complete. Returning to SCGid.")

            # Migrate and then remove temp dir, cd back to starting dir
            self.migrate_temp_dir()
            self.resetwd()

            # Return final filtered assembly to SCGid root
            return (Ok(), final_assembly)
