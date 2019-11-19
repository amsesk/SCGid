import logging
import pandas as pd
import argparse
import os
import re
import ast
from scripts.error import Error
import scripts.pkg_settings as pkg_settings
from scripts.modcomm import get_head, logger_name_gen, LoggingEntity
from scripts.library import pkl_fasta_in_out

class PathAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, os.path.abspath(values))

class SPDBTaxonomy(LoggingEntity):
    def __init__ (self, path_to_taxdb):
        self.logger = logging.getLogger ( logger_name_gen() )
        self.taxdb = {}
        with open(path_to_taxdb,'r') as f:
            '''
            for line in f:
                try:
                    ast.literal_eval(line)
                except:
                    self.logger.critical (f"Taxonomy database formatting error. Offending line: {line}")
            '''
            self.taxdb = ast.literal_eval( f.read() )
    
    def __repr__(self):
        return '\n'.join([f"{k}: {v}" for k,v in self.taxdb.items()])
    
    def add_lineage_info (self, parser):
        for p in parser.parsed_hits:
            try:
                lineage = self.taxdb[ p["sp_os"] ]
            except:
                try:
                    # Can get rid of these once local taxdb is rebuilt
                    sp_os = sp_os.replace("'","")
                    sp_os = sp_os.replace("#","")
                    lineage = self.taxdb[sp_os]
                except KeyError:
                    self.logger.warning(f"{sp_os} is missing from the taxdb. You probably specified an older -t|--taxdb file than the one associated with the SPDB blasted against. Rerun build_taxdb or specify a different taxdb to correct. Lineage set as 'Not_in_taxdb' for this run.")
                    lineage = "Not_in_taxdb"
            p["lineage"] = lineage
        return parser

class BlastoutParser(LoggingEntity):
    def __init__(self):
        self.path = None
        self.headers = pkg_settings.BLAST_HEADERS
        self.head = get_head()
        self.logger = logging.getLogger( logger_name_gen() )
        self.hits = []
        self.best_hits = {}
        self.parsed_hits = []

    def load_from_file(self, path_to_blastout):
        self.path = path_to_blastout
        with open( self.path, 'r' ) as blastout:
            for line in blastout:
                self.hits.append([x.strip() for x in line.split("\t")])

    def get_best_hits (self, write = False):
        for hit in self.hits:
            
            bitcol = self.headers["bitscore"]
            query = hit[self.headers["qseqid"]]
            bit = float( hit[self.headers["bitscore"]] )

            if query in self.best_hits.keys():
                if float(self.best_hits[query][bitcol]) < bit:
                    self.best_hits[query] = hit
            else:
                self.best_hits[query] = hit
        self.logger.info("Pulled best blast hits from blast output")
        return 0
    
    def crossref_spdb (self, nucl, prot):
        sp_fasta = pkl_fasta_in_out(self.head.config.get("spdb"), seq_type="prot", contig_info=False)
        ldict = []
        for entry in sp_fasta:
            spl = entry.label.split(" ",1)
            newrow = {
                'accession': spl[0],
                'description': spl[1]
                }
            ldict.append(newrow)
        sp_fasta = None #free
        spdb = pd.DataFrame(ldict).set_index("accession")

        search_pattern = re.compile(pkg_settings.SPDB_OS_REGEXP_PATTERN)
        for query, hit in self.best_hits.items():
            acc = hit[self.headers["sseqid"]]
            evalue = hit[self.headers["evalue"]]
            pid = query[::-1].split(".",1)[0][::-1]
            
            spl = query.split('_')
            contig_shortname = '_'.join( spl[0:2] )

            try:
                desc = spdb.loc[acc].description
                
            except:
                self.logger.critical(MissingAccessionError(accession = acc, db_path = self.head.config.get("spdb")))
            
            s = re.search(search_pattern, desc)
            if s is None:
                self.logger.critical(f"Unable to pull SPDB species information from line: \"{hit}\"")
            sp_os = s.group(1).strip()

            self.parsed_hits.append ({
                "query": contig_shortname, 
                "gc": nucl.get(contig_shortname).gc,
                "coverage": nucl.get(contig_shortname).coverage,
                "pid": pid,
                "length": prot.get(query).length,
                "sseqid": acc,
                "sp_os": sp_os,
                "desc": desc, 
                "evalue": evalue
            })
        
        self.logger.info(f"Cross-referenced blastout with database at {self.head.config.get('spdb')}")
        return 0
    
class MissingAccessionError (Error):
    def __init__(self, accession, db_path):
        super().__init__(db_path)
        self.accession = accession
        self.db_path = db_path
        self.exitcode = 4
        self.fatal = True
    def __str__(self):
        return f"[{self.__name__}] Accession `{self.accession}` missing from database at `{self.db_path}"