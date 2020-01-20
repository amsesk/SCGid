import logging
import pandas as pd
import argparse
import os
import re
import ast
import warnings
from collections import namedtuple
from ete3 import NCBITaxa
from scgid.error import ModuleError, ErrorClassNotImplemented
import scgid.pkg_settings as pkg_settings
from scgid.modcomm import get_head, logger_name_gen, LoggingEntity, ErrorHandler
from scgid.sequence import AASequenceCollection

class MalformedTaxonomyDatabaseError(ModuleError):
    def __init__(self, taxdbpath):
        super().__init__()
        self.msg = f"Taxonomy database at `{taxdbpath}` improperly formatted. Run `build_taxdb.py` to generate a new one."
        self.catch()

class MalformedLineageFileError(ModuleError):
    def __init__(self):
        super().__init__()
        self.msg = "Malformed lineage file supplied as argument to `-l|--lineages`."
        self.catch()

class MalformedDatabaseHeaderError(ModuleError):
    def __init__(self, hit):
        super().__init__()
        self.msg = f"Unable to pull SPDB species information from line: \"{hit}\"."
        self.catch()

class MissingAccessionError(ModuleError):
    def __init__(self, accession, spdbpath):
        super().__init__()
        self.msg = f"Accession `{accession}` is missing from the SPDB at `{spdbpath}`. This occurs when a different SPDB is specified from that which BLASTed against."
        self.catch()

class PathAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):

        # Deprecating for PathStore usage instead
        warnings.warn("PathAction deprecated, use PathStore instead", DeprecationWarning)

        setattr(namespace, self.dest, os.path.abspath(values))

# Check validity of Path before storing. Raises IO error if path is invalid, stores in argparser if valid.
class PathStore(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        absolute_path = os.path.abspath(values)
        if not os.path.exists(absolute_path):
            print(self.option_strings)
            #ModuleError(f"Unable to access input file: {absolute_path}")
            raise IOError(f"Error parsing arugment `{'|'.join(self.option_strings)}`. Input file does not exist: {absolute_path}")
        setattr(namespace, self.dest, absolute_path)

class OutputPathStore(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, os.path.abspath(values))

class SPDBTaxonomy(LoggingEntity, ErrorHandler):
    def __init__ (self, path_to_taxdb):
        self.logger = logging.getLogger ( logger_name_gen() )
        self.taxdb = {}
        with open(path_to_taxdb,'r') as f:
            try:
                self.taxdb = ast.literal_eval( f.read() )

            except SyntaxError:
                return MalformedTaxonomyDatabaseError(taxdbpath = self.head.config.get("taxdb"))

            except:
                return ErrorClassNotImplemented()

    def __repr__(self):
        return '\n'.join([f"{k}: {v}" for k,v in self.taxdb.items()])

    # This is a bad name since it specifically applies to InfoTable
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

    def expand(self, path_to_lineage) -> None:
        with open(path_to_lineage) as f_in:
            lines = f_in.readlines()
        
        # Split file into key,lineage for taxdb (hopefully)
        spl = [ [i.strip() for i in l.split("\t")] for l in lines]

        # Check length of splits to ensure that there are TWO columns, one key with one tab-separated lineage
        right_ncols = [len(s) == 2 for s in spl]
        
        if not all(right_ncols):
            return MalformedLineageFileError()
        
        else:
            lines_parsed = {i[0]: i[1] for i in spl} 
            self.taxdb.update(lines_parsed)

        return None

    # Takes open file object - use within with statment
    # This should be done with most functions that read/write to files I think - makes tests with StringIO easy
    def write(self, f_out) -> None:
        f_out.write("{\n")
        f_out.write(",\n".join([f"'{key}': '{lin}'" for key,lin in self.taxdb.items()]))
        f_out.write("\n}")
        return None



class BlastoutParser(LoggingEntity, ErrorHandler):
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
        sp_fasta = AASequenceCollection().from_fasta(self.head.config.get("spdb"))
        ldict = []
        for entry in sp_fasta.seqs():
            spl = entry.header.split(" ",1)
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

            except KeyError:
                return MissingAccessionError( accession=acc, spdbpath=self.head.config.get('spdb') )

            except:
                return ErrorClassNotImplemented()

            s = re.search(search_pattern, desc)
            
            if s is None:
                return MalformedDatabaseHeaderError(hit = hit)

            sp_os = s.group(1).strip()

            self.parsed_hits.append ({
                "contig": contig_shortname,
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

    def taxids(self):
        ret = []
        TaxId = namedtuple("TaxId", ["query", "taxids"])
        for q,h in self.best_hits.items():
            ret.append( TaxId(q, h[pkg_settings.BLAST_HEADERS["staxids"]] ) )
        return ret

    def ncbi_taxrpt (self, taxids):
            self.logger.info("Using taxids from blastout to pull lineage information from ncbi taxonomy database...")
            ncbi = NCBITaxa()
            ids = {}
            ldict = []
            for t in taxids:
                spl = t.taxids.split(';')
                if len(spl) > 1:
                    ids[t.query] = {}
                    for tid in spl:
                        num = int(tid)
                        if t.query in ids.keys():
                            ids[t.query].update(ncbi.get_taxid_translator(ncbi.get_lineage(num)))
                        else:
                            ids[t.query] = ncbi.get_taxid_translator(ncbi.get_lineage(num))
                else:
                    num = int(t.taxids)
                    ids[t.query] = ncbi.get_taxid_translator(ncbi.get_lineage(num))
                vals  = ids[t.query].values()
                ranks = ncbi.get_rank(ids[t.query]).values()
                row = {k:v for k,v in list(zip(ranks,vals))}
                row["query"] = t.query
                ldict.append( row )
            contig_lineages = pd.DataFrame(ldict)[["query","superkingdom","kingdom","phylum","class","order","family","genus","species"]]
            contig_lineages.set_index("query")
            return contig_lineages
