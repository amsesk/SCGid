import logging
import pandas as pd
import argparse
import os
import re
import warnings
from collections import namedtuple, OrderedDict
from ete3 import NCBITaxa
from scgid.error import ModuleError, ErrorClassNotImplemented, Ok, check_result
import scgid.pkg_settings as pkg_settings
from scgid.modcomm import get_head, logger_name_gen, LoggingEntity, ErrorHandler
from scgid.sequence import AASequenceCollection
from scgid.db import SPDBTaxonomy

class MalformedDatabaseHeaderError(ModuleError):
    def __init__(self, offender):
        super().__init__()
        self.msg = f"Unable to pull SPDB species information from line: \"{offender}\"."
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
            #ModuleError(f"Unable to access input file: {absolute_path}")
            raise IOError(f"Error parsing arugment `{'|'.join(self.option_strings)}`. Input file does not exist: {absolute_path}")
        setattr(namespace, self.dest, absolute_path)

class StoreInt(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        try:
            setattr(namespace, self.dest, int(values))

        except ValueError:
            raise ValueError(f"Value `{values}` given for {'|'.join(self.option_strings)} cannot be coerced to int.")

        except:
            return ErrorClassNotImplemented()

class OutputPathStore(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, os.path.abspath(values))

class BlastoutParser(LoggingEntity, ErrorHandler):
    def __init__(self):
        self.path = None
        self.headers = pkg_settings.BLAST_HEADERS
        self.logger = logging.getLogger( logger_name_gen() )
        self.hits = []
        self.best_hits = {}
        self.parsed_hits = []
        self.spdb_tax = None

        try:
            self.head = get_head()

        except:
            self.head = None

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

            if len(spl) != 2:
                return MalformedDatabaseHeaderError(offender = entry.header)

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
                return MalformedDatabaseHeaderError(offender = hit)

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

    @check_result
    def add_lineages (self, path_to_taxdb):

        self.spdb_tax = SPDBTaxonomy(path_to_taxdb)

        for p in self.parsed_hits:
            sp_os = p["sp_os"]
            try:
                lineage = self.spdb_tax.taxdb[ sp_os ]
            except:
                try:
                    # Can get rid of these once local taxdb is rebuilt
                    sp_os = sp_os.replace("'","")
                    sp_os = sp_os.replace("#","")
                    lineage = self.spdb_tax.taxdb[sp_os]
                except KeyError:
                    self.logger.warning(f"{sp_os} is missing from the taxdb. You probably specified an older -t|--taxdb file than the one associated with the SPDB blasted against. Rerun build_taxdb or specify a different taxdb to correct. Lineage set as 'Not_in_taxdb' for this run.")
                    lineage = "Not_in_taxdb"
            p["lineage"] = lineage
        return Ok(None)

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
