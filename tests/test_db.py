import scgid.db
import inspect
import io
import ast
import pytest

args = {
    "lineages": "/home/aimzez/dev/SCGid/tests/data/new_lineage.tsv",
    "taxdb": "/home/aimzez/dev/SCGid/tests/data/taxdb_in.taxdb"
        }
SPEXPAND_EXPECTED = {
        'Homo sapiens': 'Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Mammalia; Eutheria; Euarchontoglires; Primates; Haplorrhini; Catarrhini; Hominidae; Homo',
        'Pipcy3_1': 'Eukaryota; Fungi; Zoopagomycota; Zoopagomycotina; Piptocephalaceae; Piptocephalis'
        }
config_yaml_before = """
    default_spdb: /path/to/original/spdb
    default_taxdb: /path/to/original/taxdb
    """
config_yaml_after = """
    default_spdb: /path/to/new/spdb
    default_taxdb: /path/to/new/taxdb
    """

def test_expand_taxdb (): 
    taxdb = scgid.db.SPDBTaxonomy(args['taxdb'])
    
    taxdb.expand(args['lineages'])
    
    assert taxdb.taxdb == SPEXPAND_EXPECTED

def test_write_taxdb():
    taxdb = scgid.db.SPDBTaxonomy(args['taxdb'])
    text = None
    with io.StringIO() as buff:
        taxdb.write(buff)
        text = buff.getvalue()

    assert ast.literal_eval(text) == taxdb.taxdb

def test_expand_and_write_taxdb():
    taxdb = scgid.db.SPDBTaxonomy(args['taxdb'])

    taxdb.expand(args['lineages'])

    with io.StringIO() as buff:
        taxdb.write(buff)
        text = buff.getvalue()

    assert ast.literal_eval(text) == SPEXPAND_EXPECTED







