import io
import pytest
from scgid.infotable import InfoTable

infotable_ldict_initial = [{
	"contig": "contig_1",
	"gc": 0.5,
	"coverage": 40.5,
	"pid": "g0001",
	"length": 450,
	"sseqid": "ACC_01",
	"sp_os": "Pipcy3_1",
	"desc": "A description of the protein hit.",
	"evalue": 1e-10
	}, {
	"contig": "contig_2",
	"gc": 0.45,
	"coverage": 45.5,
	"pid": "g0002",
	"length": 760,
	"sseqid": "ACC_02",
	"sp_os": "Pipcy3_1",
	"desc": "A description of the protein hit.",
	"evalue": 1e-12
	}]

infotable_ldict_reload = [{
	"contig": "contig_1",
	"gc": 0.5,
	"coverage": 0.5,
	"pid": "g0001",
	"length": 450,
	"sseqid": "ACC_01",
	"sp_os": "Pipcy3_1",
	"desc": "A description of the protein hit.",
	"evalue": 1e-10,
	"lineage": ["Eukaryota", "Viridiplantae", "Chlorophyta", "Chlorophyceae", "Chlamydomonadales", "Chlamydomonadaceae", "Chlamydomonas"],
	"pertinent_taxlvl": "Viridiplantae",
	"parse_lineage": "nontarget"
	}, {
	"contig": "contig_2",
	"gc": 0.45,
	"coverage": 45.5,
	"pid": "g0002",
	"length": 760,
	"sseqid": "ACC_02",
	"sp_os": "Pipcy3_1",
	"desc": "A description of the protein hit.",
	"evalue": 1e-12,
	"lineage": ["Eukaryota", "Fungi", "Ascomycota", "Pezizomycetes", "Pezizales", "Pezizoaceae", "Peziza", "species"],
	"pertinent_taxlxl": "Fungi",
	"parse_lineage": "target"
	}]

infotable_ldict_reload_multiprotein = [{
	"contig": "contig_1",
	"gc": 0.5,
	"coverage": 0.5,
	"pid": "g0001",
	"length": 450,
	"sseqid": "ACC_01",
	"sp_os": "Pipcy3_1",
	"desc": "A description of the protein hit.",
	"evalue": 1e-10,
	"lineage": ["Eukaryota", "Viridiplantae", "Chlorophyta", "Chlorophyceae", "Chlamydomonadales", "Chlamydomonadaceae", "Chlamydomonas"],
	"pertinent_taxlvl": "Viridiplantae",
	"parse_lineage": "nontarget"
	}, 
	{
	"contig": "contig_1",
	"gc": 0.55,
	"coverage": 0.5,
	"pid": "g0002",
	"length": 500,
	"sseqid": "ACC_02",
	"sp_os": "Pipcy3_1",
	"desc": "A description of the protein hit.",
	"evalue": 1e-99,
	"lineage": ['Eukaryota', 'Fungi', 'Dikarya', 'Ascomycota', 'Saccharomycotina', 'Saccharomycetes', 'Saccharomycetales', 'Saccharomycetaceae', 'Saccharomyces'],
	"pertinent_taxlvl": "Fungi",
	"parse_lineage": "target"
	},
	{
	"contig": "contig_2",
	"gc": 0.45,
	"coverage": 45.5,
	"pid": "g0003",
	"length": 760,
	"sseqid": "ACC_03",
	"sp_os": "Pipcy3_1",
	"desc": "A description of the protein hit.",
	"evalue": 1e-12,
	"lineage": ["Eukaryota", "Fungi", "Ascomycota", "Pezizomycetes", "Pezizales", "Pezizoaceae", "Peziza", "species"],
	"pertinent_taxlxl": "Fungi",
	"parse_lineage": "target"
	},
	{
	"contig": "contig_2",
	"gc": 0.45,
	"coverage": 45.5,
	"pid": "g0003",
	"length": 760,
	"sseqid": "ACC_03",
	"sp_os": "Pipcy3_1",
	"desc": "A description of the protein hit.",
	"evalue": 1e-22,
	"lineage": ["Eukaryota", "Fungi", "Ascomycota", "Pezizomycetes", "Pezizales", "Pezizoaceae", "Peziza", "species"],
	"pertinent_taxlxl": "Fungi",
	"parse_lineage": "target"
	},
	{
	"contig": "contig_2",
	"gc": 0.45,
	"coverage": 45.5,
	"pid": "g0003",
	"length": 760,
	"sseqid": "ACC_03",
	"sp_os": "Pipcy3_1",
	"desc": "A description of the protein hit.",
	"evalue": 1e-22,
	"lineage": ['Bacteria', 'Proteobacteria', 'Alphaproteobacteria', 'Rhizobiales', 'Rhizobiaceae', 'Sinorhizobium/Ensifer group', 'Sinorhizobium'],
	"pertinent_taxlxl": "Proteobacteria",
	"parse_lineage": "nontarget"
	}]


initial_colnames = ["contig","gc","coverage","pid","length","sseqid","sp_os","desc","evalue"]
reload_colnames = ["contig","gc","coverage","pid","length","sseqid","sp_os","desc","evalue","lineage","pertinent_taxlvl","parse_lineage"]

infotable_reload_tsv = io.StringIO(
	'\n'.join([
		'\t'.join([str(x) for x in infotable_ldict_reload[0].values()]),
		'\t'.join([str(x) for x in infotable_ldict_reload[1].values()])
		])
	)
infotable_reload_tsv.seek(0)


'''
def test___init__(self, target_taxa = None, ident="HEAD"):
	pass

def test_spawn_child(self, ident, nested_frame):
	pass

def test_iter_descendants(self, at_head=True):
	pass

def test_target_filter(self):
	pass

def test_collect_unclassifieds(self, nucl):
	pass

def test_tidy (self, taxlvl_idx):
	pass

def test_rfilter(self, feature, window, child_label=""):
	pass

def test_rfilter_inplace_from_parent (self, feature, window):
	pass

def test_reset_from_parent (self):
	pass

def test_tnt_population (self):
	pass
'''

def test_summary_stats ():
	it = InfoTable()
	it.populate(infotable_ldict_reload)


def test_populate():
	it = InfoTable()
	it.populate(ldict = infotable_ldict_initial)
	assert it.df.shape == (2, 9)
	assert it.df.columns.to_list() == initial_colnames
	assert it.df.iloc[0,:].to_list() == list(infotable_ldict_initial[0].values())
	assert it.df.iloc[1,:].to_list() == list(infotable_ldict_initial[1].values())

def test_load():
	it = InfoTable()
	it.load(tsv = infotable_reload_tsv)
	assert it.df.shape == (2, 12)
	assert it.df.columns.to_list() == reload_colnames
	assert it.df.iloc[0,:].to_list() == list(infotable_ldict_reload[0].values())
	assert it.df.iloc[1,:].to_list() == list(infotable_ldict_reload[1].values())

def test_clear_decisions():
	it = InfoTable()
	it.keep = ["list", "of", "contigs", "to", "keep"]
	it.dump = ["list", "of", "contigs", "to", "dump"]

	it.clear_decisions()
	assert it.keep is None
	assert it.dump is None

def test_set_target_no_except ():
	it = InfoTable()
	it.set_target("Fungi, Homo,Viridiplantae")
	assert it.tar == ["Fungi", "Homo", "Viridiplantae"]
	assert it.ex == []

def test_set_target_with_except ():
	it = InfoTable()
	it.set_target("Fungi, Homo,Viridiplantae", "Bacteria, Archaea,Blastocladiomycota")
	assert it.tar == ["Fungi", "Homo", "Viridiplantae"]
	assert it.ex == ["Bacteria", "Archaea", "Blastocladiomycota"]

def test_taxon_level():
	it = InfoTable()
	it.populate(infotable_ldict_reload)
	assert it.taxon_level(2).lineage.to_list() == ["Chlorophyta", "Ascomycota"]
	assert it.taxon_level(3).lineage.to_list() == ["Chlorophyceae", "Pezizomycetes"]

def test_parse_lineage():
	it = InfoTable()
	it.populate(infotable_ldict_reload)
	it.df = it.df.drop("parse_lineage", axis=1)
	it.set_target("Viridiplantae")
	it.parse_lineage()
	assert it.df.parse_lineage.to_list() == ["target", "nontarget"]

def test_decide_inclusion_keep ():
	it = InfoTable()
	it.populate(infotable_ldict_reload_multiprotein)
	it.df = it.df.drop("parse_lineage", axis=1)
	it.set_target("Fungi")
	it.parse_lineage()
	it.decide_inclusion()
	assert it.keep == ["contig_1", "contig_2"]
	assert it.dump == []

def test_decide_inclusion_dump ():
	it = InfoTable()
	it.populate(infotable_ldict_reload_multiprotein)
	it.df = it.df.drop("parse_lineage", axis = 1)
	it.clear_decisions()
	it.set_target("Viridiplantae")
	it.parse_lineage()
	it.decide_inclusion()
	assert it.keep == []
	assert it.dump == ["contig_1", "contig_2"]
