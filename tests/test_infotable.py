import io
import pytest
from scgid.infotable import InfoTable

ldict = [{
	"contig": "contig_1",
	"gc": 0.5,
	"coverage": 40.5,
	"pid": "g0001",
	"length": 450,
	"sseqid": "ACC_01",
	"sp_os": "OS=Pipcy3_1",
	"desc": "A description of the protein hit.",
	"evalue": 1e-10
	}, {
	"contig": "contig_2",
	"gc": 0.45,
	"coverage": 45.5,
	"pid": "g0002",
	"length": 760,
	"sseqid": "ACC_02",
	"sp_os": "OS=Pipcy3_1",
	"desc": "A description of the protein hit.",
	"evalue": 1e-12
	}]
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

def test_summary_stats (self, feature):
	pass
'''
def test_populate():
	it = InfoTable()
	it.populate(ldict)
	assert it.df.shape == (2, 9)
	assert it.df.iloc[0,:].to_list() == ["contig_1", 0.5, 40.5, "g0001", 450, "ACC_01", "OS=Pipcy3_1", "A description of the protein hit.", 1e-10]
	assert it.df.iloc[1,:].to_list() == ["contig_2", 0.45, 45.5, "g0002", 760, "ACC_02", "OS=Pipcy3_1", "A description of the protein hit.", 1e-12]
'''
def test_load(self, tsv):
	pass

def test_clear_decisions(self):
	pass

def test_set_target(self, target, exceptions = None):
	pass

def test_taxon_level(self, level):
	pass

def test_parse_lineage(self):
	pass

def test_decide_inclusion (self):
	pass
	'''