import io
import pytest
from scgid.infotable import InfoTable
from scgid.sequence import DNASequenceCollection, DNASequence
import pandas as pd

nucl_fasta = {
"contig_1": DNASequence("contig_1", "ATCGGGACTGGTAGATA"),
"contig_2": DNASequence("contig_2", "CCGGGCCTAGATGGAGA"),
"contig_5": DNASequence("contig_5", "CCCGGGTAGATAGGAAA"),
}
nucl = DNASequenceCollection().from_dict(nucl_fasta)

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
	"evalue": 1e-22,
	"lineage": ["Eukaryota", "Fungi", "Ascomycota", "Pezizomycetes", "Pezizales", "Pezizoaceae", "Peziza", "species"],
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
	"evalue": 1e-22,
	"lineage": ['Bacteria', 'Proteobacteria', 'Alphaproteobacteria', 'Rhizobiales', 'Rhizobiaceae', 'Sinorhizobium/Ensifer group', 'Sinorhizobium'],
	"pertinent_taxlvl": "Proteobacteria",
	"parse_lineage": "nontarget"
	}]

infotable_ldict_raw_parsed = [
	{
	'contig': 'NODE_6371', 
	'gc': 0.5928352992573176, 
	'coverage': 92.818721, 
	'pid': 'g96', 
	'length': 521, 
	'sseqid': 'sp|P55504|Y4JD_SINFN', 
	'sp_os': 'Sinorhizobium fredii (strain NBRC 101917 / NGR234)', 
	'desc': 'Uncharacterized protein y4jD OS=Sinorhizobium fredii (strain NBRC 101917 / NGR234) OX=394 GN=NGR_a03120 PE=3 SV=1', 
	'evalue': '1.21e-79', 
	'lineage': '	Bacteria;      Proteobacteria; Alphaproteobacteria;Rhizobiales,some group;   Rhizobiaceae, a sub group;	Sinorhizobium/Ensifer group;Sinorhizobium'
	}
]

initial_colnames = ["contig","gc","coverage","pid","length","sseqid","sp_os","desc","evalue"]
reload_colnames = ["contig","gc","coverage","pid","length","sseqid","sp_os","desc","evalue","lineage","pertinent_taxlvl","parse_lineage"]

infotable_reload_tsv = io.StringIO(
	'\n'.join([
		'\t'.join([str(x) for x in infotable_ldict_reload[0].values()]),
		'\t'.join([str(x) for x in infotable_ldict_reload[1].values()])
		])
	)
infotable_reload_tsv.seek(0)

def test___init__():
	it = InfoTable()
	assert isinstance(it, InfoTable)
	assert it.df.shape == (0, 12)
	assert it.tar is None
	assert it.ex is None
	assert it.children == []
	assert it.parent is None

	it = InfoTable(target_taxa = {"target": "Fungi", "exception": "Ascomycota"})
	assert isinstance(it, InfoTable)
	assert it.df.shape == (0, 12)
	assert it.tar == "Fungi"
	assert it.ex == "Ascomycota"
	assert it.children == []
	assert it.parent is None

def test_spawn_child():
	it = InfoTable()
	it.populate(infotable_ldict_reload_multiprotein)
	nested = it.df[it.df.contig == "contig_2"]
	child = it.spawn_child("child1", nested)
	assert child.df.shape == (3, 12)
	assert child.df.iloc[0,:].contig == "contig_2"
	assert len(it.children) == 1
	assert all([isinstance(x, InfoTable) for x in it.children])

@pytest.mark.skip
def test_iter_descendants():
	pass

def test_target_filter():
	it = InfoTable()
	it.populate(infotable_ldict_reload_multiprotein)
	child = it.target_filter()
	assert it.df.shape == (5, 12)
	assert child.df.shape == (3, 12)
	assert child.df.contig.to_list() == ["contig_1", "contig_2", "contig_2"]

def test_collect_unclassifieds():
	it = InfoTable()
	it.populate(infotable_ldict_reload_multiprotein)
	unclassifieds = it.collect_unclassifieds(nucl)
	assert unclassifieds.df.shape == (1, 4)
	assert unclassifieds.df.iloc[0,:].contig == "contig_5"

def test_tidy ():
	it = InfoTable()
	it.populate(infotable_ldict_raw_parsed)
	it.tidy(taxlvl_idx = 1)
	assert it.df.iloc[0,:].lineage == ["Bacteria", "Proteobacteria", "Alphaproteobacteria", "Rhizobiales_some group", "Rhizobiaceae_a sub group", "Sinorhizobium/Ensifer group","Sinorhizobium"]
	assert it.df.iloc[0,:].pertinent_taxlvl == "Proteobacteria"

def test_rfilter():
	it = InfoTable()
	it.populate(infotable_ldict_reload_multiprotein)
	window = (0.4, 0.47)
	feature = "gc"
	child = it.rfilter(feature, window)
	assert it.df.shape == (5, 12)
	assert child.df.shape == (3, 12)
	assert child.df.gc.to_list() == [0.45, 0.45, 0.45]

def test_rfilter_inplace_from_parent ():
	it = InfoTable()
	it.populate(infotable_ldict_reload_multiprotein)
	window = (0.4, 0.47)
	feature = "gc"
	child = it.rfilter(feature, window)
	child.rfilter_inplace_from_parent(feature = feature, window = (0.49, 0.60))
	assert it.df.shape == (5, 12)
	assert child.df.shape == (2, 12)
	assert child.df.gc.to_list() == [0.5, 0.55]

def test_reset_from_parent ():
	it = InfoTable()
	it.populate(infotable_ldict_reload_multiprotein)
	window = (0.4, 0.47)
	feature = "gc"
	child = it.rfilter(feature = feature, window = window)
	assert child.df.shape == (3, 12)
	child.reset_from_parent()
	assert it.df.shape == (5,12)
	assert child.df.shape == (5, 12)

def test_reset_from_parent_fail_no_parent ():
	it = InfoTable()
	it.populate(infotable_ldict_reload_multiprotein)
	with pytest.raises(AssertionError):
		it.reset_from_parent()

def test_tnt_population ():
	it = InfoTable()
	it.populate(infotable_ldict_reload_multiprotein)
	tnt_pop = it.tnt_population()
	assert (tnt_pop.target, tnt_pop.nontarget) == (3.0, 2.0)

def test_summary_stats ():
	it = InfoTable()
	it.populate(infotable_ldict_reload)
	gc_ss = it.summary_stats("gc")
	cov_ss = it.summary_stats("coverage")
	assert (gc_ss.mean, gc_ss.std, gc_ss.min, gc_ss.max) == (0.475, 0.024999999999999994, 0.45, 0.5)
	assert (cov_ss.mean, cov_ss.std, cov_ss.min, cov_ss.max) == (23, 22.5, 0.5, 45.5)

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
