import pytest
import numpy as np
import os
import inspect
from scgid.flexwindow import FlexibleSelectionWindow, WindowManager
from scgid.infotable import InfoTable

tests_dir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))

infotable_ldict_reload_multiprotein = [{
	"contig": "contig_1",
	"gc": 0.68,
	"coverage": 78.3,
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
	"contig": "contig_2",
	"gc": 0.65,
	"coverage": 80.2,
	"pid": "g0002",
	"length": 500,
	"sseqid": "ACC_02",
	"sp_os": "Pipcy3_1",
	"desc": "A description of the protein hit.",
	"evalue": 1e-99,
	"lineage": ["Eukaryota", "Viridiplantae", "Chlorophyta", "Chlorophyceae", "Chlamydomonadales", "Chlamydomonadaceae", "Chlamydomonas"],
	"pertinent_taxlvl": "Viridiplantae",
	"parse_lineage": "nontarget"
	},
	{
	"contig": "contig_3",
	"gc": 0.49,
	"coverage": 22.5,
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
	"contig": "contig_7",
	"gc": 0.55,
	"coverage": 23.0,
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
	"contig": "contig_4",
	"gc": 0.45,
	"coverage": 22.5,
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
	"contig": "contig_5",
	"gc": 0.47,
	"coverage": 25.6,
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
	"contig": "contig_6",
	"gc": 0.59,
	"coverage": 35.0,
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

sit = InfoTable()
sit.populate(infotable_ldict_reload_multiprotein)

fit = InfoTable()
fit.load(os.path.join(tests_dir, "data", "infotable.tsv"))

def calculate_expectations(it, round_func):

	gc_stats = it.target_filter().summary_stats("gc")
	gc_step_size = 0.01*gc_stats.std
	gc_pos_steps = round_func((gc_stats.max - gc_stats.mean)/gc_step_size)
	gc_neg_steps = round_func((gc_stats.mean - gc_stats.min)/gc_step_size)
	gc_steps = [gc_neg_steps, gc_pos_steps]

	coverage_stats = it.target_filter().summary_stats("coverage")
	coverage_step_size = 0.01*coverage_stats.std
	coverage_pos_steps = round_func((coverage_stats.max - coverage_stats.mean)/coverage_step_size)
	coverage_neg_steps = round_func((coverage_stats.mean - coverage_stats.min)/coverage_step_size)
	coverage_steps = [coverage_neg_steps, coverage_pos_steps]

	stats = {
		"gc": gc_stats,
		"coverage": coverage_stats
	}

	bounds = {
		"symm": {
			"gc": (
				gc_stats.mean - min(gc_steps)*gc_step_size, 
				gc_stats.mean + min(gc_steps)*gc_step_size
				),
			"coverage": (
				coverage_stats.mean - min(coverage_steps)*coverage_step_size, 
				coverage_stats.mean + min(coverage_steps)*coverage_step_size
				),
		},
		"asymm": {
			"gc": (
				gc_stats.mean - gc_steps[0]*gc_step_size, 
				gc_stats.mean + gc_steps[1]*gc_step_size
				),
			"coverage": (
				coverage_stats.mean - coverage_steps[0]*coverage_step_size, 
				coverage_stats.mean + coverage_steps[1]*coverage_step_size
				),
		}
	}

	return stats, bounds

#stats, bounds = calculate_expectations(infotable_ldict_reload_multiprotein)
#print (bounds)

'''
def generate_windows():
'''

def test_FlexibleSelectionWindow___init__ ():
	fsw = FlexibleSelectionWindow("gc2co2")
	assert isinstance(fsw, FlexibleSelectionWindow)
	assert fsw.expPat == "gc2co2"

@pytest.mark.skip
def test_FlexibleSelectionWindow_show():
	pass

@pytest.mark.parametrize("it", [sit])
def test_FlexibleSelectionWindow_calc_1d_window_symm_excludeLimits(it):
	stats, bounds = calculate_expectations(it, round_func = np.floor)
	fsw = FlexibleSelectionWindow("gc2co2")
	for feature in ["gc", "coverage"]:
		window = fsw.calc_1d_window_symm(it, feature).window
		lowerb, upperb = bounds["symm"][feature]
		assert window == (lowerb, upperb)
		assert stats[feature].min <= window[0]
		assert stats[feature].max >= window[1]

@pytest.mark.parametrize("it", [sit])
def test_FlexibleSelectionWindow_calc_1d_window_asymm_excludeLimits(it):
	stats, bounds = calculate_expectations(it, round_func = np.floor)
	fsw = FlexibleSelectionWindow("gc2co2")
	for feature in ["gc", "coverage"]:
		window = fsw.calc_1d_window_asymm(it, feature).window
		lowerb, upperb = bounds["asymm"][feature]
		assert window == (lowerb, upperb)
		assert stats[feature].min <= window[0]
		assert stats[feature].max >= window[1]

@pytest.mark.parametrize("it", [sit])
def test_FlexibleSelectionWindow_calc_1d_window_symm_includeLimits(it):
	stats, bounds = calculate_expectations(it, round_func = np.ceil)
	fsw = FlexibleSelectionWindow("gc2co2")
	for feature in ["gc", "coverage"]:
		window = fsw.calc_1d_window_symm(it, feature, include_limits = True).window
		lowerb, upperb = bounds["symm"][feature]
		assert window == (lowerb, upperb)
		assert (stats[feature].min - stats[feature].std) <= window[0]
		assert (stats[feature].max + stats[feature].std) >= window[1]

@pytest.mark.parametrize("it", [sit])
def test_FlexibleSelectionWindow_calc_1d_window_asymm_includeLimits(it):
	stats, bounds = calculate_expectations(it, round_func = np.ceil)
	fsw = FlexibleSelectionWindow("gc2co2")
	for feature in ["gc", "coverage"]:
		window = fsw.calc_1d_window_asymm(it, feature, include_limits = True).window
		lowerb, upperb = bounds["asymm"][feature]
		assert window == (lowerb, upperb)
		assert (stats[feature].min - stats[feature].std) <= window[0]
		assert (stats[feature].max + stats[feature].std) >= window[1]

@pytest.mark.parametrize("it", [sit])
def test_FlexibleSelectionWindow_calculate(it):
	stats, bounds = calculate_expectations(it, round_func = np.floor)
	first = "gc"
	second = "coverage"
	for pat in ["gc1co1", "gc1co2"]:
		fsw = FlexibleSelectionWindow(pat)
		first_window = fsw.calc_1d_window_symm(it, first).window
		fsw.calculate(it)

		first_range = fsw.stats()[first]
		intermediate = it.rfilter(first, first_range)
		r2_range = (min(intermediate.df[second]), max(intermediate.df[second]))

		assert first_range == bounds["symm"][first]
		assert fsw.stats()[second] == fsw.calc_1d_window_symm(intermediate, second).window
		assert fsw.stats()[second][0] >= r2_range[0]
		assert fsw.stats()[second][1] <= r2_range[1]

	for pat in ["gc2co1", "gc2co2"]:
		fsw = FlexibleSelectionWindow(pat)
		first_window = fsw.calc_1d_window_asymm(it, first).window
		fsw.calculate(it)

		first_range = fsw.stats()[first]
		intermediate = it.rfilter(first, first_range)
		r2_range = (min(intermediate.df[second]), max(intermediate.df[second]))

		assert first_range == bounds["asymm"][first]
		assert fsw.stats()[second] == fsw.calc_1d_window_asymm(intermediate, second).window
		assert fsw.stats()[second][0] >= r2_range[0]
		assert fsw.stats()[second][1] <= r2_range[1]

	first = "coverage"
	second = "gc"
	for pat in ["co1gc1", "co1gc2"]:
		fsw = FlexibleSelectionWindow(pat)
		first_window = fsw.calc_1d_window_symm(it, first).window
		fsw.calculate(it)

		first_range = fsw.stats()[first]
		intermediate = it.rfilter(first, first_range)
		r2_range = (min(intermediate.df[second]), max(intermediate.df[second]))

		assert first_range == bounds["symm"][first]
		assert fsw.stats()[second] == fsw.calc_1d_window_symm(intermediate, second).window
		assert fsw.stats()[second][0] >= r2_range[0]
		assert fsw.stats()[second][1] <= r2_range[1]

	for pat in ["co2gc1", "co2gc2"]:
		fsw = FlexibleSelectionWindow(pat)
		first_window = fsw.calc_1d_window_asymm(it, first).window
		fsw.calculate(it)
		
		first_range = fsw.stats()[first]
		intermediate = it.rfilter(first, first_range)
		r2_range = (min(intermediate.df[second]), max(intermediate.df[second]))

		assert first_range == bounds["asymm"][first]
		assert fsw.stats()[second] == fsw.calc_1d_window_asymm(intermediate, second).window
		assert fsw.stats()[second][0] >= r2_range[0]
		assert fsw.stats()[second][1] <= r2_range[1]

@pytest.mark.parametrize("it", [sit])
def test_FlexibleSelectionWindow_stats(it):
	fsw = FlexibleSelectionWindow("gc2co2")
	fsw.calculate(it)
	window_stats = fsw.stats()
	assert window_stats == {
		'expPat': fsw.expPat,
        'gc': fsw.gc_range,
        'coverage': fsw.coverage_range,
        'tp': fsw.tp,
        'ntp': fsw.ntp,
        'gc_width': fsw.gc_range[1] - fsw.gc_range[0],
        'co_width': fsw.coverage_range[1] - fsw.coverage_range[0],
        }

@pytest.mark.skip
def test_FlexibleSelectionWindow_to_pdf():
	pass

patterns = [
            'gc0co0',
            'gc0co1',
            'gc0co2',
            'gc1co1',
            'gc1co2',
            'gc2co1',
            'gc2co2',
            'co0gc1',
            'co0gc2',
            'co1gc1',
            'co1gc2',
            'co2gc1',
            'co2gc2'
            ]

@pytest.mark.parametrize("it", [sit])
def test_WindowManager_pick_below_thresh (it):
	s = 0.05
	wm = WindowManager(infotable = it, patterns = patterns)
	assert isinstance(wm, WindowManager)
	assert len(wm.windows) == 13

	below_thresh = wm.window_frame[wm.window_frame.ntp <= float(s)]
	max_target = below_thresh[below_thresh.tp == below_thresh.tp.max()]
	largest = max_target[max_target.sqfootage == max_target.sqfootage.max()]
	best_indep = largest.iloc[0,:]

	best_window = wm.pick(0.05)

	assert best_window.ntp == best_indep.ntp
	assert best_window.tp == best_indep.tp
	assert best_window.gc_range[1] - best_window.gc_range[0] == best_indep.gc_width
	assert best_window.coverage_range[1] - best_window.coverage_range[0] == best_indep.co_width
	assert best_window.expPat == best_indep.expPat

@pytest.mark.parametrize("it", [sit])
def test_WindowManager_pick_above_thresh (it):
	s = -1.0
	wm = WindowManager(infotable = it, patterns = patterns)
	assert isinstance(wm, WindowManager)
	assert len(wm.windows) == 13

	best_window = wm.pick(0.05)

	sqfootage = (best_window.coverage_range[1] - best_window.coverage_range[0]) * (best_window.gc_range[1] - best_window.gc_range[0])

	assert np.divide(best_window.tp, best_window.ntp) == wm.window_frame.quotient.max()
	assert sqfootage == wm.window_frame[wm.window_frame.quotient == wm.window_frame.quotient.max()].sqfootage.max()