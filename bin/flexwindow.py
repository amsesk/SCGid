import pandas as pd
import numpy as np
import operator
import sys
from collections import namedtuple

### Definitions
PltOut = namedtuple("PltOut", ["axis","mean","window","points"])
Points = namedtuple("Points",["step","value","tp","ntp","tradeoff"])


def get_numeric_range(pd_series):
    return (pd_series.min(), pd_series.max())
#%%
def calc_gc_window_1tailed (info_table, target_taxa, inc_factor = 0.01, plot = False):
    assert type(info_table) is pd.core.frame.DataFrame, "Input info_table to parse_infotable() must be a pandas dataframe."

    ### Get summary stats needed for determining gc window from info_table ###
    target = info_table.loc[info_table['parse_lineage'] == 'target']
    target_mean_gc = np.mean(target.gc)
    target_std_gc = np.std(target.gc)
    target_gc_range = (target.gc.min(),target.gc.max())
    all_target = float(target.shape[0])
    all_nontarget = float(info_table.shape[0] - all_target)

    ## Account for case where there is no nontarget (should be really rare), avoid division by zero error
    if all_nontarget == 0:
        all_nontarget = 1

    increment = float(target_std_gc)*inc_factor
    points = {
            "target": [],
            "nontarget": []
            }

    while (target_mean_gc - increment >= target_gc_range[0]) and (target_mean_gc + increment <= target_gc_range[1]):
        sel = (target_mean_gc - increment, target_mean_gc + increment)
        window = info_table[(info_table['gc'] >= sel[0]) & (info_table['gc'] <= sel[1])]
        target_in_window = float(window.loc[window['parse_lineage'] == 'target'].shape[0])
        nontarget_in_window = float(window.shape[0] - target_in_window)
        points["target"].append(target_in_window/all_target)
        points["nontarget"].append(nontarget_in_window/all_nontarget)
        increment += target_std_gc*inc_factor

    diff = map(operator.sub,points["target"],points["nontarget"])
    tradeoffs = map(operator.mul,points["target"],diff)
    points["tradeoffs"] = tradeoffs

    num_maxes = [i for i in points["tradeoffs"] if i == max(points["tradeoffs"])]

    if len(num_maxes) > 1:
        ## get max tradeoff at highest step (ie for asymtotic tradeoff curves)
        idx = 0
        indices_of_max = []
        for val in points["tradeoffs"]:
            if val == max(points["tradeoffs"]):
                indices_of_max.append(idx)
            idx+=1
        best_step_idx = max(indices_of_max)
    else:
        best_step_idx = points["tradeoffs"].index(max(points["tradeoffs"]))

    opt_tail = target_std_gc * inc_factor * (best_step_idx+1)
    final_window = (target_mean_gc - opt_tail, target_mean_gc + opt_tail)

    if plot is True:
        return PltOut("gc", final_window, points)

    else:
        return final_window
#%%
def calc_gc_window_2tailed (info_table, target_taxa, inc_factor = 0.01, plot = False):
    target = info_table.loc[info_table['parse_lineage'] == 'target']
    target_mean_gc = np.mean(target.gc)
    target_std_gc = np.std(target.gc)
    target_gc_range = (target.gc.min(),target.gc.max())
    all_target = float(target.shape[0])
    all_nontarget = float(info_table.shape[0] - all_target)

    ## Account for case where there is no nontarget (should be really rare), avoid division by zero error
    if all_nontarget == 0:
        all_nontarget = 1

    increment = float(target_std_gc)*inc_factor
    points = {
            "neg":{
                    "target": [],
                    "nontarget": []
                    },
            "pos":{
                    "target": [],
                    "nontarget": []
                    }
            }
    while target_mean_gc - increment >= target_gc_range[0]:
        sel = (target_mean_gc - increment, target_mean_gc)
        window = info_table.loc[(info_table['gc'] >= sel[0]) & (info_table['gc'] <= sel[1])]
        target_in_window = float(window.loc[window['parse_lineage'] == 'target'].shape[0])
        nontarget_in_window = float(window.shape[0] - target_in_window)
        points["neg"]["target"].append(target_in_window/all_target)
        points["neg"]["nontarget"].append(nontarget_in_window/all_nontarget)
        increment += target_std_gc*inc_factor

    ### Reset increment ###
    increment = float(target_std_gc)*inc_factor

    while target_mean_gc + increment <= target_gc_range[1]:
        sel = (target_mean_gc, target_mean_gc + increment)
        window = info_table[(info_table['gc'] >= sel[0]) & (info_table['gc'] <= sel[1])]
        target_in_window = float(window.loc[window['parse_lineage'] == 'target'].shape[0])
        nontarget_in_window = float(window.shape[0] - target_in_window)
        points["pos"]["target"].append(target_in_window/all_target)
        points["pos"]["nontarget"].append(nontarget_in_window/all_nontarget)
        increment += target_std_gc*inc_factor


    diff_neg = map(operator.sub,points["neg"]["target"],points["neg"]["nontarget"])
    diff_pos = map(operator.sub,points["pos"]["target"],points["pos"]["nontarget"])

    tradeoffs_neg = map(operator.mul,points["neg"]["target"],diff_neg)
    tradeoffs_pos = map(operator.mul,points["pos"]["target"],diff_pos)

    points["neg"]["tradeoffs"] = tradeoffs_neg
    points["pos"]["tradeoffs"] = tradeoffs_pos

    num_maxes_neg = [i for i in points["neg"]["tradeoffs"] if i == max(points["neg"]["tradeoffs"])]
    num_maxes_pos = [i for i in points["pos"]["tradeoffs"] if i == max(points["pos"]["tradeoffs"])]

    ## make sure that the following FOR loop is necessary in the negative direction (ie asymptotic tradeoff curves)
    if len(num_maxes_neg) > 1:
        ## get max tradeoff at highest step in NEGATIVE direction (ie for asymtotic tradeoff curves)
        idx = 0
        indices_of_max_neg = []
        for val in points["neg"]["tradeoffs"]:
            if val == max(points["neg"]["tradeoffs"]):
                indices_of_max_neg.append(idx)
            idx+=1
        best_step_idx_neg = max(indices_of_max_neg)
    else:
        best_step_idx_neg = points["neg"]["tradeoffs"].index(max(points["neg"]["tradeoffs"]))
    

    if len(num_maxes_pos) > 1:
        ## get max tradeoff at highest step in POSITIVE direction (ie for asymtotic tradeoff curves)
        idx = 0
        indices_of_max_pos = []
        for val in points["pos"]["tradeoffs"]:
            if val == max(points["pos"]["tradeoffs"]):
                indices_of_max_pos.append(idx)
            idx+=1
        best_step_idx_pos = max(indices_of_max_pos)
    else:
        best_step_idx_pos = points["pos"]["tradeoffs"].index(max(points["pos"]["tradeoffs"]))
    
    opt_tail_neg = target_std_gc * inc_factor * (best_step_idx_neg+1)
    opt_tail_pos = target_std_gc * inc_factor * (best_step_idx_pos+1)

    opt_window_neg = target_mean_gc - opt_tail_neg
    opt_window_pos = target_mean_gc + opt_tail_pos

    final_window = (opt_window_neg, opt_window_pos)

    if plot is True:
        return PltOut("gc", final_window, points)

    else:
        return final_window
#%%

### When should there NOT be a coverage cut-off and how do I code that? ###

def calc_coverage_window_1tailed (info_table, target_taxa, inc_factor = 0.01, plot = False):
    target = target = info_table.loc[info_table['parse_lineage'] == 'target']
    target_mean_cov = np.mean(target.coverage)
    target_std_cov = np.std(target.coverage)
    target_cov_range = (target.coverage.min(),target.coverage.max())
    all_target = float(target.shape[0])
    all_nontarget = float(info_table.shape[0] - all_target)
    increment = float(target_std_cov) * inc_factor

    ## Account for case where there is no nontarget (should be really rare), avoid division by zero error
    if all_nontarget == 0:
        all_nontarget = 1

    points = {
            "target": [],
            "nontarget": []
            }
    while (target_mean_cov - increment >= target_cov_range[0]) and (target_mean_cov + increment <= target_cov_range[1]):
        sel = (target_mean_cov - increment, target_mean_cov + increment)
        window = info_table[(info_table['coverage'] >= sel[0]) & (info_table['coverage'] <= sel[1])]
        target_in_window = float(window.loc[window['parse_lineage'] == 'target'].shape[0])
        nontarget_in_window = float(window.shape[0] - target_in_window)
        points["target"].append(target_in_window/all_target)
        points["nontarget"].append(nontarget_in_window/all_nontarget)
        increment += target_std_cov*inc_factor

    diff = map(operator.sub,points["target"],points["nontarget"])
    tradeoffs = map(operator.mul,points["target"],diff)
    points["tradeoffs"] = tradeoffs

    num_maxes = [i for i in points["tradeoffs"] if i == max(points["tradeoffs"])]

    if len(num_maxes) > 1:
        ## get max tradeoff at highest step (ie for asymtotic tradeoff curves)
        idx = 0
        indices_of_max = []
        for val in points["tradeoffs"]:
            if val == max(points["tradeoffs"]):
                indices_of_max.append(idx)
            idx+=1
        best_step_idx = max(indices_of_max)
    else:
        best_step_idx = points["tradeoffs"].index(max(points["tradeoffs"]))

    opt_tail = target_std_cov * inc_factor * (best_step_idx+1)
    final_window = (target_mean_cov - opt_tail, target_mean_cov + opt_tail)

    if plot is True:
        return PltOut("coverage", final_window, points)

    else:
        return final_window
    
#%%         
def calc_coverage_window_2tailed (info_table, target_taxa, inc_factor = 0.01, plot = False):    
    target = info_table.loc[info_table['parse_lineage'] == 'target']
    target_mean_cov = np.mean(target.coverage)
    target_std_cov = np.std(target.coverage)
    target_cov_range = (target.coverage.min(),target.coverage.max())
    increment = float(target_std_cov)*inc_factor
    all_target = float(target.shape[0])
    all_nontarget = float(info_table.shape[0] - all_target)

    ## Account for case where there is no nontarget (should be really rare), avoid division by zero error
    if all_nontarget == 0:
        all_nontarget = 1

    points = {
            "neg":{
                    "target": [],
                    "nontarget": []
                    },
            "pos":{
                    "target": [],
                    "nontarget": []
                    }
            }

    while target_mean_cov - increment >= target_cov_range[0]:
        sel = (target_mean_cov - increment, target_mean_cov)
        window = info_table[(info_table['coverage'] >= sel[0]) & (info_table['coverage'] <= sel[1])]
        target_in_window = float(window.loc[window['parse_lineage'] == 'target'].shape[0])
        nontarget_in_window = float(window.shape[0] - target_in_window)
        points["neg"]["target"].append(target_in_window/all_target)
        points["neg"]["nontarget"].append(nontarget_in_window/all_nontarget)
        increment += target_std_cov*inc_factor
        #print target_mean_cov-increment, points["neg"]["target"][-1]

    ### Reset increment ###
    increment = float(target_std_cov)*inc_factor

    while target_mean_cov + increment <= target_cov_range[1]:
        sel = (target_mean_cov, target_mean_cov + increment)
        window = info_table[(info_table['coverage'] >= sel[0]) & (info_table['coverage'] <= sel[1])]
        target_in_window = float(window.loc[window['parse_lineage'] == 'target'].shape[0])
        nontarget_in_window = float(window.shape[0] - target_in_window)
        ## Account for case where there is no nontarget (should be really rare), avoid division by zero error
        if all_nontarget == 0:
            all_nontarget = 1

        points["pos"]["target"].append(target_in_window/all_target)
        points["pos"]["nontarget"].append(nontarget_in_window/all_nontarget)


        increment += target_std_cov*inc_factor


    diff_neg = map(operator.sub,points["neg"]["target"],points["neg"]["nontarget"])
    diff_pos = map(operator.sub,points["pos"]["target"],points["pos"]["nontarget"])
    #print len(points["pos"]["target"])
    #print len(diff_pos)

    tradeoffs_neg = map(operator.mul,points["neg"]["target"],diff_neg)
    tradeoffs_pos = map(operator.mul,points["pos"]["target"],diff_pos)

    points["neg"]["tradeoffs"] = tradeoffs_neg
    points["pos"]["tradeoffs"] = tradeoffs_pos

    ## how many points represent the maximum trade-off?
    num_maxes_neg = [i for i in points["neg"]["tradeoffs"] if i == max(points["neg"]["tradeoffs"])]
    num_maxes_pos = [i for i in points["pos"]["tradeoffs"] if i == max(points["pos"]["tradeoffs"])]

    ## make sure that the following FOR loop is necessary in the negative direction (ie asymptotic tradeoff curves)
    if len(num_maxes_neg) > 1:
        ## get max tradeoff at highest step in NEGATIVE direction (ie for asymptotic tradeoff curves)
        idx = 0
        indices_of_max_neg = []
        for val in points["neg"]["tradeoffs"]:
            if val == max(points["neg"]["tradeoffs"]):
                indices_of_max_neg.append(idx)
            idx+=1
        best_step_idx_neg = max(indices_of_max_neg)

    else:
        best_step_idx_neg = points["neg"]["tradeoffs"].index(max(points["neg"]["tradeoffs"]))

    ## get max tradeoff at highest step in POSITIVE direction (ie for asymtotic tradeoff curves)
    if len(num_maxes_pos) > 1:
        ## get max tradeoff at highest step in POSITIVE direction (ie for asymptotic tradeoff curves)
        idx = 0
        indices_of_max_pos = []
        for val in points["pos"]["tradeoffs"]:
            if val == max(points["pos"]["tradeoffs"]):
                indices_of_max_pos.append(idx)
            idx+=1
        best_step_idx_pos = max(indices_of_max_pos)

    else:
        best_step_idx_pos = points["pos"]["tradeoffs"].index(max(points["pos"]["tradeoffs"]))
    #print best_step_idx_pos

    #best_step_idx_neg = points["neg"]["tradeoffs"].index(max(points["neg"]["tradeoffs"]))
    #best_step_idx_pos = points["pos"]["tradeoffs"].index(max(points["pos"]["tradeoffs"]))

    opt_tail_neg = target_std_cov * inc_factor * (best_step_idx_neg+1)
    opt_tail_pos = target_std_cov * inc_factor * (best_step_idx_pos+1)

    #print opt_tail_pos

    opt_window_neg = target_mean_cov - opt_tail_neg
    opt_window_pos = target_mean_cov + opt_tail_pos

    #print opt_window_pos

    final_window = (opt_window_neg, opt_window_pos)

    if plot is True:
        return PltOut("coverage", final_window, points)

    else:
        return final_window
#%%
def get_window_table (info_table, window, param):
    assert type(info_table) is pd.core.frame.DataFrame, "Input info_table to parse_infotable() must be a pandas dataframe."
    assert type(window) is tuple, "Arg window must be a tuple."
    assert param in info_table.columns, "Specified parameter is not present in input info_table."

    #print window[0],window[1]
    window_table = info_table[(info_table[param] >= window[0]) & (info_table[param] <= window[1])]
    return window_table


#%%
def calc_1d_window_symm (it, axis, inc_factor = 0.01, plot = False):
    
    flextbl = it.spawn_child("flextable", it.df)
    target = it.target_filter()
    axis_ss = target.summary_stats(axis)
    pdist = axis_ss.max - axis_ss.mean
    ndist = axis_ss.mean - axis_ss.min
    dist = min( [pdist,ndist] )
    step_size = inc_factor * axis_ss.std
    population_whole = it.tnt_population()
    steps_to_take = int(np.floor(dist/step_size))
    
    points = np.ndarray(shape=(steps_to_take+1, 6), dtype='float', order='C')
    
    final_window = ()
    
    for s in range(0,steps_to_take+1):
        step_window = (axis_ss.mean - s*step_size, axis_ss.mean + s*step_size)
        flextbl.rfilter_inplace_from_parent(axis, step_window)
        population_now = flextbl.tnt_population()
        tp = population_now.target/population_whole.target
        ntp = population_now.nontarget/population_whole.nontarget
        pointrow = np.array([
                s,
                step_window[0],
                step_window[1],
                tp,
                ntp,
                tp*(tp-ntp)
                ])
        points[s] = pointrow
    points = pd.DataFrame(points, columns=["step",
                                           "lower_{}".format(axis), 
                                           "upper_{}".format(axis),
                                           "tp",
                                           "ntp",
                                           "tradeoff"]) 
        
    maxes = points.loc[points.tradeoff == points.tradeoff.max()]
    if maxes.shape[0] > 1:
            final_window = (
                    maxes.loc[maxes["step"] == maxes["step"].max()]["lower_{}".format(axis)].item(),
                    maxes.loc[maxes["step"] == maxes["step"].max()]["upper_{}".format(axis)].item()
                    )
    else:
        final_window = (
                maxes["lower_{}".format(axis)].item(),
                maxes["upper_{}".format(axis)].item()
                )
    
    return PltOut(axis, axis_ss.mean, final_window, points)

#%%
def calc_1d_window_asymm (it, axis, inc_factor = 0.01, plot = False):
    
    flextbl = it.spawn_child("flextable", it.df)
    target = it.target_filter()
    axis_ss = target.summary_stats(axis)
    pdist = axis_ss.max - axis_ss.mean
    ndist = axis_ss.mean - axis_ss.min
    step_size = inc_factor * axis_ss.std
    population_whole = it.tnt_population()
    steps_to_take = {
            -1: int(np.floor(ndist/step_size)),
            1: int(np.floor(pdist/step_size))
            }
    points = {
            -1: np.ndarray(shape=(steps_to_take[-1]+1, 6), dtype='float', order='C'), #negative direction
            1: np.ndarray(shape=(steps_to_take[1]+1, 6), dtype='float', order='C'), #positive direction
            }
    point_column = {
            -1: "lower",
            1: "upper"
            }
    final_window = {}
    
    for d in points.keys():
        flextbl.reset_from_parent()
        for s in range(0,steps_to_take[d]+1):
            step_window = (axis_ss.mean, axis_ss.mean + d * s*step_size)
            flextbl.rfilter_inplace_from_parent(axis, step_window)
            population_now = flextbl.tnt_population()
            tp = population_now.target/population_whole.target
            ntp = population_now.nontarget/population_whole.nontarget
            pointrow = np.array([
                    s,
                    min(step_window),
                    max(step_window),
                    tp,
                    ntp,
                    tp*(tp-ntp)
                    ])
            points[d][s] = pointrow
        points[d] = pd.DataFrame(points[d], columns=["step", 
                                                      "lower_{}".format(axis),
                                                      "upper_{}".format(axis),
                                                      "tp", 
                                                      "ntp", 
                                                      "tradeoff"]) 
        
        maxes = points[d].loc[points[d].tradeoff == points[d].tradeoff.max()]
        if maxes.shape[0] > 1:
            final_window[d] = maxes.loc[maxes["step"] == maxes["step"].max()]["{}_{}".format(point_column[d], axis)].item()
        else:
            final_window[d] = maxes["{}_{}".format(point_column[d], axis)].item()
    
    return PltOut(axis, axis_ss.mean, (final_window[-1], final_window[1]), points)
#%%

WindowFunc = {
        1: calc_1d_window_symm,
        2: calc_1d_window_asymm,
        }
code_to_col = {
        'gc': 'gc',
        'co': 'coverage'
        }

class FlexibleSelectionWindow(object):
    def __init__ (self, expPat):
        self.expPat = expPat
        
        self.step1 = None
        self.firstaxis = None
        self.firstsymm = None
        
        self.step = None
        self.secondaxis = None
        self.secondsymm = None
        
        self.means = {
                "gc": None,
                "coverage": None
                }

        self.gc_range = ()
        self.gc_points = None
        self.coverage_range = ()
        self.coverage_points = None

        self.tp = None
        self.ntp = None
        self.Wtable = None
    def show(self):
        outstr = "Type: {}\np(target): {}\np(nontarget): {}\nGC Range: {}\nCoverage Range: {}\n".format(
                self.expPat, 
                round(self.tp, 4), 
                round(self.ntp,4), 
                [round(x,4) for x in self.gc_range],
                [round(x,4) for x in self.coverage_range]
                )
        return outstr

    def calculate(self, it, inc_factor=0.01, plot=False):
        self.step1 = self.expPat[0:3]
        self.firstaxis = code_to_col[self.step1[0:2]]
        self.firstsymm = int(self.step1[2])
        
        self.step2 = self.expPat[3:6]
        self.secondaxis = code_to_col[self.step2[0:2]]
        self.secondsymm = int(self.step2[2])

        ## Step 1
        axis1_range = it.summary_stats(self.firstaxis).range
        if self.firstsymm == 0:
            axis1_window = axis1_range
            axis1_table = it.rfilter(self.firstaxis, axis1_range, self.step1)
        else:
            self.firstaxis, self.means[self.firstaxis], axis1_window, axis1_points = WindowFunc[self.firstsymm](it, self.firstaxis, inc_factor, plot = True)
                
            axis1_table = it.rfilter(self.firstaxis, axis1_window, self.step1)

        ## Step 2
        axis2_range = axis1_table.summary_stats(self.secondaxis).range
        if self.secondsymm == 0:
            axis2_window = axis2_range
            self.Wtable = it.rfilter(self.secondaxis, axis2_range, self.step2) #window = axis2_range
        else:
            self.secondaxis, self.means[self.secondaxis], axis2_window, axis2_points = WindowFunc[self.secondsymm](axis1_table, self.secondaxis, inc_factor, plot = True)
                
            self.Wtable = axis1_table.rfilter(self.secondaxis, axis2_window, self.step2)
        
        limits2D = {
                self.firstaxis: axis1_window,
                self.secondaxis: axis2_window,
                }
        if plot:
            points = {
                    self.firstaxis: axis1_points,
                    self.secondaxis: axis2_points
                    }
            self.gc_points = points["gc"]
            self.coverage_points = points["coverage"]

        self.gc_range = (limits2D['gc'][0], limits2D['gc'][1])
        self.coverage_range = (limits2D['coverage'][0], limits2D['coverage'][1])

        ### Assess target/nontarget proportions
        whole_pop = it.tnt_population()
        Wtable_pop = self.Wtable.tnt_population()
        
        all_target = whole_pop.target
        all_nontarget = whole_pop.nontarget
        
        target_in_window = Wtable_pop.target
        nontarget_in_window = Wtable_pop.nontarget
        
        if all_nontarget == 0.0:
            all_nontarget = 1.0

        self.tp = target_in_window / all_target
        self.ntp = nontarget_in_window / all_nontarget


