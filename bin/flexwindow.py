import pandas as pd
import numpy as np
import operator
import sys

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
        print "Final window:",final_window[0],"-->",final_window[1]
        return points

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
        window = info_table[(info_table['gc'] >= sel[0]) & (info_table['gc'] <= sel[1])]
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
        print "Final window:",opt_window_neg,"-->",opt_window_pos
        return points

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
        print "Final window:",final_window[0],"-->",final_window[1]
        return points

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
        print "Final window:",opt_window_neg,"-->",opt_window_pos
        return points

    else:
        return final_window

def get_window_table (info_table, window, param):
    assert type(info_table) is pd.core.frame.DataFrame, "Input info_table to parse_infotable() must be a pandas dataframe."
    assert type(window) is tuple, "Arg window must be a tuple."
    assert param in info_table.columns, "Specified parameter is not present in input info_table."

    #print window[0],window[1]
    window_table = info_table[(info_table[param] >= window[0]) & (info_table[param] <= window[1])]
    return window_table

WindowFunc = {
        'gc0': None,
        'gc1': calc_gc_window_1tailed,
        'gc2': calc_gc_window_2tailed,
        'co0': None,
        'co1': calc_coverage_window_1tailed,
        'co2': calc_coverage_window_2tailed
        }
code_to_col = {
        'gc': 'gc',
        'co': 'coverage'
        }

class FlexibleSelectionWindow(object):
    def __init__ (self, expPat):
        self.expPat = expPat
        self.firstaxis = code_to_col[self.expPat[0:2]]
        self.secondaxis = code_to_col[self.expPat[3:5]]

        self.gc = ()
        self.coverage = ()

        self.tp = None
        self.ntp = None
        self.Wtable = None
    def show(self):
        outstr = "Type: {}\np(target): {}\np(nontarget): {}\nGC Range: {}\nCoverage Range: {}\n".format(self.expPat, self.tp, self.ntp, self.gc, self.coverage)
        return outstr

    def calculate(self, info_table, target_taxa, inc_factor=0.01, plot=False):
        self.step1 = self.expPat[0:3]
        self.step2 = self.expPat[3:6]

        ## Step 1
        if WindowFunc[self.step1] is None:
            axis1 = get_numeric_range(info_table.df[self.firstaxis])
        else:
            axis1 = WindowFunc[self.step1](info_table.df, target_taxa, inc_factor, plot)
        axis1_table = get_window_table(info_table.df, axis1, self.firstaxis)

        ## Step 2
        if WindowFunc[self.step2] is None:
            axis2 = get_numeric_range(info_table.df[self.secondaxis])
        else:
            axis2 = WindowFunc[self.step2](axis1_table, target_taxa, inc_factor, plot)
        self.Wtable = get_window_table(axis1_table, axis2, self.secondaxis)

        limits2D = {
                self.firstaxis: axis1,
                self.secondaxis: axis2
                }

        self.gc = (round(limits2D['gc'][0],4), round(limits2D['gc'][1],4))
        self.coverage = (round(limits2D['coverage'][0],4), round(limits2D['coverage'][1],4))

        ### Assess target/nontarget proportions
        all_target = float(info_table.df[info_table.df.parse_lineage == 'target'].shape[0])
        target_in_window = float(self.Wtable[self.Wtable.parse_lineage == 'target'].shape[0])
        nontarget_in_window = float(self.Wtable.shape[0] - target_in_window)
        if nontarget_in_window == 0.0:
            all_nontarget = 1.0
        else:
            all_nontarget = float(info_table.df.shape[0] - all_target)

        self.tp = round( (target_in_window / all_target),4 )
        self.ntp = round( (nontarget_in_window / all_nontarget),4 )



