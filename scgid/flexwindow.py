import pandas as pd
import numpy as np
import operator
import sys
import os
import logging
from collections import namedtuple
from scgid.modcomm import LoggingEntity, ErrorHandler, get_head, logger_name_gen
from scgid.library import subprocessP
from scgid.error import ModuleError

### Definitions
PltOut = namedtuple("PltOut", ["axis","mean","window","points"])
Points = namedtuple("Points",["step","value","tp","ntp","tradeoff"])

def generate_windows(it, inc_factor=0.01):
    # factorial window expansion combinations MINUS redundancy ie gc0co0 = co0gc0
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

    return WindowManager(it, patterns, inc_factor)

code_to_col = {
        'gc': 'gc',
        'co': 'coverage'
        }

class NoTargetError(ModuleError):
    def __init__(self, error_catch = True):
        super().__init__()
        self.msg = f"No target-annotated points in plot data, so nothing to do. This affects all windows. Try changing target-except groupings."
        self.errno = 27
        if error_catch:
            self.catch()

class FlexibleSelectionWindow(LoggingEntity, ErrorHandler):
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
        self.WindowFunc = {
            1: self.calc_1d_window_symm,
            2: self.calc_1d_window_asymm,
        }

        self.zero_width_axis = False

        self.gc_range = ()
        self.gc_points = None
        self.coverage_range = ()
        self.coverage_points = None

        self.n_target = None
        self.n_nontarget = None
        self.tp = None
        self.ntp = None
        self.Wtable = None

        self.head = get_head()
        self.logger = logging.getLogger( logger_name_gen() )

    def show(self):
        outstr = "Type: {}\np(target): {}\np(nontarget): {}\nGC Range: {}\nCoverage Range: {}\n".format(
                self.expPat,
                round(self.tp, 4),
                round(self.ntp,4),
                [round(x,4) for x in self.gc_range],
                [round(x,4) for x in self.coverage_range]
                )
        return outstr

    #%%
    def calc_1d_window_symm (self, it, axis, inc_factor = 0.01, plot = False, include_limits = False):

        if include_limits:
            round_func = np.ceil
        else:
            round_func = np.floor

        flextbl = it.spawn_child("flextable", it.df)
        target = it.target_filter()

        if target.df.shape[0] == 0:
            #print(axis, "AH!")
            raise NoTargetError

        axis_ss = target.summary_stats(axis)
        pdist = axis_ss.max - axis_ss.mean
        ndist = axis_ss.mean - axis_ss.min
        dist = min( [pdist,ndist] )
        step_size = inc_factor * axis_ss.std
        population_whole = it.tnt_population()

        try:
            steps_to_take = int(round_func(dist/step_size))

        except ValueError:
            self.logger.warning(f"Abnormal data. Standard deviation of target points along {axis}-axis equals 0. Setting value to 1 to avoid ZeroDivisionError.")
            steps_to_take = int(round_func(dist/1))

        points = np.ndarray(shape=(steps_to_take+1, 6), dtype='float', order='C')

        final_window = ()

        for s in range(0,steps_to_take+1):
            step_window = (axis_ss.mean - s*step_size, axis_ss.mean + s*step_size)
            flextbl.rfilter_inplace_from_parent(axis, step_window)
            population_now = flextbl.tnt_population()

            try:
                tp = population_now.target/population_whole.target
            except ZeroDivisionError:
                raise NoTargetError

            try:
                ntp = population_now.nontarget/population_whole.nontarget
            except ZeroDivisionError:
                self.logger.debug("No nontarget in window. Setting value to 1 to avoid ZeroDivsionError")
                ntp = population_now.nontarget/1

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
                        next( iter(maxes.loc[maxes["step"] == maxes["step"].max()]["lower_{}".format(axis)]), "no match"),
                        next( iter(maxes.loc[maxes["step"] == maxes["step"].max()]["upper_{}".format(axis)]), "no match")
                        )
        else:
            final_window = (
                    maxes["lower_{}".format(axis)].item(),
                    maxes["upper_{}".format(axis)].item()
                    )
        #print(axis_ss.max, axis_ss.min, steps_to_take)
        #print(axis, axis_ss.mean, final_window, points)
        return PltOut(axis, axis_ss.mean, final_window, points)

    def calc_1d_window_asymm (self, it, axis, inc_factor = 0.01, plot = False, include_limits = False):

        if include_limits:
            round_func = np.ceil
        else:
            round_func = np.floor

        flextbl = it.spawn_child("flextable", it.df)
        target = it.target_filter()

        if target.df.shape[0] == 0:
            raise NoTargetError

        axis_ss = target.summary_stats(axis)
        pdist = axis_ss.max - axis_ss.mean
        ndist = axis_ss.mean - axis_ss.min
        step_size = inc_factor * axis_ss.std
        population_whole = it.tnt_population()

        try:
            steps_to_take = {
                -1: int(round_func(ndist/step_size)),
                1: int(round_func(pdist/step_size))
                }
        except ValueError:
            print(step_size)
            self.logger.warning(f"Abnormal data. Standard deviation of target points along {axis}-axis equals 0. Setting value to 1 to avoid ZeroDivisionError.")
            steps_to_take = {
                -1: int(round_func(ndist/1)),
                1: int(round_func(pdist/1))
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

                try:
                    tp = population_now.target/population_whole.target
                except ZeroDivisionError:
                    raise NoTargetError

                try:
                    ntp = population_now.nontarget/population_whole.nontarget
                except ZeroDivisionError:
                    self.logger.debug("No nontarget in window. Setting value to 1 to avoid ZeroDivsionError")
                    ntp = population_now.nontarget/1

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
                final_window[d] = next( iter(maxes.loc[maxes["step"] == maxes["step"].max()]["{}_{}".format(point_column[d], axis)]), "no match")
            else:
                final_window[d] = maxes["{}_{}".format(point_column[d], axis)].item()

        return PltOut(axis, axis_ss.mean, (final_window[-1], final_window[1]), points)

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
            self.firstaxis, self.means[self.firstaxis], axis1_window, axis1_points = self.WindowFunc[self.firstsymm](it, self.firstaxis, inc_factor, plot = True)

            axis1_table = it.rfilter(self.firstaxis, axis1_window, self.step1)

        if axis1_window[1] - axis1_window[0] == 0:
            self.logger.warning(f"Axis 1 ({self.firstaxis}) of {self.expPat} window is 0 width. Skipping...")
            self.zero_width_axis = True
            return None

        ## Step 2
        axis2_range = axis1_table.summary_stats(self.secondaxis).range
        if self.secondsymm == 0:
            axis2_window = axis2_range
            self.Wtable = it.rfilter(self.secondaxis, axis2_range, self.step2) #window = axis2_range
        else:
            self.secondaxis, self.means[self.secondaxis], axis2_window, axis2_points = self.WindowFunc[self.secondsymm](axis1_table, self.secondaxis, inc_factor, plot = True)

            self.Wtable = axis1_table.rfilter(self.secondaxis, axis2_window, self.step2)

        limits2D = {
                self.firstaxis: axis1_window,
                self.secondaxis: axis2_window,
                }

        if axis2_window[1] - axis2_window[0] == 0:
            self.logger.warning(f"Axis 2 ({self.secondaxis}) of {self.expPat} window is 0 width. Skipping...")
            self.zero_width_axis = True
            return None

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

        if all_target == 0.0:
            raise NoTargetError

        if all_nontarget == 0.0:
            self.logger.warning(f"No nontarget in final window `{self.expPat}`. Setting to 1 to avoid ZeroDivisionError.")
            all_nontarget = 1.0

        self.n_target = int(target_in_window)
        self.n_nontarget = int(nontarget_in_window)

        self.tp = target_in_window / all_target
        self.ntp = nontarget_in_window / all_nontarget

    def stats(self):
        return {
            'expPat': self.expPat,
            'gc': self.gc_range,
            'coverage': self.coverage_range,
            'tp': self.tp,
            'ntp': self.ntp,
            'gc_width': self.gc_range[1] - self.gc_range[0],
            'co_width': self.coverage_range[1] - self.coverage_range[0],
            }

    def to_pdf(self, outdir):
        cmd = [
            os.path.join(self.head.config.get("path_to_Rscript"), "Rscript"),
            "--vanilla",
            os.path.join(self.head.config.SCGID_SCRIPTS, "gc_cov.plot.R"),
            f"{self.head.config.get('prefix')}.infotable.tsv",
            f"{self.head.config.get('prefix')}.unclassified.infotable.tsv",
            ','.join(map(str,self.gc_range)),
            ','.join(map(str,self.coverage_range)),
            os.path.join(outdir,
                f"{self.head.config.get('prefix')}.{self.expPat}.pdf"
                )
            ]

        self.logger.info(' '.join(cmd))
        subprocessP(cmd, self.head.logger)

        return 0

class WindowManager(LoggingEntity):
    def __init__(self, infotable, patterns, inc_factor = 0.01):
        self.windows = {p: FlexibleSelectionWindow(p) for p in patterns}
        self.window_frame = None
        self.infotable = infotable

        self.logger = logging.getLogger( logger_name_gen() )

        for w in self.windows.values():
            #print(f"Starting {w.expPat} with inc_factor = {inc_factor}")
            w.calculate(self.infotable, inc_factor)
            #print(f"Done with {w.expPat}")

        self.windows = {p:w for p,w in self.windows.items() if not w.zero_width_axis}

        self.logger.info(f"{len(self.windows.values())}/13 total possible windows calculated passed >0 axis width requirement. This is temporary, see fix in SCGid ToDo list on repo.")

        ldict = [w.stats() for w in self.windows.values()]
        colnames = ldict[0].keys()

        self.window_frame = pd.DataFrame(ldict, columns = colnames)
        self.window_frame = self.window_frame.assign(sqfootage=(self.window_frame.gc_width * self.window_frame.co_width))

        self.logger = logging.getLogger ( logger_name_gen() )
        self.head = get_head()

    def print_all_pdf(self, outdir):
        for w in self.windows.values():
            w.to_pdf(outdir)

    def print_all_tsv(self, outdir):
        with open(outdir, 'w') as f:
            total_target = self.infotable.df[self.infotable.df["parse_lineage"] == "target"].shape[0]
            total_nontarget = self.infotable.df[self.infotable.df["parse_lineage"] == "nontarget"].shape[0]

            lines = []
            row = "{:<10} {:<10} {:<10} {:<10} {:<10} {:<10} {:<10} {:<20} {:<20}"

            lines.append( row.format(
                "Pattern",
                "#target",
                "#total",
                "p(target)",
                "#ntarget",
                "total",
                "p(ntarget)",
                "GC Range",
                "Coverage Range"
                )
            )
            for win in self.windows.values():
                gc_range_repr = f"[{round(win.gc_range[0],4)}, {round(win.gc_range[1],4)}]"
                cov_range_repr = f"[{round(win.coverage_range[0],4)}, {round(win.coverage_range[1],4)}]"
                lines.append( row.format(
                    win.expPat,
                    win.n_target,
                    total_target,
                    round(win.tp, 4),
                    win.n_nontarget,
                    total_nontarget,
                    round(win.ntp, 4),
                    gc_range_repr,
                    cov_range_repr
                    )
                )
            output = "\n".join(lines)
            f.write(output)

            return 0

    def pick(self, stringency):

        # Calculate tp/ntp quotient and add to frame
        self.window_frame = self.window_frame.assign(
                quotient = lambda x: x.tp / x.ntp
                )

        # Sort descending by tp/ntp quotient and sqfootage
        self.window_frame.sort_values(
            by = ["quotient", "sqfootage"],
            ascending = False,
            inplace = True
            )

        below_thresh = self.window_frame[self.window_frame.ntp <= float(stringency)]

        if below_thresh.shape[0] == 0:

            self.logger.info( f"No usable window at set stringency threshold, `s = {stringency}`. Determining next best option..." )

            # Pick best window by maximizing tp/ntp quotient and window size (sorted above)
            best = self.window_frame.iloc[0,:]

            # Set Gct.config.stringency to new "best" stringency level determined above by sorting by quotient and sqfootage descending
            setattr( self.head.config, "stringency", best["ntp"] )

            if self.window_frame.quotient.value_counts()[ best.quotient ] > 1 and self.window_frame.sqfootage.value_counts()[ best.sqfootage ] > 1:
                self.logger.warning("More than one best window in terms of target:nontarget ratio and size. Arbitrarily picking 1st one.")

            return self.windows[ best.expPat ]

        else:

            # Filter for window(s) with the highest proportion of target
            max_target = below_thresh[below_thresh.tp == below_thresh.tp.max()]

            # Filter for largest window(s)
            largest = max_target[max_target.sqfootage == max_target.sqfootage.max()]

            if largest.shape[0] > 1:

                self.logger.warning( f"More than one best window under set stringency level `s = {self.head.config.get('stringency')}` and of maximal size. Arbitrarily picking 1st one." )

            best = largest.iloc[0,:]
            return self.windows [ best.expPat ]


''' DEPRECATED STUFF THAT HAS BEEN REPLICATED IN CLASSES ABOVE
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

    try:
        steps_to_take = int(np.floor(dist/step_size))
    except ZeroDivisionError:
        steps_to_take = int(np.floor(dist/1))

    points = np.ndarray(shape=(steps_to_take+1, 6), dtype='float', order='C')

    final_window = ()

    for s in range(0,steps_to_take+1):
        step_window = (axis_ss.mean - s*step_size, axis_ss.mean + s*step_size)
        flextbl.rfilter_inplace_from_parent(axis, step_window)
        population_now = flextbl.tnt_population()

        try:
            tp = population_now.target/population_whole.target
        except ZeroDivisionError:
            tp = population_now.target/1

        try:
            ntp = population_now.nontarget/population_whole.nontarget
        except:
            population_now.nontarget/1

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

            try:
                tp = population_now.target/population_whole.target
            except ZeroDivisionError:
                tp = population_now.target/1

            try:
                ntp = population_now.nontarget/population_whole.nontarget
            except:
                population_now.nontarget/1

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
'''

'''
WindowFunc = {
        1: calc_1d_window_symm,
        2: calc_1d_window_asymm,
        }
'''
