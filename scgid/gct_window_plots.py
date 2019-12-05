#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 21:05:52 2019

@author: kevinamses
"""

saved_points = points

from matplotlib import pyplot as plt
from scgid.flexwindow import calc_gc_window_1tailed, calc_gc_window_2tailed, calc_coverage_window_1tailed, calc_coverage_window_2tailed, calc_1d_window_asymm, calc_1d_window_symm, get_window_table, FlexibleSelectionWindow
from scgid.lib import generate_windows

from scgid.infotable import infotable
import numpy as np
import pandas as pd
import operator

target_taxa = {
        'target': ["Eukaryota"],
        'exception': []
        }

it = infotable(target_taxa)
it.load("/Users/kevinamses/Documents/Flux_Downlaods/dissertation/scgid/stylopage_scgid_output/blob/stylopage_info_table.tsv")

#windows = generate_windows(it, inc_factor=0.5)
#it.iter_descendants()

covMean = it.target_filter().summary_stats("coverage").mean
axis1W = calc_1d_window_asymm (it, "coverage", plot=False)
axis1T = it.rfilter("coverage", axis1W)
gcMean = axis1T.target_filter().summary_stats("gc").mean
axis2W = calc_1d_window_asymm (axis1T, "gc", plot=False)
print axis1W, axis2W

figureWindow = FlexibleSelectionWindow("co2gc2")
figureWindow.calculate(it, inc_factor=0.01, plot=True)
print figureWindow.show()

forR = (
        figureWindow.means["gc"],
        figureWindow.gc_range,
        figureWindow.means["coverage"],
        figureWindow.coverage_range
        )

for d in [-1,1]:
    figureWindow.gc_points[d].to_csv("/Users/kevinamses/Documents/Papers/scgid/figures/gct_window_calc/gc_{}_points.csv".format(d), sep=",", index=False)
    figureWindow.coverage_points[d].to_csv("/Users/kevinamses/Documents/Papers/scgid/figures/gct_window_calc/cov_{}_points.csv".format(d), sep=",", index=False)

testTable = infotable(target_taxa)
testTable.df = it.df
tabs = [testTable]
for i in range(0,20):
    now = tabs[i].spawn_child("gen{}".format(i), testTable.df)
    tabs.append(now)
testTable.iter_descendants()


    


print "\n".join([x.ident for x in it.children])
print [x.children for x in it.children]

print [(c.ident,c.depth) for c in it.children]

for c in it.children:
    for c2 in c.children:
        print c2.ident,c2.depth

dim1axis, dim1win, points = calc_1d_window_symm(it, "coverage", plot=True)

step1_it = it.rfilter(dim1axis, dim1win)
dim2axis, dim2win, gc_points = calc_1d_window_symm(step1_it, "gc", plot=True)

#calc_1d_window_symm(step1_it, "gc")

print dim1axis, dim1win
print dim2axis, dim2win

print ""

dim1axis, dim1win, points = calc_coverage_window_1tailed(it.df, target_taxa, plot=True)

step1_it = get_window_table(it.df, dim1win, "coverage")
dim2axis, dim2win, gc_points = calc_gc_window_1tailed(step1_it, target_taxa, plot=True)

print dim1axis, dim1win
print dim2axis, dim2win

testWindow = FlexibleSelectionWindow("co1gc1")
testWindow.calculate(it)
print testWindow.show()



### Output coverage points
pos = {
       'target': points['pos']['target'],
       'nontarget': points['pos']['nontarget'],
       'tradeoffs': points['pos']['tradeoffs']
       }
neg = {
       'target': points['neg']['target'],
       'nontarget': points['neg']['nontarget'],
       'tradeoffs': points['neg']['tradeoffs']
       }
pos_frame = pd.DataFrame(pos)
pos_frame['step'] = [ cov_stats.mean + (x * 0.01*cov_stats.std ) for x in range(1,len(points['pos']['target'])+1)]
pos_frame['direction'] = "positive"

neg_frame = pd.DataFrame(neg)
neg_frame['step'] = [cov_stats.mean-(x * 0.01*cov_stats.std) for x in range(1,len(points['neg']['target'])+1)]
neg_frame['direction'] = "negative"

pts_frame = pd.concat([pos_frame, neg_frame], axis=0, sort=False)
pts_frame.to_csv("/Users/kevinamses/Documents/Papers/scgid/figures/gct_window_calc/cov2_pts.csv", sep=",", index=False)

##output GC points
wt = get_window_table(it.df,(29.8953856255,1370.08335001),"coverage")
gc_mean = np.mean(wt.loc[wt.parse_lineage == "target"].gc)
gc_std = np.std(wt.loc[wt.parse_lineage == "target"].gc)
pos = {
       'target': gc_points['pos']['target'],
       'nontarget': gc_points['pos']['nontarget'],
       'tradeoffs': gc_points['pos']['tradeoffs']
       }
neg = {
       'target': gc_points['neg']['target'],
       'nontarget': gc_points['neg']['nontarget'],
       'tradeoffs': gc_points['neg']['tradeoffs']
       }
pos_frame = pd.DataFrame(pos)
pos_frame['step'] = [ gc_mean + (x * 0.01 * gc_std) for x in range(1,len(gc_points['pos']['target'])+1)]
pos_frame['direction'] = "positive"

neg_frame = pd.DataFrame(neg)
neg_frame['step'] = [gc_mean - (x * 0.01 * gc_std) for x in range(1,len(gc_points['neg']['target'])+1)]
neg_frame['direction'] = "negative"

pts_frame = pd.concat([pos_frame, neg_frame], axis=0, sort=False)
pts_frame.to_csv("/Users/kevinamses/Documents/Papers/scgid/figures/gct_window_calc/gc2_pts.csv", sep=",", index=False)

