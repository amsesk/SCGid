#!/usr/bin/env python
# coding: utf-8

import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
from scgid.infotable import InfoTable
from scgid.library import random_colors
from scgid.module import Module
from scgid.modcomm import get_head
from scgid.error import ModuleError
import sys
import os
import plotly

class PlotlyPlotter (Module):
    def __init__(self, infotable = None):
        self.head = get_head()

        self.taxon_colors = {}
        self.display_taxa = None
        
        self._set_data (infotable)
        self._scale()
    
    def _set_data (self, path):
        if self.head is None:
            if path is None:
                raise ValueError ("Data as SCGid.infotable.InfoTable must be provided to PlotlyPlotter when calling it from outside of SCGid runtime system.")

            else:
                self.path = path

        else:
            if path is None:
                self.path = self.head.path

            else:
                self.path = path

        self.infotable = InfoTable()
        self.infotable.load(path)

        return None


    def _scale (self):
        self.infotable.df.loc[self.infotable.df.evalue == 0, "evalue"] = 2.225074e-308
        self.infotable.df = self.infotable.df.assign(scaled_evalue = np.log(self.infotable.df.evalue))

    def collapse_top_n (self, top_n):

        choices = self.infotable.df.pertinent_taxlvl.value_counts()
        top_n_indices = choices.nlargest(top_n).index.to_list()

        self.infotable.df.loc[self.infotable.df.pertinent_taxlvl.isin(top_n_indices) == False, "pertinent_taxlvl"] = "Other"
        self.display_taxa = self.infotable.df.pertinent_taxlvl.value_counts().index.to_list()

    def set_colors(self, palette = None):
        if palette is None:
            palette = [f"rgb({','.join([str(v) for v in x])})" for x in random_colors(len(self.display_taxa), minimum_distance=20)]

        self.taxon_colors = dict(zip(self.display_taxa, palette))
        return None


if len(sys.argv) > 1:
    it_path = sys.argv[1]

else:
    it_path = "/home/aimzez/work/allomyces/Burma1F_EDFdb_scgid_output/gct/Burma1F_EDFdb.infotable.tsv"

prefix = os.path.split(it_path)[1].split(".")[0]

plotter = PlotlyPlotter(infotable = it_path)

plotter.collapse_top_n(10)
plotter.set_colors()


tax_to_color=plotter.taxon_colors

fig = make_subplots(rows=2, cols=2, column_widths=[0.8,0.2], row_heights=[0.2,0.8], 
                    shared_xaxes=True, shared_yaxes=True, vertical_spacing=0.01, horizontal_spacing=0.01)
for lin in plotter.display_taxa:
    
    # Get values for this trace
    subframe = plotter.infotable.df[plotter.infotable.df.pertinent_taxlvl == lin]
    gc_counts,gc_div = np.histogram(plotter.infotable.df[plotter.infotable.df.pertinent_taxlvl == lin].gc, bins=10000)
    cov_counts,cov_div = np.histogram(plotter.infotable.df[plotter.infotable.df.pertinent_taxlvl == lin].coverage, bins=10000)
    
    # Add dots trace for each individual category
    fig.add_trace(
        go.Scatter(
            x=np.log(subframe.coverage), 
            y=subframe.gc,
            name=lin,
            mode="markers",
            marker=dict(
                color=tax_to_color[lin], 
                size=np.sqrt(abs(subframe.scaled_evalue))*1.5, 
                opacity=0.4),
            hovertext = subframe.contig
        ),
        row=2,
        col=1
    )
    
    # Add coverage histogram for category on top of main panel
    fig.add_trace( 
        go.Scatter(
            x=np.log(cov_div), 
            y=np.sqrt(cov_counts), 
            line=dict(color=tax_to_color[lin]),
            showlegend=False
        ), 
        row=1, 
        col=1,
    )
    
    # Add GC histogram for category to right side of main panel
    fig.add_trace( 
        go.Scatter(
            x=gc_counts, 
            y=gc_div, 
            line=dict(color=tax_to_color[lin]),
            showlegend=False
        ), 
        row=2, 
        col=2,
    )
    
fig.update_layout(legend = {"itemsizing":"constant"})

fig.update_xaxes(title="log(Coverage)", row=2, col=1)
fig.update_yaxes(title="GC Content", row=2, col=1)

fig.update_xaxes(row=2, col=2, showticklabels = False)
fig.update_yaxes(row=1, col=1, showticklabels = False)

#fig.show()
plotly.offline.plot(fig, filename=f"{prefix}.gctplt.html", auto_open=False)
