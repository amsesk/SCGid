#!/usr/bin/env python
# coding: utf-8

import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
from scgid.infotable import InfoTable
from scgid.library import random_colors
import sys
import plotly

if len(sys.argv) > 1:
    it_path = sys.argv[1]
else:
    it_path = "/home/aimzez/work/allomyces/Burma1F_EDFdb_scgid_output/gct/Burma1F_EDFdb.infotable.tsv"

prefix = os.path.split(it_path).split(".")[0]

it = InfoTable()
it.load(it_path)
it.df = it.df.assign(scaled_evalue = np.log(it.df.evalue))

tnt_to_color = {
    "target": "blue",
    "nontarget": "red"
}

tax_to_color = {}
choices = list(set(it.df.pertinent_taxlvl))
colpal = [x for x in random_colors(len(choices), minimum_distance=1)]
tax_to_color = {choices[i]: f"rgb({','.join([str(c) for c in colpal[i]])})" for i in range(0,len(choices))}

fig = make_subplots(rows=2, cols=2, column_widths=[0.8,0.2], row_heights=[0.2,0.8], 
                    shared_xaxes=True, shared_yaxes=True, vertical_spacing=0.01, horizontal_spacing=0.01)
for lin in choices:
    
    # Get values for this trace
    subframe = it.df[it.df.pertinent_taxlvl == lin]
    gc_counts,gc_div = np.histogram(it.df[it.df.pertinent_taxlvl == lin].gc, bins=10000)
    cov_counts,cov_div = np.histogram(it.df[it.df.pertinent_taxlvl == lin].coverage, bins=10000)
    
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
                opacity=0.4)
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
    
fig.update_layout(showlegend=True)

fig.update_xaxes(title="log(Coverage)", row=2, col=1)
fig.update_yaxes(title="GC Content", row=2, col=1)

fig.update_xaxes(row=2, col=2, showticklabels = False)
fig.update_yaxes(row=1, col=1, showticklabels = False)

#fig.show()
plotly.offline.plot(fig, filename=f"{prefix}.gctplt.html", autoopen=False)
