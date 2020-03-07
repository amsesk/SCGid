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
from scgid.error import ModuleError, Ok, check_result
import sys
import os
import plotly

def change_colors(plotter, level, extra_data):
    #print (len(plotter.infotable.df.lineage))
    #plotter.infotable.df.pertinent_taxlvl = [len(str(x).replace("[","").replace("]","").split(",")) for x in plotter.infotable.df.lineage]
    new_perts = []
    for i in plotter.infotable.df.lineage:
        try:
            new_perts.append(str(i[level]))
        except IndexError:
            new_perts.append("Unclassified")

    extra_data = extra_data.replace({"pertinent_taxlvl": new_perts})
    print(extra_data)
    #print(new_perts)
    #print(plotter.infotable.df.pertinent_taxlvl)
    plotter.collapse_top_n(10)
    plotter.set_colors()
    return [{'marker.color': [plotter.infotable.df.color.to_list()], 'customdata': extra_data.to_numpy()}]

class PlotlyPlotter (Module):
    def __init__(self, infotable = None, n = 10):
        self.head = get_head()

        self.infotable = None
        self.n = n

        self.taxon_colors = {}
        self.display_taxa = None
        
        self._set_data (infotable)
        self._scale()
    
    def _set_data (self, infotable):
        # Get infotable from module head or function kwarg or raise ValueError if both are None
        if self.head is None and infotable is None:
            raise ValueError ("Data as SCGid.infotable.InfoTable must be provided to PlotlyPlotter when calling it from outside of SCGid runtime system.")
        elif self.head is None and self.infotable is not None:
            self.infotable = infotable
        else:
            self.infotable = self.head.infotable

        ##############################################################################
        # If infotable is not and instance of scgid.infotable.Infotable              #
        # (i.e., from a head module), it should be a path to one that can be loaded. #
        ##############################################################################=
        if not isinstance(infotable, InfoTable):
            self.infotable = InfoTable()
            self.infotable.load(path)

        self.infotable.df["log_coverage"] = np.log(self.infotable.df.coverage)
        self.infotable.df.evalue = self.infotable.df.evalue.astype("float64")

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

        self.infotable.df["color"] = [self.taxon_colors[x] for x in self.infotable.df["pertinent_taxlvl"]]
        #print(self.infotable.df["color"])
        return None

    @check_result
    def plot(self, outpath):
        self.collapse_top_n(self.n)
        self.set_colors()

        fig = make_subplots(
            rows=2, 
            cols=2, 
            column_widths=[0.8,0.2], 
            row_heights=[0.2,0.8], 
            shared_xaxes=True, 
            shared_yaxes=True, 
            vertical_spacing=0.01, 
            horizontal_spacing=0.01
            )

        for lin in self.display_taxa:
    
            # Get values for this trace
            subframe= self.infotable.df[self.infotable.df.pertinent_taxlvl == lin]
            gc_counts,gc_div = np.histogram(self.infotable.df[self.infotable.df.pertinent_taxlvl == lin].gc, bins=1000)
            cov_counts,cov_div = np.histogram(self.infotable.df[self.infotable.df.pertinent_taxlvl == lin].log_coverage, bins=1000)
            extra_data = subframe[["contig", "coverage", "evalue", "pertinent_taxlvl", "lineage", "pid"]]
            extra_data.loc[:,"lineage"] = extra_data["lineage"].apply(lambda x: "<br>".join(x))
         
            points = go.Scatter(
                        x=subframe.log_coverage, 
                        y=subframe.gc,
                        name=lin,
                        mode="markers",
                        marker=dict(
                            color=self.taxon_colors[lin], 
                            size=np.sqrt(abs(subframe.scaled_evalue))*1.5, 
                            opacity=0.4),
                        customdata = extra_data,
                        hovertemplate = "<b>%{customdata[0]}; %{customdata[5]}</b><br>GC Content: %{y:.2f}<br>Coverage: %{customdata[1]:.2f}<br>ln(Coverage): %{x:.2f}<br>Taxonomy: %{customdata[3]}<br>e-value: %{customdata[2]:.2e}<extra><b>Lineage:</b><br>%{customdata[4]}</extra>"
                    )
            # Add dots trace for each individual category
            fig.add_trace(
                    points,
                    row=2,
                    col=1
                )

            # Add coverage histogram for category on top of main panel
            fig.add_trace( 
                go.Scatter(
                    x=cov_div, 
                    y=np.sqrt(cov_counts), 
                    name=lin,
                    line=dict(color=self.taxon_colors[lin]),
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
                    name=lin,
                    line=dict(color=self.taxon_colors[lin]),
                    showlegend=False
                ), 
                row=2, 
                col=2,
            )
         
        fig.update_layout(legend = {"itemsizing":"constant"})

        fig.update_xaxes(title="ln(Coverage)", row=2, col=1)
        fig.update_yaxes(title="GC Content", row=2, col=1)

        fig.update_xaxes(row=2, col=2, showticklabels = False)
        fig.update_yaxes(row=1, col=1, showticklabels = False)

        plotly.offline.plot(fig, filename=outpath, auto_open=False)

        return Ok()



'''
if len(sys.argv) > 1:
    it_path = sys.argv[1]

extra_data = plotter.infotable.df[["contig", "coverage", "evalue", "pertinent_taxlvl"]]
points = go.Scatter(
            x=plotter.infotable.df.log_coverage, 
            y=plotter.infotable.df.gc,
            name="All",
            mode="markers",
            marker=dict(
                color=plotter.infotable.df.color, 
                size=np.sqrt(abs(plotter.infotable.df.scaled_evalue))*1.5, 
                opacity=0.4
                ),
            customdata = extra_data.to_numpy(),
            hovertemplate = "Contig Name:  %{customdata[0]}<br>GC Content: %{y:.2f}<br>Coverage: %{customdata[1]:.2f}<br>ln(Coverage): %{x:.2f}<br>Taxonomy: %{customdata[3]}<br>e-value: %{customdata[2]:.2e}"
            )
    # Add dots trace for each individual category
fig.add_trace(
        points,
        row=2,
        col=1
    )

print(dir(points))


for lin in self.display_taxa:
    
    # Get values for this trace
    subframe= self.infotable.df[self.infotable.df.pertinent_taxlvl == lin]
    gc_counts,gc_div = np.histogram(self.infotable.df[self.infotable.df.pertinent_taxlvl == lin].gc, bins=1000)
    cov_counts,cov_div = np.histogram(self.infotable.df[self.infotable.df.pertinent_taxlvl == lin].log_coverage, bins=1000)
    extra_data = subframe[["contig", "coverage", "evalue", "pertinent_taxlvl", "lineage", "pid"]]
    extra_data.loc[:,"lineage"] = extra_data["lineage"].apply(lambda x: "<br>".join(x))
 
    points = go.Scatter(
                x=subframe.log_coverage, 
                y=subframe.gc,
                name=lin,
                mode="markers",
                marker=dict(
                    color=tax_to_color[lin], 
                    size=np.sqrt(abs(subframe.scaled_evalue))*1.5, 
                    opacity=0.4),
                customdata = extra_data,
                hovertemplate = "<b>%{customdata[0]}; %{customdata[5]}</b><br>GC Content: %{y:.2f}<br>Coverage: %{customdata[1]:.2f}<br>ln(Coverage): %{x:.2f}<br>Taxonomy: %{customdata[3]}<br>e-value: %{customdata[2]:.2e}<extra><b>Lineage:</b><br>%{customdata[4]}</extra>"
            )
    # Add dots trace for each individual category
    fig.add_trace(
            points,
            row=2,
            col=1
        )
    #print (points)

    # Add coverage histogram for category on top of main panel
    fig.add_trace( 
        go.Scatter(
            x=cov_div, 
            y=np.sqrt(cov_counts), 
            name=lin,
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
            name=lin,
            line=dict(color=tax_to_color[lin]),
            showlegend=False
        ), 
        row=2, 
        col=2,
    )

    fig.update_layout(
    updatemenus=[
        dict(
            buttons=list([
                dict(
                    args=change_colors(plotter, 1, extra_data),
                    label="Kingdom",
                    method="restyle"
                ),
                dict(
                    args=change_colors(plotter, 3, extra_data),
                    label="Phylum",
                    method="restyle"
                )
            ]),
            direction="down",
            pad={"r": 10, "t": 10},
            showactive=True,
            x=1.0,
            xanchor="left",
            y=1.1,
            yanchor="top"
        ),
    ]
)

fig.update_layout(legend = {"itemsizing":"constant"})

fig.update_xaxes(title="ln(Coverage)", row=2, col=1)
fig.update_yaxes(title="GC Content", row=2, col=1)

fig.update_xaxes(row=2, col=2, showticklabels = False)
fig.update_yaxes(row=1, col=1, showticklabels = False)
print(dir(fig))

#fig.show()
plotly.offline.plot(fig, filename=f"{prefix}.gctplt.html", auto_open=True)

'''