
import igraph as ig
import networkx as nx
import seaborn as sb
import plotly as ply
import altair as alt
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import Network_Metrics as NM
from pathlib import Path
from PIL import Image
import io as io
import cairocffi as cfi
from matplotlib import pyplot as plt
from pathlib import Path
import re as re
import types as types
from typing import Union
import plotly.graph_objs as go
import os as os
ig.config["plotting.backend"] = "matplotlib"
# ig.config["plotting.backend"] = "cairo"

class SetDirectories:
    def __init__(self):
        self.root_dir = Path("/Users/GBS/Master/Figures/iGraph")

    def get_output_dir(self):
        return self.root_dir


def mcf7_10_raw_lowres_graphs():
    root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/lowres_mcf7_mcf10/raw")
    graph_creator = NM.CreateGraphsFromDirectory(root_dir)
    graph_creator.from_edgelists()
    mcf7_10_graphs = graph_creator.graph_dict
    return mcf7_10_graphs

def mcf7_10_norm_lowres_graphs():
    root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/lowres_mcf7_mcf10/norm")
    graph_creator = NM.CreateGraphsFromDirectory(root_dir)
    graph_creator.from_edgelists()
    mcf7_10_graphs = graph_creator.graph_dict
    return mcf7_10_graphs


def mcf7_chr18_1mb():
    graph_filter = NM.FilterGraphs(mcf7_10_norm_lowres_graphs())
    filtered_graphs = graph_filter.filter_graphs(cell_lines=["mcf10"], chromosomes=["chr18"], resolutions=["500000"])
    # graph_filter.print_filtered_edges(filtered_graphs)
    return filtered_graphs
mcf7_chr18_1mb()

# TODO: Make plotting class that takes any graph dict from any class in Network_Metrics and plots it as a network
#   need to compare the LCC to the full graphs, because the number of communities and merges are the same but the sizes are different.


class plot_graph:

    def __init__(self, graph_dict, output_dir):
        self.graph_dict = graph_dict
        self.output_dir = output_dir

    @staticmethod
    def abbreviate_label(self, label, resolution):
        chrom, region = label.split(':')
        start, end = region.split('-')

        if resolution >= 1000000:
            start_unit = int(start) // 1000000
            end_unit = int(end) // 1000000
            unit = "MB"
        elif resolution >= 1000:
            start_unit = int(start) // 1000
            end_unit = int(end) // 1000
            unit = "KB"
        else:
            start_unit = int(start)
            end_unit = int(end)
            unit = "B"

        # Use to label graphs like this:
        # resolution = int(re.search(r"_([^_]+)_", graph_name).group(1))
        # abbreviated_labels = [self.abbreviate_label(label, resolution) for label in graph.vs["name"]]
        return f"{chrom}:{start_unit}-{end_unit} {unit}"

    def show_graph(self):
        for graph_name, graph in self.graph_dict.items():
            ig.plot(graph)
            plt.show()

    def show_graph_with_lcc(self, largest_component_obj):
        for graph_name, graph in self.graph_dict.items():
            lcc_membership = largest_component_obj.lcc_membership()[graph_name]

            node_colors = ["red" if membership == 1 else "black" for membership in lcc_membership]
            edge_colors = ["red" if lcc_membership[edge.source] == 1 and lcc_membership[edge.target] else "black" for edge in graph.es]
            ig.plot(graph, vertex_color=node_colors, edge_color=edge_colors)
            plt.show()

    def save_graph(self):
        for graph_name, graph in self.graph_dict.items():
            visual_style = {"vertex_size": 1, "vertex_label": graph.vs["name"], "layout": graph.layout("kk")}
            # Remove special characters from graph name to use as file name:
            safe_graph_name = "".join(e for e in graph_name if e.isalnum() or e == "_")
            output_filename = self.output_dir / f"{safe_graph_name}.png"
            ig.plot(graph, output_filename, **visual_style)
            print(f"Saved plot to {output_filename}")

def plot_full():
    dir_manager = SetDirectories()
    output_dir = dir_manager.get_output_dir()
    plot = plot_graph(mcf7_chr18_1mb(), output_dir)
    plot.show_graph()
# plot_full()

def plot_lcc():
    dir_manager = SetDirectories()
    output_dir = dir_manager.get_output_dir()
    graph_dict = mcf7_chr18_1mb()
    plot = plot_graph(graph_dict, output_dir)
    largest_component_obj = NM.LargestComponent(graph_dict)
    plot.show_graph_with_lcc(largest_component_obj)
plot_lcc()


# TODO: Make stacked bar plots showing ratio of nodes in LCC to total nodes for each cell line and each resolution, takes graph dicts as input.

class plot_lcc_ratio:
    def __init__(self, graph_dict, output_dir):
        self.graph_dict = graph_dict
        self.output_dir = output_dir

    def get_lcc_ratio(self):
        for graph_name, graph in self.graph_dict.items():
            lcc = NM.LCC_ratio(graph).calculate_lcc_ratio()

    def plot_lcc_ratio_bar(self):
        for graph_name, graph in self.graph_dict.items():
            lcc = graph.components(mode="weak").giant()
            lcc_ratio = len(lcc.vs) / len(graph.vs)
            print(f"{graph_name}: {lcc_ratio}")



# TODO: Make degree distribution plots for each cell lines LCC for each chromosome, per resolution? Takes graph dict as input.

# TODO: Make betweenness plot that colors the nodes by betwenness centrality, takes a graph dict as input. Make

# TODO: Make closeness plot that colors the nodes by closeness, takes graph dict as input.



# Networkx plot

class PlotNetworkxGraphs:

    def __init__(self, networkx_graph_dict):
        self.networkx_graph_dict = networkx_graph_dict

    def plot(self):
        for graph_key, nx_graph in self.networkx_graph_dict.items():
            plt.figure()
            plt.title(f"{graph_key}")
            nx.draw(nx_graph, with_labels=True, node_color="skyblue", font_weight="bold", node_size=1000)
            plt.show()

# plotter = PlotNetworkxGraphs(ConvertIgraphToNetworkx(chr18_inc_norm_graphs_50kb()).convert())
# plotter.plot()








