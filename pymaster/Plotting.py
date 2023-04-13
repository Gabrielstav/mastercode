
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

def mcf7_10_lowres_graphs():
    root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/lowres_mcf7_mcf10")
    graph_creator = NM.CreateGraphsFromDirectory(root_dir)
    graph_creator.from_edgelists()
    mcf7_10_graphs = graph_creator.graph_dict
    return mcf7_10_graphs

def mcf7_chr18_1mb():
    graph_filter = NM.FilterGraphs(mcf7_10_lowres_graphs())
    filtered_graphs = graph_filter.filter_graphs(cell_lines=["mcf7"], chromosomes=["chr18", "chrX"], resolutions=["1000000"])
    # graph_filter.print_filtered_edges(filtered_graphs)
    return filtered_graphs
mcf7_chr18_1mb()

# TODO: Make plotting class that takes any graph dict from any class in Network_Metrics and plots it as a network
#   need to compare the LCC to the full graphs, because the number of communities and merges are the same but the sizes are different.

class plot_graph:
    def __init__(self, graph_dict):
        self.graph_dict = graph_dict

    def plot_graph(self):
        for graph in self.graph_dict:
            ig.plot(graph)
            plt.show()

plot_graph(mcf7_chr18_1mb())


# TODO: Make LCC plot, where the network is plotted with the LCC highlighted, takes graph dict as input.

# TODO: Make stacked bar plots showing ratio of nodes in LCC to total nodes for each cell line and each resolution, takes graph dicts as input.

# TODO: Make degree distribution plots for each cell lines LCC for each chromosome, per resolution? Takes graph dict as input.

# TODO: Make betweenness plot that colors the nodes by betwenness centrality, takes a graph dict as input. Make

# TODO: Make closeness plot that colors the nodes by closeness, takes graph dict as input.





