# Import modules
import math
import networkx as nx
import igraph as ig
from matplotlib import pyplot as plt
from pathlib import Path
import re as re
ig.config["plotting.backend"] = "matplotlib"


import numpy as np
# from igraph import plot as iplot
from collections import Counter
import altair as alt
import pandas as pd
import seaborn as sb
import os as os
import scipy.stats as stats
import pandas as pd


# import pycairo as pc


# TODO: Eventually make class to read in edgelists, pass those to another class that makes graph objects, and then pass those to a metric class that calculates the metrics we want to plot.
#       This class then passes the metrics to a plotting class that plots the metrics. Which can be saved to an output directory in the end.
#       Also need to include overlapping nodes into networks, since we do not want to compare nodes not overlapping, but for this we need filtering classes first.

# TODO: Read in data, just set as hardcoded paths to begin with (later on we can make a class to read in data)

class Directories:
    """
    Set input directories containing edge lists here,
    and output directories to save figures to.
    """

    root_dir = Path("/Users/GBS/Master/HiC-Data/Pipeline_out")
    norm_dir = Path("/Users/GBS/Master/HiC-Data/Pipeline_out/chr18_INC/chr18_norm/edgelists")
    raw_dir = Path("/Users/GBS/Master/HiC-Data/Pipeline_out/chr18_INC/chr18_raw/edgelists")
    fig_dir = Path("/Users/GBS/Master/Figures/Network_metrics/chr18_inc")

    @classmethod
    def root_directory(cls, root_dir=None):
        if root_dir is not None:
            cls.root_dir = root_dir
        return cls.root_dir

    @classmethod
    def normalized_directory(cls, norm_dir=None):
        if norm_dir is not None:
            cls.norm_dir = norm_dir
        return cls.norm_dir

    @classmethod
    def raw_directory(cls, raw_dir=None):
        if raw_dir is not None:
            cls.raw_dir = raw_dir
        return cls.raw_dir

    @classmethod
    def figure_directory(cls, fig_dir=None):
        if fig_dir is not None:
            cls.fig_dir = fig_dir
        return cls.fig_dir


class Create_graphs:

    @staticmethod
    def get_files():
        edgelist_files = (file_path for file_path in Directories.root_directory().rglob('**/*edgelist*.txt') if file_path.is_file() and "weighted" not in file_path.name)
        edgelist_files_list = [str(file_path.resolve()) for file_path in edgelist_files]
        for file_path in edgelist_files_list:
            print(f"Processing file: {file_path}")
        return edgelist_files_list

    @staticmethod
    def get_weighted_files():
        weighted_files = (file_path for file_path in Directories.root_directory().rglob('**/*weighted*.txt') if file_path.is_file())
        weighted_files_list = [str(file_path.resolve()) for file_path in weighted_files]
        for file_path in weighted_files_list:
            print(f"Processing file: {file_path}")
        return weighted_files_list

    @staticmethod
    def from_edgelists():
        graph_objects = {}
        for edgelist in Create_graphs.get_files():
            # Extract graph name from file path
            pattern = r".*/(.+)_edgelist\.txt$"
            match = re.match(pattern, str(edgelist))
            graph_name = match.group(1)
            graph_objects[graph_name] = ig.Graph.Load(edgelist, format="ncol")

        return graph_objects

    @staticmethod
    def from_weighted_edgelists():
        graph_objects = {}
        for edgelist in Create_graphs.get_weighted_files():
            # Extract graph name from file path
            pattern = r".*/(.+)_weighted_edgelist\.txt$"
            match = re.match(pattern, str(edgelist))
            graph_name = match.group(1)
            graph_objects[graph_name] = ig.Graph.Load(edgelist, format="ncol")

        return graph_objects


# TODO: This is messy RN; automate so each method takes the previous method as input, and then the last method returns the final output (plots).
#    calculate largest connected ocmponent, and then calculate metrics for that component only. Compare with metrics for whole network?
#



# Need to make the graph objects for each file read in, and then pass them to a class





# {'chr18_lowres_500000': <igraph.Graph object at 0x12122e340>, ...}

# g= ig.Graph.Load("/Users/GBS/Master/HiC-Data/Pipeline_out/chr18_INC/chr18_raw/edgelists/chr18_lowres_50000_edgelist.txt", format="ncol")
# ig.plot(g)
# # ig.plot(g, vertex_size=2, vertex_label=g.vs["name"], vertex_label_size=1, vertex_label_dist=1.5, vertex_label_color="black", layout=g.layout("kk")) # TODO: Figure out iptimal plot layout
# plt.show()

def plot_20kb_hires_raw():
    h = ig.Graph.Load("/Users/GBS/Master/HiC-Data/Pipeline_out/chr18_INC/chr18_raw/edgelists/chr18_hires_20000_edgelist.txt", format="ncol", directed=False) # very messy
    ig.plot(h, edge_width=0.07, node_size=0.5, node_color="red")
    plt.show()
    plt.close()

def plot_20kb_hires_norm():
    g = ig.Graph.Load("/Users/GBS/Master/HiC-Data/Pipeline_out/chr18_INC/chr18_norm/edgelists/chr18_hires_iced_20000_edgelist.txt", format="ncol", directed=False) # every node connected
    ig.plot(g, edge_width=0.07, node_size=0.5, node_color="red")
    plt.show()
    plt.close()

def plot_500kb_raw():
    g = ig.Graph.Load("/Users/GBS/Master/HiC-Data/Pipeline_out/chr18_INC/chr18_raw/edgelists/chr18_lowres_500000_edgelist.txt", format="ncol", directed=False)
    ig.plot(g, edge_width=0.07, node_size=0.5, node_color="red")
    plt.show()
    plt.close()

def plot_500kb_norm():
    g = ig.Graph.Load("/Users/GBS/Master/HiC-Data/Pipeline_out/chr18_INC/chr18_norm/edgelists/chr18_lowres_iced_500000_edgelist.txt", format="ncol", directed=False)
    ig.plot(g, edge_width=0.07, node_size=0.5, node_color="red")
    plt.show()
    plt.close()

def plot_50kb_raw():
    g = ig.Graph.Load("/Users/GBS/Master/HiC-Data/Pipeline_out/chr18_INC/chr18_raw/edgelists/chr18_lowres_50000_edgelist.txt", format="ncol", directed=False)
    ig.plot(g, edge_width=0.1, node_size=0.9, node_color="red")
    plt.show()
    plt.close()

def plot_50kb_norm():
    g = ig.Graph.Load("/Users/GBS/Master/HiC-Data/Pipeline_out/chr18_INC/chr18_norm/edgelists/chr18_lowres_iced_50000_edgelist.txt", format="ncol", directed=False)
    ig.plot(g, edge_width=0.1, node_size=0.9, node_color="red")
    plt.show()
    plt.close()

def plot_120kb_raw():
    g = ig.Graph.Load("/Users/GBS/Master/HiC-Data/Pipeline_out/chr18_INC/chr18_raw/edgelists/chr18_hires_120000_edgelist.txt", format="ncol", directed=False)
    ig.plot(g, edge_width=0.1, node_size=0.9, node_color="red")
    plt.show()
    plt.close()

def plot_120kb_norm():
    g = ig.Graph.Load("/Users/GBS/Master/HiC-Data/Pipeline_out/chr18_INC/chr18_norm/edgelists/chr18_hires_iced_120000_edgelist.txt", format="ncol", directed=False)
    ig.plot(g, edge_width=0.1, node_size=0.9, node_color="red")
    plt.show()
    plt.close()

def plot_250kb_raw():
    g = ig.Graph.Load("/Users/GBS/Master/HiC-Data/Pipeline_out/chr18_INC/chr18_raw/edgelists/chr18_lowres_250000_edgelist.txt", format="ncol", directed=False)
    ig.plot(g, edge_width=0.1, node_size=0.9, node_color="red")
    plt.show()
    plt.close()

def plot_250kb_norm():
    g = ig.Graph.Load("/Users/GBS/Master/HiC-Data/Pipeline_out/chr18_INC/chr18_norm/edgelists/chr18_lowres_iced_250000_edgelist.txt", format="ncol", directed=False)
    ig.plot(g, edge_width=0.1, node_size=0.9, node_color="red")
    plt.show()
    plt.close()

plot_50kb_raw()

kb50_raw = ig.Graph.Load("/Users/GBS/Master/HiC-Data/Pipeline_out/chr18_INC/chr18_raw/edgelists/chr18_lowres_50000_edgelist.txt", format="ncol", directed=False)
kb50_norm = ig.Graph.Load("/Users/GBS/Master/HiC-Data/Pipeline_out/chr18_INC/chr18_norm/edgelists/chr18_lowres_iced_50000_edgelist.txt", format="ncol", directed=False)

def fg_kb50_raww():
    fg_kb50_raw = kb50_raw.community_fastgreedy()
    communities_kb50_raw = fg_kb50_raw.as_clustering()
    print(communities_kb50_raw.modularity)

def fg_kb50_normm():
    fg_kb50_norm = kb50_norm.community_fastgreedy()
    communities_kb50_norm = fg_kb50_norm.as_clustering()
    print(communities_kb50_norm.modularity)

# fg_kb50_normm()

def fg_kb500_raww():
    g = ig.Graph.Load("/Users/GBS/Master/HiC-Data/Pipeline_out/chr18_INC/chr18_norm/edgelists/chr18_lowres_iced_500000_edgelist.txt", format="ncol", directed=False)
    fg_kb50_raw = g.community_fastgreedy()
    communities_kb50_raw = fg_kb50_raw.as_clustering()
    print(communities_kb50_raw.modularity)

def fg_kb500_normm():
    g = ig.Graph.Load("/Users/GBS/Master/HiC-Data/Pipeline_out/chr18_INC/chr18_norm/edgelists/chr18_lowres_iced_500000_edgelist.txt", format="ncol", directed=False)
    fg_kb50_norm = g.community_fastgreedy()
    communities_kb50_norm = fg_kb50_norm.as_clustering()
    print(communities_kb50_norm.modularity)

def fg_kb20_raww():
    h = ig.Graph.Load("/Users/GBS/Master/HiC-Data/Pipeline_out/chr18_INC/chr18_raw/edgelists/chr18_hires_20000_edgelist.txt", format="ncol", directed=False) # very messy
    fg_kb50_raw = h.community_fastgreedy()
    communities_kb50_raw = fg_kb50_raw.as_clustering()
    print(communities_kb50_raw.modularity)

def fg_kb20_normm():
    g = ig.Graph.Load("/Users/GBS/Master/HiC-Data/Pipeline_out/chr18_INC/chr18_norm/edgelists/chr18_hires_iced_20000_edgelist.txt", format="ncol", directed=False) # every node connected
    fg_kb50_norm = g.community_fastgreedy()
    communities_kb50_norm = fg_kb50_norm.as_clustering()
    print(communities_kb50_norm.modularity)

def fg_kb120_raww():
    g = ig.Graph.Load("/Users/GBS/Master/HiC-Data/Pipeline_out/chr18_INC/chr18_raw/edgelists/chr18_hires_120000_edgelist.txt", format="ncol", directed=False)
    fg_kb50_raw = g.community_fastgreedy()
    communities_kb50_raw = fg_kb50_raw.as_clustering()
    print(communities_kb50_raw.modularity)

def fg_kb120_normm():
    g = ig.Graph.Load("/Users/GBS/Master/HiC-Data/Pipeline_out/chr18_INC/chr18_norm/edgelists/chr18_hires_iced_120000_edgelist.txt", format="ncol", directed=False)
    fg_kb50_norm = g.community_fastgreedy()
    communities_kb50_norm = fg_kb50_norm.as_clustering()
    print(communities_kb50_norm.modularity)

def fg_kb250_raww():
    g = ig.Graph.Load("/Users/GBS/Master/HiC-Data/Pipeline_out/chr18_INC/chr18_raw/edgelists/chr18_lowres_250000_edgelist.txt", format="ncol", directed=False)
    fg_kb50_raw = g.community_fastgreedy()
    communities_kb50_raw = fg_kb50_raw.as_clustering()
    print(communities_kb50_raw.modularity)

def fg_kb250_normm():
    g = ig.Graph.Load("/Users/GBS/Master/HiC-Data/Pipeline_out/chr18_INC/chr18_norm/edgelists/chr18_lowres_iced_250000_edgelist.txt", format="ncol", directed=False)
    fg_kb50_norm = g.community_fastgreedy()
    communities_kb50_norm = fg_kb50_norm.as_clustering()
    print(communities_kb50_norm.modularity)

fg_kb250_raww()
fg_kb250_normm()

# kb50_norm.isomorphic(kb50_raw)
#
# degree_distribution_raw_50kb = kb50_raw.degree_distribution()
# print(degree_distribution_raw_50kb)
# degree_distribution_norm_50kb = kb50_norm.degree_distribution()
#
# plt.scatter(degree_distribution_raw_50kb, degree_distribution_raw_50kb, color="red")
# plt.show()


# dd = g.degree_distribution()
# degree_distribution = dd.bins()
# ig.plot(dd)
# plt.show()
# plt.close()

# How to make undirected?
# c = nx.read_edgelist("/Users/GBS/Master/HiC-Data/Pipeline_out/chr18_INC/chr18_raw/edgelists/chr18_lowres_50000_edgelist.txt", directed=False, edgetype=str, nodetype=str, data=(('weight', float),), create_using=nx.DiGraph())
# nx.draw(c, node_size=10)
# plt.show()


# ex = ig.Graph.Famous("petersen")
# ig.plot(ex)
#
# plt.plot(ex)
#
# y = nx.complete_graph(5)
# nx.draw(y)


# class Degree_distribution:
#
#     @staticmethod
#     def get_degree(graph_objects):
#         degree = {}
#         for graph_name, graph in graph_objects.items():
#             degree[graph_name] = graph.degree()
#         return degree

# print(Degree_distribution.get_degree(Create_graphs.from_edgelists()))

# Ugly af method
# @staticmethod
# def plot_degree_distribution():
#     graphs = Create_graphs.from_edgelists()
#     for graph_name, graph in graphs.items():
#         degrees = graph.degree()
#         plt.figure()
#         plt.scatter(degrees, range(len(degrees)), s=5, color='red')
#         plt.xscale('log')
#         plt.xlim(1, max(degrees))
#         plt.ylim(0, len(degrees))
#         plt.title(graph_name)
#         plt.xlabel('Degree')
#         plt.ylabel('Frequency')
#         plot_path = Directories.figure_directory()
#         plot_path /= f"{graph_name}_degree_distribution.png"
#         print(plot_path)
#         plt.savefig(str(plot_path))

# First test:
# @staticmethod
# def plot_degree_distribution(degree_distributions):
#     for graph_name, degrees in degree_distributions.items():
#         freqs = [degrees.count(d) for d in degrees]
#         plt.scatter(degrees, freqs)
#         plt.xscale('log')
#         plt.yscale('-log')
#         plt.xlabel('Degree')
#         plt.ylabel('Frequency')
#         plt.title(f'Degree distribution for {graph_name}')
#         plt.savefig(str(Directories.figure_directory() / f"{graph_name}_degree_distribution.png"))
#         plt.close()


#
# graphs = Create_graphs.from_edgelists()
# degree_distributions = Degree_distribution.get_degree(graphs)
# Degree_distribution.plot_degree_distribution(degree_distributions)


# TODO: Make degree to map to edgelists, so we can plot the degree distribution for each bin. Make class for this?
# @staticmethod
# def map_degree_to_edgelists(degree, edgelists):
#     mapped_degree = {}
#     for edgelists, degree in zip(edgelists, degree):
#         mapped_degree[edgelists] = degree
#     return mapped_degree


# TODO: Create graph objects from data per resolution
# Maybe make function to split data into different resolutions per root dir provided? Quick fix.

# TODO: Figure out what metrics to calculate, and what to plot.
# think about what we want to show in the end.
# Maybe we can make a function that takes in a graph object and calculates the metrics we want to plot? IDk about speed tho, but definiatey for commong stuff like degree and betweenness etc.

# TODO: Then do plotting quick for now, that is, just plot the metrics we want to show without making the code automatic and flexible.


# TODO: Do community detection for all graphs, and plot the communities for each graph.

# TODO: Jaccard index for all graphs, and plot the Jaccard index for each graph.

# TODO: Calculate modulariry for each graph

# TODO: Calcualte


# Quick look at new method on INC data


# Check if we can create graph objects directly without having to create an edgelist format file each time and
# import it from a different directory.

# K562 chromosome 18 (just testing iGraph stuff)
# K562_chr18 = ig.Graph.Load("/Users/GBS/Master/HiC-Data/Processed_Data_Edgelist/HiC_from_Jonas/Chr18/K562_processed_chr18.txt",
#                            format="ncol")
#
# HUVEC_chr18 = ig.Graph.Load("/Users/GBS/Master/HiC-Data/Processed_Data_Edgelist/HiC_from_Jonas/Chr18/HUVEC_processed_chr18.txt", format="ncol")
#
# IMR90 = ig.Graph.Load("/Users/GBS/Master/HiC-Data/Processed_Data_Edgelist/HiC_from_Jonas/FullGenome/IMR90_processed.txt", format="ncol")

# K562_chr18.save("K562_chr18", format = "ncol") Saving doesn't work?
# Finding degree of graph
# degree_K562_chr18 = K562_chr18.degree()
# print(degree_K562_chr18)

# Finding betweenness of edges (same as using .es but our graph is not directed)
# edge_betweenness_K562_chr18 = K562_chr18.edge_betweenness()
# print(edge_betweenness_K562_chr18)

# Finding betweenness of vertices (.vs: nodes)
# nodes_betweenness_K562_chr18 = K562_chr18.vs.betweenness()
# print(nodes_betweenness_K562_chr18)

# Finding adjacency matrix for graph
# adjacency_K562_chr18 = K562_chr18.get_adjacency()
# print(adjacency_K562_chr18)

# Is our graph directed:
# print("Our graph is directed:", K562_chr18.is_directed())

# Testing if the Cairo package works (it works)
# g = ig.Graph.Famous("petersen")
# plot(g)
# And it works on my data
# plot(K562_chr18, layout = "fr", directed="false") # 2D layouts can be: circle, drl, fr, kk, lgl, random, rt


#


# plot(IMR90, layout = "")
# plot(HUVEC_chr18, layout = "fr")


# Frequency distribution test

# Må ha frequency på Y-aksen og Degree på X-aksen
# Må ogsp ha dotplot og ikke histogram, det ser ut som de har log-transformert grafen i paperet.
# Sjekk paper for å se hvordan det skal se ut: https://pubmed.ncbi.nlm.nih.gov/27618581/

# Plot this using Altair? Looked nice on OMgenomics.

# Something is wrong with iGraph, troubleshoot later: https://github.com/scverse/scanpy/issues/961

# ig.add_vertices(5)
# ig.add_edges([(0,1), (0,2), (0,3), (1,2), (2,3), (3,4)])

# Compute the degree distribution using the Counter class
# degree_values = g.degree()
# degree_counts = Counter(degree_values)

# Create a dataframe with the degree and frequency data
# df = pd.DataFrame({"degree": list(degree_counts.keys()), "frequency": list(degree_counts.values())})

# Plot the degree distribution using the Chart and Scatter classes
# chart = alt.Chart(df).mark_circle(size=60).encode(
#     x=alt.X("degree:Q", scale=alt.Scale(type="log")),
#     y=alt.Y("frequency:Q", scale=alt.Scale(type="log")),
# )
# chart.show()

# bins = 20
# # log_data = np.log(K562_chr18)
# plt.scatter(K562_chr18.degree(), bins)
# plt.semilogx()
# plt.xscale("log")
# plt.ylim(0.1, 100)
# plt.xlim(9, max(xscale))
# plt.ylabel("Degree")
# plt.xlabel("Frequency")
# plt.show()

# Something like this for degree distribution?

# plt.scatter(degree_counts.keys(), degree_counts.freq(), linewidth=0)
# plt.semilogx(degree_counts.keys(), degree_counts.values())
# plt.semilogy(degree_counts.keys, degree_counts.values())
# plt.xlabel("Degree")
# plt.ylabel("Frequency")
# plt.show()


# After we standardize nodes, we can run community detection (e.g fast greedy first), then look at
# clustering, degree distritution, betweeness centrality and other network metrics.

# communities = g.communnity_fastgreedy()
# print(communities)
# communitites_list = communitites.as_clustering() # communitites as list of vertices
# The index represents the ID of each edge:
# membership = communitites.membership
# print(membership)

# Plot the graph, coloring the vertices by their community
# plot_options = {
#     "vertex_color": communities.membership,
#     "vertex_label": g.vs["name"],
#     "vertex_size": 30,
#     "edge_color": "lightgray",
#     "margin": 20
# }
# ig.plot(g, **plot_options)
