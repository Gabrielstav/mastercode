# Import modules
import igraph as ig
from matplotlib import pyplot as plt
from pathlib import Path
import re as re
import types as types
from typing import Union

# Set backend of igraph to matplotlib:
ig.config["plotting.backend"] = "matplotlib"


# Maybe used later?
# import math
# import networkx as nx
# import numpy as np
# from collections import Counter
# import altair as alt
# import pandas as pd
# import seaborn as sb
# import os as os
# import scipy.stats as stats
# import pandas as pd
# import netZooPy as nz


# import pycairo as pc


# TODO: Eventually make class to read in edgelists, pass those to another class that makes graph objects, and then pass those to a metric class that calculates the metrics we want to plot.
#       This class then passes the metrics to a plotting class that plots the metrics. Which can be saved to an output directory in the end.
#       Also need to include overlapping nodes into networks, since we do not want to compare nodes not overlapping, but for this we need filtering classes first.
#       Big unknown is the overlap between the networks, since we do not know how the overlap from full ganome data yet. But we can use synthetic data to test this.


# TODO: The chromosome filtering wont work on full genome data, because now the graph name contains the chromosome, but it wont for full genome data.
#      So we need to make a class that can filter the graphs based on the chromosome, and then we can use the graph name to filter the graphs.
#      Need to figure out the filtering when fg data is available, can probably just ook at node names, since they are always containing "chr".


class CreateGraphsFromDirectory:
    """
    Class to create igraph objects from edgelists in a directory.
    Set exclude_weighted = False in from_edgelist method to create graphs from weighted edgelists.
    """

    def __init__(self, input_directory):
        self.input_directory = input_directory
        self.graph_dict = {}

    def get_files(self, pattern):
        files = (file_path for file_path in self.input_directory.rglob(pattern) if file_path.is_file())
        files_list = [str(file_path.resolve()) for file_path in files]
        return files_list

    def from_edgelists(self, pattern='**/*edgelist*.txt', exclude_weighted=True):
        if exclude_weighted:
            pattern = f"**/*[^weighted]*edgelist*.txt"
        edgelist_files = self.get_files(pattern)

        for edgelist in edgelist_files:
            graph_name = re.findall(r".*/(.+)_edgelist\.txt$", str(edgelist))[0]
            self.graph_dict[graph_name] = ig.Graph.Load(edgelist, format="ncol", directed=False)

    def from_weighted_edgelists(self, pattern='**/*weighted*.txt'):
        weighted_files = self.get_files(pattern)

        for edgelist in weighted_files:
            graph_name = re.findall(r".*/(.+)_weighted_edgelist\.txt$", str(edgelist))[0]
            self.graph_dict[graph_name] = ig.Graph.Load(edgelist, format="ncol", directed=False)

    # rename to with_filter
    def filter_graphs(self, chromosomes=None, resolutions=None):
        if chromosomes:
            if isinstance(chromosomes, str):
                chromosomes = [chromosomes]
            chromosomes = set(chromosomes)
        if resolutions:
            if isinstance(resolutions, str):
                resolutions = [resolutions]
            resolutions = set(resolutions)

        filtered_graph_dict = {}

        for graph_name, graph in self.graph_dict.items():
            if chromosomes is not None and not any(chrom in graph_name for chrom in chromosomes):
                continue

            if resolutions is not None:
                graph_name_resolution = re.search(r"_([^_]+)$", graph_name).group(1)
                if graph_name_resolution not in resolutions:
                    continue

            filtered_graph_dict[graph_name] = graph

        return filtered_graph_dict


################################
# Creating graphs from edgelists
################################

# Creating graph objects from raw and normalized edge lists, as well as filtered on resolution and chromosome:
# Move this to another module later.

def chr18_inc_norm_graphs():
    root_dir = Path("/Users/GBS/Master/HiC-Data/Pipeline_out/chr18_INC/chr18_norm/edgelists")
    graph_creator = CreateGraphsFromDirectory(root_dir)
    graph_creator.from_edgelists()
    chr18_inc_graphs_norm = graph_creator.graph_dict
    return chr18_inc_graphs_norm


# print(chr18_inc_norm_graphs())

def chr18_inc_raw_graphs():
    root_dir = Path("/Users/GBS/Master/HiC-Data/Pipeline_out/chr18_INC/chr18_raw/edgelists")
    graph_creator = CreateGraphsFromDirectory(root_dir)
    graph_creator.from_edgelists()
    chr18_inc_graphs_raw = graph_creator.graph_dict
    return chr18_inc_graphs_raw


# Filter on resolution and chromosome:
def chr18_inc_raw_graphs_50kb():
    root_dir = Path("/Users/GBS/Master/HiC-Data/Pipeline_out/chr18_INC/chr18_raw/edgelists")
    graph_creator = CreateGraphsFromDirectory(root_dir)
    graph_creator.from_edgelists()
    filtered_graphs = graph_creator.filter_graphs(resolutions=["50000"], chromosomes="chr18")
    return filtered_graphs


# TODO: Make "find largest connected component" here so it can be passed to NetworkMetrics class.
#      Or just do this in the NetworkMetrics class? Write upsides/downsides with this.

class LargestComponent:
    """
    Class to find the LCC from a given graph object.
    Pass any graph obj dict and return largest conected component.
    """

    def __init__(self, graph_dict_or_function):
        if isinstance(graph_dict_or_function, types.FunctionType):
            self.graph_dict = graph_dict_or_function()
        else:
            self.graph_dict = graph_dict_or_function

    def largest_component(self):
        largest_component_dict = {}
        for graph_name, graph in self.graph_dict.items():
            largest_component = graph.connected_components().giant()  # or g.clusters().giant() or graph.components().giant()?
            largest_component_dict[graph_name] = largest_component
        return largest_component_dict







# print(LargestComponent(chr18_inc_raw_graphs_50kb()).largest_component())


class NetworkMetrics:
    """
    Class to calculate network metrics for a given graph object.
    """

    def __init__(self, graph_dict_or_function, metrics=None):
        if isinstance(graph_dict_or_function, types.FunctionType):
            self.graph_dict = graph_dict_or_function()
        else:
            self.graph_dict = graph_dict_or_function
        self.metrics_dict = metrics if metrics is not None else {}

    available_metrics = {
        "size": lambda g: g.vcount(),
        "density": lambda g: g.density(),
        "fg_communities": lambda g: g.community_fastgreedy(),
        "community": lambda g: g.community_multilevel(),
        "assortativity": lambda g: g.assortativity_degree(),
        "betweenness": lambda g: g.betweenness(),
        "degree": lambda g: g.degree(),
        "closeness": lambda g: g.closeness(),
        "eigen_centrality": lambda g: g.eigenvector_centrality(),
        "shortest_path_length": lambda g: g.shortest_paths(),
        "radius": lambda g: g.radius(),
        "diameter": lambda g: g.diameter(),
        "average_path_length": lambda g: g.average_path_length(),
        "clustering_coefficient": lambda g: g.transitivity_undirected(),
        "jaccard_coefficient": lambda g: g.similarity_jaccard(),
        "pagerank": lambda g: g.pagerank(),
        "fg_modularity": lambda g: g.community_fastgreedy().as_clustering(),
        "average_neighbor_degree": lambda g: g.average_neighbor_degree(),
        "average_degree_connectivity": lambda g: g.average_degree_connectivity(),
        "average_clustering": lambda g: g.transitivity_avglocal_undirected(),
        "average_local_clustering": lambda g: g.transitivity_local_undirected(),
        "transitivity_undirected": lambda g: g.transitivity_undirected(),
    }

    # Calculate metrics for all graphs passed to class
    def calculate_metrics(self, selected_metrics=None):
        if selected_metrics is None:
            selected_metrics = self.available_metrics.keys()

        calculated_metrics = {}

        for graph_name, graph in self.graph_dict.items():
            metrics_for_graph = {}
            for metric in selected_metrics:
                metric_function = self.available_metrics[metric]
                metrics_for_graph[metric] = metric_function(graph)
            calculated_metrics[graph_name] = metrics_for_graph

        return calculated_metrics

    @classmethod
    def get_metrics(cls, graph_dict_or_function=None, chromosome=None, resolution=None, metrics: Union[None, dict, list] = None, root_dir=None):  # ,**kwargs):
        """
        Filter metrics from dict, function returning dict or from root directory containing edge lists
        :param graph_dict_or_function: input as dictionary or function returning dictionary
        :param chromosome: chromosome to filter on
        :param resolution: specific resolution to filter on
        :param metrics: network metrics to calculate
        :param root_dir: root dir containing edge lists (if calculating metrics from edge lists)
        :return: dict containing graph objects and metrics
        """

        # Create an instance of NetworkMetrics
        # network_metrics = cls(graph_dict_or_function)

        graph_dict = {}

        # If graph_dict_or_function is a function, call it to get the graph_dict
        if callable(graph_dict_or_function):
            graph_dict = graph_dict_or_function()
        # If graph_dict_or_function is a dictionary, use it directly
        elif isinstance(graph_dict_or_function, dict):
            graph_dict = graph_dict_or_function

        # If root_dir is provided, create graph_dict from the directory
        if root_dir is not None:
            graph_creator = CreateGraphsFromDirectory(root_dir)
            graph_creator.from_edgelists()
            graph_creator.filter_graphs(chromosomes=chromosome, resolutions=resolution)
            graph_dict = graph_creator.graph_dict

        # If metrics is None, use all available metrics
        if metrics is None:
            metrics = cls.available_metrics
        # If metrics is a list, filter the default_metrics dictionary
        elif isinstance(metrics, list):
            metrics = {metric.lower(): cls.available_metrics[metric.lower()] for metric in metrics}

        metrics_data = {}
        for graph_name, graph in graph_dict.items():
            metrics_data[graph_name] = {}
            for metric_name, metric_function in metrics.items():
                metric_value = metric_function(graph)
                metrics_data[graph_name][metric_name] = metric_value

        return metrics_data

    @classmethod
    def print_metrics(cls, metrics_dict, metric_names=None):
        if metric_names is None:
            metric_names = []

        for graph_name, metrics in metrics_dict.items():
            print(f"Metrics for {graph_name}:")
            for metric_name, metric_value in metrics.items():
                if metric_name in metric_names and callable(metric_value):
                    try:
                        metric_value = metric_value()
                    except TypeError:
                        if metric_name == "fg_modularity":
                            membership = metric_value.community_multilevel().membership
                            metric_value = metric_value.modularity(membership)
                        else:
                            raise ValueError(f"Unsupported metric '{metric_name}' with additional arguments")
                print(f"  {metric_name}: {metric_value}")
            print()


###################
# Calculate metrics
###################

# From function returning dictionary containing graph objects (or from dictionary):
def chr18_50kb_metrics():
    metrics = NetworkMetrics.get_metrics(
        graph_dict_or_function=chr18_inc_raw_graphs_50kb,
        metrics=[
            "size", "fg_communities", "assortativity", "clustering_coefficient", "radius", "diameter", "average_path_length", "fg_modularity"
        ])
    return NetworkMetrics.print_metrics(metrics)


chr18_50kb_metrics()


# Or calculate metrics from root directory containing edge lists:
def chr18_size_mod_from_directory():
    metrics = NetworkMetrics.get_metrics(chromosome="chr18", resolution="50000", metrics=["size", "fg_communities"], root_dir=Path("/Users/GBS/Master/HiC-Data/Pipeline_out/chr18_INC/chr18_raw/edgelists"))
    return NetworkMetrics.print_metrics(metrics)


# chr18_size_mod_from_directory()


class largest_connected_component:
    """
    This class should be used to calculate metrics for the largest connected components of a graph.
    Can atomatically do community detection, and find the largest connected component for all graphs in CreateGraphsFromDirectory?
    Or you can pass graph objects to the class, and it will calculate metrics for the largest connected component.
    Maybe it could filter out all smaller communities, name them based on size, and do metrics for each community?
    Another question is what community detection algorithms to use?
    """


# TODO: Make differential community detection class that takes two graphs, and calculates the difference in community detection metrics
#   Class methods can be the different metrics (Jaccard, Alpaca, genomic(?)) and the class can be initialized with the two graphs?
#   Or, we make separate class for each metric? This allows for more flexibility, but also more code duplication.


class JaccardIndex:
    """
    Calculate Jaccard index for two graphs
    """

    def __init__(self, graph1, graph2):
        self.graph1 = graph1
        self.graph2 = graph2

    def calculate(self):
        pass


class DifferentialCommunityDetection:
    """
    Calculate the difference in community detection metrics between two graphs using Alpaca R package
    """

    def __init__(self, graph1, graph2):
        self.graph1 = graph1
        self.graph2 = graph2

    def calculate(self):
        pass


# TODO: Make a class for plotting graphs, and compare largest connected components with whole network (move to network_plots.py later)

class plot_lcc_to_graph_ratio:
    """
    Plot the ratio of the size of the largest connected component to the size of the whole graph
    pass graph object and get out stacked bar chart with the ratio for each network
    """
    def something(self):
        pass


# def plot_chr18_50kb():
#     h = ig.Graph.Load("/Users/GBS/Master/HiC-Data/Pipeline_out/chr18_INC/chr18_raw/edgelists/chr18_hires_20000_edgelist.txt", format="ncol", directed=False)  # very messy
#     ig.plot(h, edge_width=0.07, node_size=0.5, node_color="red")
#     plt.show()
#     plt.close()


# PLotting class

# class NetworkPlots:
#     @staticmethod
#     def plot_degree_distribution(graph_dict, title=None):
#         for graph_name, graph in graph_dict.items():
#             degree_distribution = graph.degree_distribution()
#             x = [left + (width / 2) for left, _, width in degree_distribution.bins()]
#             y = [count for _, count, _ in degree_distribution.bins()]
#             plt.plot(x, y, marker='o', linestyle='-', label=graph_name)
#
#         plt.xlabel("Degree")
#         plt.ylabel("Frequency")
#         plt.legend()
#         if title:
#             plt.title(title)
#         plt.show()
#
# # Usage
# network_plots = NetworkPlots()
# network_plots.plot_degree_distribution(chr18_inc_graphs_norm, title="Degree Distribution")

#
# # TODO: This is messy RN; automate so each method takes the previous method as input, and then the last method returns the final output (plots).
# #    calculate largest connected ocmponent, and then calculate metrics for that component only. Compare with metrics for whole network?
# #
#


# Need to make the graph objects for each file read in, and then pass them to a class


# {'chr18_lowres_500000': <igraph.Graph object at 0x12122e340>, ...}

# g= ig.Graph.Load("/Users/GBS/Master/HiC-Data/Pipeline_out/chr18_INC/chr18_raw/edgelists/chr18_lowres_50000_edgelist.txt", format="ncol")
# ig.plot(g)
# # ig.plot(g, vertex_size=2, vertex_label=g.vs["name"], vertex_label_size=1, vertex_label_dist=1.5, vertex_label_color="black", layout=g.layout("kk")) # TODO: Figure out iptimal plot layout
# plt.show()

def plot_20kb_hires_raw():
    h = ig.Graph.Load("/Users/GBS/Master/HiC-Data/Pipeline_out/chr18_INC/chr18_raw/edgelists/chr18_hires_20000_edgelist.txt", format="ncol", directed=False)  # very messy
    ig.plot(h, edge_width=0.07, node_size=0.5, node_color="red")
    plt.show()
    plt.close()


def plot_20kb_hires_norm():
    g = ig.Graph.Load("/Users/GBS/Master/HiC-Data/Pipeline_out/chr18_INC/chr18_norm/edgelists/chr18_hires_iced_20000_edgelist.txt", format="ncol", directed=False)  # every node connected
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


# plot_50kb_raw()

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
    h = ig.Graph.Load("/Users/GBS/Master/HiC-Data/Pipeline_out/chr18_INC/chr18_raw/edgelists/chr18_hires_20000_edgelist.txt", format="ncol", directed=False)  # very messy
    fg_kb50_raw = h.community_fastgreedy()
    communities_kb50_raw = fg_kb50_raw.as_clustering()
    print(communities_kb50_raw.modularity)


def fg_kb20_normm():
    g = ig.Graph.Load("/Users/GBS/Master/HiC-Data/Pipeline_out/chr18_INC/chr18_norm/edgelists/chr18_hires_iced_20000_edgelist.txt", format="ncol", directed=False)  # every node connected
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

# fg_kb250_raww()
# fg_kb250_normm()

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
