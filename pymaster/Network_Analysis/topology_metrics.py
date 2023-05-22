# Import modules
import igraph as ig
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sn
from Graph_Processing import graph_generator as gg
from Graph_Processing import graph_metrics as gm
from Graph_Processing import graph_instances as gi
from Network_Analysis import centrality_metrics as cm


# TODO: Calculate topology using fg, louvian and infomap

# Each class takes a graph, and return the graph object with the topology calculated as an attribute?
# Or it makes a new graph? It needs to be easility plottable.
# Alpaca takes two graphs as input?

def method_lists():
    return ["fg", "lv", "im", "ld", "lp"]


class CommunityDetection:

    def __init__(self, graph_dict):
        self.methods = None
        self.graph_dict = graph_dict

    def calculate_fastreedy(self):
        for graph_name, graph in self.graph_dict.items():
            graph.vs["fg"] = graph.community_fastgreedy().as_clustering().membership
        return self.graph_dict

    def calculate_louvian(self):
        for graph_name, graph in self.graph_dict.items():
            graph.vs["louvian"] = graph.community_multilevel().membership
        return self.graph_dict

    def calculate_infomap(self):
        for graph_name, graph in self.graph_dict.items():
            graph.vs["infomap"] = graph.community_infomap().membership
        return self.graph_dict

    def calculate_leiden(self):
        for graph_name, graph in self.graph_dict.items():
            graph.vs["leiden"] = graph.community_leiden().membership
        return self.graph_dict

    def calculate_label_propagation(self):
        for graph_name, graph in self.graph_dict.items():
            graph.vs["label_propagation"] = graph.community_label_propagation().membership
        return self.graph_dict

    def print_communities(self, method=None):
        if method in method_lists():
            for graph_name, graph in self.graph_dict.items():
                print(graph_name)
                print(graph.vs[method])


class PlotTopology:

    def __init__(self, graph_dict):
        self.graph_dict = graph_dict

    def plot_community(self, method=None):
        for graph_name, graph in self.graph_dict.items():

            if method not in method_lists():
                print(f"Method not found. Please choose from the following: {method_lists()}")

            fig, ax = plt.subplots()
            ig.plot(graph, vertex_color=graph.vs[method], ax=ax)
            ax.set_title(f"{graph_name} - {method}")
            plt.show()


class CommunityDetection_1:
    def __init__(self, graph_dict):
        self.graph_dict = graph_dict
        self.methods = {
            'fg': self._fastgreedy,
            'lv': self._louvain,
            'im': self._infomap,
            'ld': self._leiden,
            'lp': self._label_propagation,
        }

    def calculate_communities(self, method):
        if method not in self.methods:
            raise ValueError(f"Invalid method. Choose from: {', '.join(self.methods.keys())}")
        for graph in self.graph_dict.values():
            self.methods[method](graph)

    @staticmethod
    def _fastgreedy(graph):
        graph.vs["fastgreedy"] = graph.community_fastgreedy().as_clustering().membership
        graph.vs["fastgreedy"] = graph.cluster_fast_greedy().membership

    @staticmethod
    def _louvain(graph):
        graph.vs["louvain"] = graph.community_multilevel().membership

    @staticmethod
    def _infomap(graph):
        graph.vs["infomap"] = graph.community_infomap().membership

    @staticmethod
    def _leiden(graph):
        graph.vs["leiden"] = graph.community_leiden().membership

    @staticmethod
    def _label_propagation(graph):
        graph.vs["label_propagation"] = graph.community_label_propagation().membership

    def print_communities(self, method=None):
        if method not in self.methods:
            raise ValueError(f"Invalid method. Choose from: {', '.join(self.methods.keys())}")
        for graph_name, graph in self.graph_dict.items():
            if method not in graph.vs.attributes():
                print(f"No communities calculated for {graph_name} with method {method}")
                continue
            print(graph_name)
            print(graph.vs[method])


class PlotTopology_1:
    def __init__(self, graph_dict):
        self.graph_dict = graph_dict

    def plot_community(self, method=None):
        if method not in CommunityDetection(None).methods:
            raise ValueError(f"Invalid method. Choose from: {', '.join(CommunityDetection(None).methods.keys())}")
        for graph_name, graph in self.graph_dict.items():
            fig, ax = plt.subplots()
            ig.plot(graph, vertex_color=graph.vs[method], ax=ax)
            ax.set_title(f"{graph_name} - {method}")
            plt.show()


def fg_comm_mcf10():
    graph_dict = gi.intra_1mb_graphs()
    filter_instance = gm.FilterGraphs(graph_dict)
    filtered = filter_instance.filter_graphs(chromosomes=["chr1"], resolutions=[1000000], cell_lines=["mcf10"], condition="intra-split-raw")
    cd = CommunityDetection(filtered)
    print(cd)
    cd.calculate_fastreedy()
    # cd.calculate_communities(method="fg")
    # print(cd.calculate_communities(method="fg"))
    cd.print_communities(method="fg")
    pt = PlotTopology_1(filtered)
    pt.plot_community(method="fg")
    return print(graph_dict)


fg_comm_mcf10()
# TODO: Ok so my implementation works but the one where I try to make the struct nices does not.
#   spend 30 minutes trying to fix it, if not, just use the working one, and make the plots of networks.

# TODO: Then tomorrow make the heatmaps of communities, to compare across cell lines.
# TODO: Then implement Alpaca and do differential community detection, and make similar plots.