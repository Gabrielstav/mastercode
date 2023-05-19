# Import modules
from Graph_Processing import graph_instances as gi
from Graph_Processing import graph_metrics as gm
from Graph_Processing import graph_generator as gg
import igraph as ig
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pathlib as path


# TODO: Degree class

# TODO: Betweenness class

# TODO: Closeness class

# TODO: Comparison class (for metrics between two graphs, KS, Jaccard?)

# Move betweenness and closeness to graph analysis? Since this is not necessarily comparative. "But it is comparative in the sense that it is comparing the same graph to itself" -copilot 24/4 23
# TODO: How to compare across cell lines? I mean one thing is to plot the betweenness centrality for each cell line, by coloring nodes and edges.
#   But the next thing is then to compare centrality metrics, so maybe distribution plots? Or boxplot of each metric for each cell line? (and compare resolution?)
#   - Needs to write out as node/edge attributes for plotting in graphs, and as pandas df with title as network for Jaccard index.


def root_directory():
    root_dir = path.Path("/Users/GBS/Master/Network_Analysis")
    return root_dir

class Degree:
    def __init__(self, graph_dict):
        self.graph_dict = graph_dict

    def calculate_degree(self):
        for graph_name, graph in self.graph_dict.items():
            graph.vs["degree"] = graph.degree()

    def normalize_degree(self):
        for graph_name, graph in self.graph_dict.items():
            max_degree = max(graph.vs["degree"])
            min_degree = min(graph.vs["degree"])
            graph.vs["normalized_degree"] = [(d - min_degree) / (max_degree - min_degree) for d in graph.vs["degree"]]


    def plot_degree_distribution(self, plot_type='scatter', save_as=None):
        # Define colors
        colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']

        for idx, (graph_name, graph) in enumerate(self.graph_dict.items()):
            # Degree sequence
            degree_sequence = graph.degree()

            # Unique degrees and their counts
            degrees, counts = np.unique(degree_sequence, return_counts=True)

            # Convert counts to frequencies
            frequencies = counts / counts.sum()

            # Convert frequencies to negative log scale
            log_freq = -np.log10(frequencies)

            # Split graph_name for legend
            legend_name = '_'.join(graph_name.split('_')[1:3])

            # Check plot type and plot accordingly
            if plot_type == 'scatter':
                plt.scatter(degrees, log_freq, color=colors[idx % len(colors)], alpha=0.5, label=legend_name)
            elif plot_type == 'line':
                plt.plot(degrees, log_freq, color=colors[idx % len(colors)], label=legend_name)

        # Set log scale on x-axis
        plt.xscale('log')

        # Set y ticks to represent frequency in a 10^n format with more granularity
        plt.yticks(ticks=-np.arange(0, np.ceil(np.max(log_freq)) + 0.1, 0.1),
                   labels=[f'$10^{{-{i:.1f}}}$' for i in np.arange(0, np.ceil(np.max(log_freq)) + 0.1, 0.1)])

        plt.title("Degree Distribution")
        plt.xlabel("Degree (log scale)")
        plt.ylabel("Frequency (negative log10 scale)")

        # Add a legend
        plt.legend(loc='upper right')

        # Add max value text below figure
        max_x = np.max(degrees)
        max_y = -np.min(log_freq)
        plt.text(0.05, 0.05, f'Max Degree: {max_x}\nMax Frequency: $10^{{-{max_y:.1f}}}$', transform=plt.gca().transAxes)

        # Show the plot
        plt.show()

        # Save the plot conditionally
        if save_as:
            if not os.path.exists(save_dir):
                os.makedirs(save_dir)
            plt.savefig(os.path.join(save_dir, f'{save_as}.png'))

    def plot_graph(self):
        for graph_name, graph in self.graph_dict.items():
            fig, ax = plt.subplots()
            visual_style = {"vertex_color": [plt.cm.Reds(deg) for deg in graph.vs["normalized_degree"]], "vertex_size": 20}
            ig.plot(graph, **visual_style, target=ax)


def degree():
    graph_dict = gi.norm_mcf7_graphs()
    degree_instance = Degree(graph_dict)
    degree_instance.calculate_degree()
    degree_instance.normalize_degree()
    degree_instance.logplot_degree_distribution()


degree()


class ClosenessCentrality:

    def __init__(self, graph_dict):
        self.graph_dict = graph_dict

    def calculate_closeness(self):
        # Store closeness centrality as node attributes?
        pass

    def normalize_closeness(self):
        # Normalize metrics for plotting in network to color nodes and edges?
        pass

    def plot_closeness(self):
        pass


class BetweennessCentrality:

    def __init__(self, graph_dict):
        self.graph_dict = graph_dict

    def calculate_betweenness(self):
        # store betweenness centrality as node attributes?
        pass

    def normalize_betweenness(self):
        # Normalize metrics for plotting in network to color nodes and edges?
        pass

    def plot_betweenness(self):
        pass


class GraphCentralityMetrics(gm.NetworkMetrics):

    def __init__(self, graph_dict):
        self.graph_dict = graph_dict

    def calculate_degree(self):
        for graph_name, graph in self.graph_dict.items():
            graph.vs["degree"] = graph.degree()
            graph.es["degree"] = graph.degree()

    def calculate_betweenness(self):
        pass

    def calculate_closeness(self):
        pass

    def as_dataframe(self):
        pass


class CompareDegree:

    def __init__(self, graph_dict):
        self.graph_dict = graph_dict

    def compare_degree(self):
        pass


class CompareBetweenness:

    def __init__(self, graph_dict):
        self.graph_dict = graph_dict

    def compare_betweenness(self):
        pass


class CompareCloseness:

    def __init__(self, graph_dict):
        self.graph_dict = graph_dict

    def compare_closeness(self):
        pass


class BetweennessCentrality:

    def __init__(self, graph_dict):
        self.graph_dict = graph_dict

    def calculate_betweenness_centrality(self):
        betweenness_centrality_dict = {}
        for graph_name, graph in self.graph_dict.items():
            betweenness_centrality_dict[graph_name] = graph.betweenness()
        return betweenness_centrality_dict

    def calculate_betweenness_centrality_lcc(self):
        betweenness_centrality_dict = {}
        for graph_name, graph in self.graph_dict.items():
            betweenness_centrality_dict[graph_name] = graph.betweenness()
        return betweenness_centrality_dict

    def print_betweenness_centrality(self):
        for graph_name, graph in self.calculate_betweenness_centrality().items():
            print(f"Betweenness centrality for: {graph_name} \n {graph}")


class ClosenessCentrality:

    def __init__(self, graph_dict):
        self.graph_dict = graph_dict

    def calculate_closeness_centrality(self):
        closeness_centrality_dict = {}
        for graph_name, graph in self.graph_dict.items():
            closeness_centrality_dict[graph_name] = graph.closeness()
        return closeness_centrality_dict

    def print_closeness_centrality(self):
        for graph_name, graph in self.calculate_closeness_centrality().items():
            print(f"Closeness centrality for: {graph_name} \n {graph}")


class Degree:
    pass


class FindUniqueNodes:
    """
    This class takes a dictionary of graph objects and finds the unique nodes in each graph by comparing all nodes
    between all graphs.
    """

    def __init__(self, graph_dict):
        self.graph_dict = graph_dict

    def find_shared_nodes(self):
        shared_nodes = set()
        graph_names = list(self.graph_dict.keys())

        for i, (graph_name1, graph1) in enumerate(self.graph_dict.items()):
            for graph_name2, graph2 in list(self.graph_dict.items())[i + 1:]:
                shared_nodes.update(set(graph1.vs) & set(graph2.vs))

        return shared_nodes

    def find_unique_nodes(self):
        unique_nodes_dict = {}
        graph_names = list(self.graph_dict.keys())

        for graph_name, graph in self.graph_dict.items():
            other_graphs = [g for g in self.graph_dict.values() if g != graph]
            unique_nodes = set(graph.vs) - set().union(*[g.vs for g in other_graphs])
            unique_nodes_dict[graph_name] = unique_nodes

        return unique_nodes_dict

    def print_shared_nodes(self):
        shared_nodes = self.find_shared_nodes()
        print("Shared Nodes:")
        for node in shared_nodes:
            print(node)

    def print_unique_nodes(self):
        unique_nodes_dict = self.find_unique_nodes()
        print("Unique Nodes:")
        for graph_name, unique_nodes in unique_nodes_dict.items():
            print(f"Graph: {graph_name}")
            for node in unique_nodes:
                print(node)
