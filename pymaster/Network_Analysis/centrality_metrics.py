# Import modules
from igraph import plot

from Graph_Processing import graph_instances as gi
from Graph_Processing import graph_metrics as gm
from Graph_Processing import graph_generator as gg
import igraph as ig
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
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


class Directories:
    base_path = path.Path("/Users/GBS/Master/Figures")
    degree_path = base_path / "Degree"

    if not degree_path.exists():
        degree_path.mkdir(parents=True)


class DegreeCentrality:

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


class PlotDegreeDistribution(DegreeCentrality):

    def plot_degree_distribution(self, plot_type='scatter', save_as=None, normalize=False):
        # Define colors
        colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
        max_frequency = -np.inf
        max_log_frequency = -np.inf
        min_log_frequency = np.inf
        maximum_x = -np.inf
        minimum_y = np.inf

        for idx, (graph_name, graph) in enumerate(self.graph_dict.items()):
            # Decide whether to normalize
            if normalize:
                degree_sequence = graph.vs["normalized_degree"]
            else:
                degree_sequence = graph.degree()

            # Unique degrees and their counts
            degrees, counts = np.unique(degree_sequence, return_counts=True)

            # Convert counts to frequencies
            frequencies = counts / counts.sum()

            # Log transform frequencies
            log_freq = np.log10(frequencies)

            # Split graph_name for legend
            legend_name = '_'.join(graph_name.split('_')[1:3])

            # Check plot type and plot accordingly
            if plot_type == 'scatter':
                plt.scatter(degrees, log_freq, color=colors[idx % len(colors)], alpha=0.5, label=legend_name)
            elif plot_type == 'line':
                plt.plot(degrees, log_freq, color=colors[idx % len(colors)], label=legend_name)

            # Update max frequency and max log-frequency
            max_frequency = max(max_frequency, np.max(frequencies))
            max_log_frequency = max(max_log_frequency, np.max(log_freq))
            min_log_frequency = min(min_log_frequency, np.min(log_freq))
            maximum_x = max(maximum_x, np.max(degrees))
            minimum_y = min(minimum_y, np.min(log_freq))

        # Apply locator to y-axis
        majorlocator = MultipleLocator(1)
        minorlocator = AutoMinorLocator(10)
        # TODO: Why does my minor ticks not show up? wtf.. (Are they outside the plot?)
        plt.gca().yaxis.set_major_locator(majorlocator)
        plt.gca().yaxis.set_minor_locator(minorlocator)
        plt.gca().tick_params(axis='y', which='minor', bottom=True)

        # Set y ticks to represent frequency in a 10^n format with more granularity
        plt.yticks(ticks=np.arange(max_log_frequency, min_log_frequency - 1, -1),
                   labels=[f'$10^{{{int(tick)}}}$' for tick in np.arange(max_log_frequency, min_log_frequency - 1, -1)])

        # Set log scale on x-axis
        plt.xscale('log')

        plt.title("Degree Distribution")
        plt.xlabel("Degree")
        plt.ylabel("Frequency")

        # Add a legend
        plt.legend(loc='upper right')

        # Add max value text below figure
        max_x = maximum_x
        min_y = minimum_y
        max_y = max_log_frequency
        plt.text(0.03, 0.03, f'Max Degree: {max_x}'
                             f'\nMax frequency: $10^{{{max_y:.1f}}}$'
                             f'\nMin Frequency: $10^{{{min_y:.1f}}}$',
                 transform=plt.gca().transAxes)

        # Save the plot conditionally
        if save_as:
            plt.savefig(save_as, dpi=300, format='png')

        # Show the plot
        plt.show()


def degree():
    graph_dict = gi.intra_1mb_graphs()
    degree_instance = PlotDegreeDistribution(graph_dict)
    degree_instance.calculate_degree()
    degree_instance.normalize_degree()
    degree_instance.plot_degree_distribution(save_as=Directories.degree_path / "degree_distribution.png", normalize=False)

# degree()


class PlotDegreeNetwork(DegreeCentrality):

    def plot_graph(self):
        for graph_name, graph in self.graph_dict.items():
            graph_name = '_'.join(graph_name.split('_')[1:3])
            graph = graph.simplify()
            visual_style = {"vertex_size": 0.7, "bbox": (1000, 1000), "margin": 20, "vertex_color": [DegreeCentrality.normalize_degree(graph)]}  # [degree / max(graph.degree()) for degree in graph.degree()]}
            ig.plot(graph, **visual_style, target=f"Degree Network of {graph_name}.png")
            fig, ax = plt.subplots()
            ig.plot(graph, bbox=(0, 0, 300, 300), target=ax, vertex_size=0.1, edge_width=0.5)  # node_color=filtered_graph.vs["chromosome"])
            plt.show()

def plot_degree_network():
    graphs = gi.intra_1mb_graphs()
    filter_instance = gm.FilterGraphs(graphs)
    filter_instance.filter_graphs(cell_lines=["MCF10"], interaction_type="intra")
    graph_dict = filter_instance.graph_dict
    lcc_instance = gm.LargestComponent(graph_dict)
    lcc_instance.find_lcc()
    graph_dict = lcc_instance.graph_dict
    degree_instance = PlotDegreeNetwork(graph_dict)
    degree_instance.calculate_degree()
    degree_instance.normalize_degree()
    degree_instance.plot_graph()
plot_degree_network()

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
