# Import modules
from igraph import plot

from Graph_Processing import graph_instances as gi
from Graph_Processing import graph_metrics as gm
from Graph_Processing import graph_generator as gg
import igraph as ig
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
import matplotlib.patches as mpatches
import matplotlib.colors as clrs
from matplotlib.colors import Normalize
import matplotlib.cm as cm
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

class PlotDegreeNetwork(DegreeCentrality):

    def calculate_color(self, edge):
        for graph_name, graph in self.graph_dict.items():
            source_degree = graph.degree(edge.source)
            target_degree = graph.degree(edge.target)

            # Use the largest of the source and target degrees
            if source_degree > target_degree:
                degrees = source_degree
            else:
                degrees = target_degree

            # Use this instead to keep source and target degrees separate:
            # degrees = graph.degree(edge.source)

            # Map degree to a color
            color = plt.cm.viridis(degrees / max(graph.degree()))
            return color

    def plot_graph(self, normalize=False, color_edges=False, save_as=None, layout=None):
        for graph_name, graph in self.graph_dict.items():
            fig, ax = plt.subplots()
            graph_name = '_'.join(graph_name.split('_')[1:3])
            if normalize:
                degrees = [deg / max(graph.degree()) for deg in graph.degree()]  # Normalize degrees
            else:
                max_degree = max(graph.degree())
                degrees = [deg / max_degree for deg in graph.degree()]  # Normalize degrees for color mapping

            colors = [list(color) for color in plt.cm.viridis(degrees)]  # Convert numpy.ndarray to list of lists

            # node_size = 0.8
            # edge_width = 1
            node_size = 80 / graph.vcount()
            edge_width = 500 / graph.ecount()
            edge_colors = [self.calculate_color(edge) for edge in graph.es] if color_edges else 'gray'
            visual_style = {
                "layout": layout,
                "vertex_size": node_size,
                "edge_width": edge_width,
                "bbox": (1000, 1000),
                "margin": 20,
                "vertex_color": colors,
                "edge_color": edge_colors,
                "vertex_label_size": 2,
                "vertex_label_dist": 10,
                "vertex_label_angle": 100,
                "vertex_label_color": "black",
                "vertex_frame_color": colors,  # Remove black outline, use in circle layout
                # "vertex_label": graph.vs["location"],  # Add location as label, can be used in small graphs
            }
            ig.plot(graph, **visual_style, target=ax)
            plt.title(graph_name)  # Add a title

            # Create a colorbar
            if normalize:
                norm = Normalize(vmin=min(degrees), vmax=max(degrees))
            else:
                norm = Normalize(vmin=min(graph.degree()), vmax=max(graph.degree()))
            sm = plt.cm.ScalarMappable(cmap=plt.cm.viridis, norm=norm)
            sm.set_array([])
            plt.colorbar(sm, ax=ax, orientation='vertical', label='Degree')

            # Save the plot conditionally
            if save_as:
                plt.savefig(save_as, dpi=300, format='png')

            plt.show()

def plot_degree_network():
    graphs = gi.intra_50kb_graphs()
    filter_instance = gm.FilterGraphs(graphs)
    filter_instance.filter_graphs(cell_lines=["mcf10"], interaction_type="intra")  # , chromosomes=["chr1"])
    graph_dict = filter_instance.graph_dict
    print(graph_dict)
    lcc_instance = gm.LargestComponent(graph_dict)
    lcc = lcc_instance.find_lcc()
    print(lcc)
    degree_instance = PlotDegreeNetwork(graph_dict) # Use LCC for withtin-chromosome plots
    degree_instance.calculate_degree()
    degree_instance.normalize_degree()
    degree_instance.plot_graph(save_as=None, normalize=False, color_edges=True, layout="fr")
# plot_degree_network()


class ClosenessCentrality:

    def __init__(self, graph_dict):
        self.graph_dict = graph_dict

    def calculate_closeness(self):
        for graph_name, graph in self.graph_dict.items():
            graph.vs["closeness"] = graph.closeness()

    def normalize_closeness(self):
        for graph_name, graph in self.graph_dict.items():
            max_closeness = max(graph.closeness())
            graph.vs["closeness"] = [closeness / max_closeness for closeness in graph.closeness()]

class PlotClosenessNetwork(ClosenessCentrality):

    def calculate_color(self, edge):
        for graph_name, graph in self.graph_dict.items():
            source_closeness = graph.closeness(edge.source)
            target_closeness = graph.closeness(edge.target)

            # Use the largest of the source and target degrees
            if source_closeness > target_closeness:
                closeness = source_closeness
            else:
                closeness = target_closeness

            # Use this instead to keep source and target degrees separate:
            # degrees = graph.degree(edge.source)

            # Map degree to a color
            color = plt.cm.viridis(closeness / max(graph.closeness()))
            return color

    def plot_closeness(self, normalize=False, color_edges=False, save_as=None, layout=None):
        for graph_name, graph in self.graph_dict.items():
            fig, ax = plt.subplots()
            graph_name = '_'.join(graph_name.split('_')[1:3])
            if normalize:
                closeness = [closeness / max(graph.closeness()) for closeness in graph.closeness()]
            else:
                closeness = graph.closeness()

            colors = [list(color) for color in plt.cm.viridis(closeness)]  # Convert numpy.ndarray to list of lists

            node_size = 50 / graph.vcount()
            edge_width = 350 / graph.ecount()
            edge_colors = [self.calculate_color(edge) for edge in graph.es] if color_edges else 'gray'
            visual_style = {
                "layout": layout,
                "vertex_size": node_size,
                "edge_width": edge_width,
                "bbox": (1000, 1000),
                "margin": 20,
                "vertex_color": colors,
                "edge_color": edge_colors,
                "vertex_label_size": 2,
                "vertex_label_dist": 10,
                "vertex_label_angle": 100,
                "vertex_label_color": "black",
                "vertex_frame_color": colors,  # Remove black outline, use in circle layout
                # "vertex_label": graph.vs["location"],  # Add location as label, can be used in small graphs
            }
            ig.plot(graph, **visual_style, target=ax)
            plt.title(graph_name)  # Add a title

            # Create a colorbar
            if normalize:
                norm = Normalize(vmin=min(closeness), vmax=max(closeness))
            else:
                norm = Normalize(vmin=min(graph.closeness()), vmax=max(graph.closeness()))
            sm = plt.cm.ScalarMappable(cmap=plt.cm.viridis, norm=norm)
            sm.set_array([])
            plt.colorbar(sm, ax=ax, orientation='vertical', label='Closeness')

            # Save the plot conditionally
            if save_as:
                plt.savefig(save_as, dpi=300, format='png')

            plt.show()


def plot_closeness_network():
    graphs = gi.intra_1mb_graphs()
    filter_instance = gm.FilterGraphs(graphs)
    filter_instance.filter_graphs(cell_lines=["mcf7"], interaction_type="intra", chromosomes=["chr1"])
    graph_dict = filter_instance.graph_dict
    print(graph_dict)
    lcc_instance = gm.LargestComponent(graph_dict)
    lcc = lcc_instance.find_lcc()
    print(lcc)
    closeness_instance = PlotClosenessNetwork(lcc) # Use LCC for withtin-chromosome plots
    closeness_instance.plot_closeness(save_as=None, normalize=False, color_edges=True, layout="fr")
# plot_closeness_network()


class BetweennessCentrality:

    def __init__(self, graph_dict):
        self.graph_dict = graph_dict

    def calculate_betweenness(self):
        for graph_name, graph in self.graph_dict.items():
            graph.vs["betweenness"] = graph.betweenness()

    def normalize_betweenness(self):
        for graph_name, graph in self.graph_dict.items():
            max_betweenness = max(graph.betweenness())
            graph.vs["betweenness"] = [betweenness / max_betweenness for betweenness in graph.betweenness()]

class PlotBetweennessNetwork(BetweennessCentrality):

    def calculate_color(self, edge):
        for graph_name, graph in self.graph_dict.items():
            source_betweenness = graph.betweenness(edge.source)
            target_betweenness = graph.betweenness(edge.target)

            # Use the largest of the source and target degrees
            if source_betweenness > target_betweenness:
                betweenness = source_betweenness
            else:
                betweenness = target_betweenness

            # Use this instead to keep source and target degrees separate:
            # degrees = graph.degree(edge.source)

            # Map degree to a color
            color = plt.cm.viridis(betweenness / max(graph.betweenness()))
            return color

    def plot_betweenness(self, normalize=False, color_edges=False, save_as=None, layout=None):
        for graph_name, graph in self.graph_dict.items():
            fig, ax = plt.subplots()
            graph_name = '_'.join(graph_name.split('_')[1:3])
            if normalize:
                betweenness = [betweenness / max(graph.betweenness()) for betweenness in graph.betweenness()]
            else:
                betweenness = graph.betweenness()

            colors = [list(color) for color in plt.cm.viridis(betweenness)]

            node_size = 50 / graph.vcount()
            edge_width = 225 / graph.ecount()
            edge_colors = [self.calculate_color(edge) for edge in graph.es] if color_edges else 'gray'
            visual_style = {
                "layout": layout,
                "vertex_size": node_size,
                "edge_width": edge_width,
                "bbox": (1000, 1000),
                "margin": 20,
                "vertex_color": colors,
                "edge_color": edge_colors,
                "vertex_label_size": 2,
                "vertex_label_dist": 10,
                "vertex_label_angle": 100,
                "vertex_label_color": "black",
                "vertex_frame_color": colors,  # Remove black outline, use in circle layout
                # "vertex_label": graph.vs["location"],  # Add location as label, can be used in small graphs
            }

            ig.plot(graph, **visual_style, target=ax)
            plt.title(graph_name)  # Add a title

            # Create a colorbar
            if normalize:
                norm = Normalize(vmin=min(betweenness), vmax=max(betweenness))
            else:
                norm = Normalize(vmin=min(graph.betweenness()), vmax=max(graph.betweenness()))

            sm = plt.cm.ScalarMappable(cmap=plt.cm.viridis, norm=norm)
            sm.set_array([])
            plt.colorbar(sm, ax=ax, orientation='vertical', label='Betweenness')

            # Save the plot conditionally
            if save_as:
                plt.savefig(save_as, dpi=300, format='png')

            plt.show()


def plot_betweenness_network():
    graphs = gi.intra_1mb_graphs()
    filter_instance = gm.FilterGraphs(graphs)
    filter_instance.filter_graphs(cell_lines=["gsm2824367"], interaction_type="intra", chromosomes=["chr1"])
    graph_dict = filter_instance.graph_dict
    print(graph_dict)
    lcc_instance = gm.LargestComponent(graph_dict)
    lcc = lcc_instance.find_lcc()
    print(lcc)
    betweenness_instance = PlotBetweennessNetwork(lcc) # Use LCC for withtin-chromosome plots
    betweenness_instance.plot_betweenness(save_as=None, normalize=True, color_edges=False, layout="fr")
# plot_betweenness_network()


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

            # TODO: If normalized, use normalized min and max degree as x- and y-limits
            # Update max frequency and max log-frequency
            max_frequency = max(max_frequency, np.max(frequencies))
            max_log_frequency = max(max_log_frequency, np.max(log_freq))
            min_log_frequency = min(min_log_frequency, np.min(log_freq))
            maximum_x = max(maximum_x, np.max(degrees))
            minimum_y = min(minimum_y, np.min(log_freq))

        # Apply locator to y-axis
        majorlocator = MultipleLocator(1)
        minorlocator = AutoMinorLocator(10)
        # TODO: Fix - Why do minor ticks not show up?
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

def plot_degree():
    graph_dict = gi.intra_1mb_graphs()
    degree_instance = PlotDegreeDistribution(graph_dict)
    degree_instance.calculate_degree()
    degree_instance.normalize_degree()
    degree_instance.plot_degree_distribution(save_as=Directories.degree_path / "degree_distribution.png", normalize=False)

# plot_degree()


class ClosenessDistribution(ClosenessCentrality):

    def __init__(self, graph_dict):
        super().__init__(graph_dict)
        self.closeness_distribution_dict = {}

    def calculate_closeness_distribution(self, num_bins=20):
        for graph_name, graph in self.graph_dict.items():
            closeness_scores = graph.closeness()
            frequencies, bin_edges = np.histogram(closeness_scores, bins=num_bins)
            bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
            self.closeness_distribution_dict[graph_name] = (bin_centers, frequencies)

    def plot_closeness_distribution(self, plot_type='scatter', save_as=None):

        colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']

        for idx, (graph_name, (bin_centers, frequencies)) in enumerate(self.closeness_distribution_dict.items()):

            frequencies = frequencies / frequencies.sum()  # Convert counts to frequencies
            log_frequencies = np.log10(frequencies + 1e-10)  # Add a small constant to avoid -inf values from log(0)
            legend_name = '_'.join(graph_name.split('_')[1:3])

            # Exclude zero-frequency bins
            nonzero_indices = np.where(frequencies > 0)
            bin_centers = bin_centers[nonzero_indices]
            log_frequencies = log_frequencies[nonzero_indices]  # Exclude zero-frequency bins

            if plot_type == 'scatter':
                plt.scatter(bin_centers, log_frequencies, color=colors[idx % len(colors)], alpha=0.5, label=legend_name)
            elif plot_type == 'line':
                plt.plot(bin_centers, log_frequencies, color=colors[idx % len(colors)], label=legend_name)

            max_log_freq = np.max(log_frequencies)
            min_log_freq = np.min(log_frequencies)
            plt.yticks(ticks=np.arange(max_log_freq, min_log_freq - 1, -1),
                       labels=[f'$10^{{{int(tick)}}}$' for tick in np.arange(max_log_freq, min_log_freq - 1, -1)])

        plt.xlabel('Closeness')
        plt.ylabel('Log Frequency')
        plt.legend()
        plt.title('Log closeness distribution')

        if save_as:
            plt.savefig(save_as, dpi=300, format='png')
        plt.show()

def closeness_plot():
    graph_dict = gi.intra_1mb_graphs()
    closeness_instance = ClosenessDistribution(graph_dict)
    closeness_instance.calculate_closeness()
    closeness_instance.calculate_closeness_distribution()
    closeness_instance.normalize_closeness()
    # closeness_instance.plot_closeness_distribution(save_as=None)
    closeness_instance.plot_closeness_distribution(save_as=None)
# closeness_plot()

class BetweennessDistribution(BetweennessCentrality):

    def __init__(self, graph_dict):
        super().__init__(graph_dict)
        self.betweenness_distribution_dict = {}

    def calculate_betweenness_distribution(self, num_bins=100):
        for graph_name, graph in self.graph_dict.items():
            betweenness_scores = graph.betweenness()
            frequencies, bin_edges = np.histogram(betweenness_scores, bins=num_bins)
            bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
            self.betweenness_distribution_dict[graph_name] = (bin_centers, frequencies)

    def plot_betweenness_distribution(self, plot_type='scatter', save_as=None):

        colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']

        for idx, (graph_name, (bin_centers, frequencies)) in enumerate(self.betweenness_distribution_dict.items()):

            frequencies = frequencies / frequencies.sum()  # Convert counts to frequencies
            log_frequencies = np.log10(frequencies + 1e-10)  # Add a small constant to avoid -inf values from log(0)
            legend_name = '_'.join(graph_name.split('_')[1:3])

            # Exclude zero-frequency bins
            nonzero_indices = np.where(frequencies > 0)
            bin_centers = bin_centers[nonzero_indices]
            log_frequencies = log_frequencies[nonzero_indices]  # Exclude zero-frequency bins

            if plot_type == 'scatter':
                plt.scatter(bin_centers, log_frequencies, color=colors[idx % len(colors)], alpha=0.5, label=legend_name)
            elif plot_type == 'line':
                plt.plot(bin_centers, log_frequencies, color=colors[idx % len(colors)], label=legend_name)

            max_log_freq = np.max(log_frequencies)
            min_log_freq = np.min(log_frequencies)
            plt.yticks(ticks=np.arange(max_log_freq, min_log_freq - 1, -1),
                       labels=[f'$10^{{{int(tick)}}}$' for tick in np.arange(max_log_freq, min_log_freq - 1, -1)])

        plt.xscale("log")
        plt.xlabel('Betweenness')
        plt.ylabel('Log Frequency')
        plt.legend()
        plt.title('Log betweenness distribution')

        if save_as:
            plt.savefig(save_as, dpi=300, format='png')
        plt.show()

def betweenness_plot():
    graph_dict = gi.intra_1mb_graphs()
    betweenness_instance = BetweennessDistribution(graph_dict)
    betweenness_instance.calculate_betweenness()
    betweenness_instance.calculate_betweenness_distribution()
    betweenness_instance.normalize_betweenness()
    betweenness_instance.plot_betweenness_distribution(save_as=None)

# betweenness_plot()


class CentralityCorrelation(BetweennessCentrality, ClosenessCentrality, DegreeCentrality):

    def __init__(self, graph_dict):
        super().__init__(graph_dict)
        self.centrality_correlation_dict = {}

    def calculate_centrality_correlation(self):
        for graph_name, graph in self.graph_dict.items():
            betweenness_scores = graph.betweenness()
            closeness_scores = graph.closeness()
            degree_scores = graph.degree()
            self.centrality_correlation_dict[graph_name] = {
                'betweenness': betweenness_scores,
                'closeness': closeness_scores,
                'degree': degree_scores
            }

    def plot_centrality_correlation(self, metric1, metric2, save_as=None):
        colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
        for idx, (graph_name, metrics) in enumerate(self.centrality_correlation_dict.items()):
            scores1 = metrics[metric1]
            scores2 = metrics[metric2]

            # Calculate the correlation coefficient
            correlation = np.corrcoef(scores1, scores2)[0, 1]

            legend_name = '_'.join(graph_name.split('_')[1:3])
            plt.scatter(scores1, scores2, color=colors[idx % len(colors)], alpha=0.5, label=f'{legend_name} (corr={correlation:.2f})')

            # Fit a line to the data
            m, b = np.polyfit(scores1, scores2, 1)

            # Plot the line
            plt.plot(scores1, m * np.array(scores1) + b, color=colors[idx % len(colors)])

        plt.xlabel(metric1.capitalize())
        plt.ylabel(metric2.capitalize())
        plt.legend()
        plt.title(f'{metric1.capitalize()} vs. {metric2.capitalize()}')
        if save_as:
            plt.savefig(save_as, dpi=300, format='png')
        plt.show()

def centrality_correlation_plot():
    graph_dict = gi.intra_1mb_graphs()
    centrality_instance = CentralityCorrelation(graph_dict)
    centrality_instance.calculate_centrality_correlation()
    centrality_instance.plot_centrality_correlation(save_as=None, metric1="closeness", metric2="betweenness")
centrality_correlation_plot()



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
