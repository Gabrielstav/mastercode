# Import modules
from Graph_Processing import graph_instances as gi
from Graph_Processing import graph_metrics as gm
from Graph_Processing import graph_generator as gg
import igraph as ig
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
from matplotlib.colors import Normalize
import matplotlib.cm as cm
import numpy as np
from scipy.stats import pearsonr, spearmanr, ks_2samp
import pandas as pd
import pathlib as path
import seaborn as sns
from itertools import combinations
from collections import defaultdict
from scipy.spatial.distance import jensenshannon
from collections import Counter
import random as rd
rd.seed = 42

# TODO: Combine graphs (inter + intra, and mcf7 + mcf10) and visualize shared nodes and edges, and inter vs intra edges ?
# TODO: Make Jaccard edge similarity visualizaiton (Heatmap or barplot?)

class Directories:

    root_path = path.Path("/Users/GBS/Master")
    base_path_figures = path.Path("/Users/GBS/Master/Figures")

    degree_path = base_path_figures / "Degree"
    closeness_path = base_path_figures / "Closeness"
    betweenness_path = base_path_figures / "Betweenness"
    comms_path = base_path_figures / "Communities"
    bed_path = root_path / "HiC-Data/bed"


    if not degree_path.exists():
        degree_path.mkdir(parents=True)

    if not closeness_path.exists():
        closeness_path.mkdir(parents=True)

    if not betweenness_path.exists():
        betweenness_path.mkdir(parents=True)

    if not bed_path.exists():
        bed_path.mkdir(parents=True)


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
            node_size = 20 / graph.vcount()
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
            edge_width = 400 / graph.ecount()
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
            plt.scatter(scores1, scores2, color=colors[idx % len(colors)], alpha=0.35, label=f'{legend_name} (corr={correlation:.2f})')

            # Fit a line to the data
            m, b = np.polyfit(scores1, scores2, 1)

            # Plot the line
            plt.plot(scores1, m * np.array(scores1) + b, color=colors[idx % len(colors)], alpha=0.8, linewidth=2, zorder=10)

        plt.xlabel(metric1.capitalize())
        plt.ylabel(metric2.capitalize())
        plt.legend()
        plt.title(f"{metric1.capitalize()} vs. {metric2.capitalize()}")
        if save_as:
            plt.savefig(save_as, dpi=300, format='png')
        plt.show()


class JaccardSimilarity:

    def __init__(self, graph_dict):
        self.graph_dict = graph_dict
        self.jaccard_similarities = {}

    def calculate_similarity(self, inter_similarity=None):
        node_positions_dict = defaultdict(lambda: defaultdict(set))
        for graph_name, graph in self.graph_dict.items():
            # Extract node location and add it to the set
            for node in graph.vs:
                location = node["name"]
                chromosome = location.split(':')[0]
                node_positions_dict[graph_name][chromosome].add(location)

        if not inter_similarity:
            # Calculate Jaccard similarity between every pair of graphs
            for graph_name1, graph_name2 in combinations(node_positions_dict.keys(), 2):
                for chromosome in node_positions_dict[graph_name1].keys():
                    positions1 = node_positions_dict[graph_name1][chromosome]
                    positions2 = node_positions_dict[graph_name2][chromosome]
                    jaccard_sim = self._jaccard_similarity(positions1, positions2)
                    self.jaccard_similarities[(graph_name1, graph_name2, chromosome)] = jaccard_sim

        if inter_similarity:
            # Calculate Jaccard similarity between every pair of graphs
            for graph_name1, graph_name2 in combinations(node_positions_dict.keys(), 2):
                for chromosome1 in node_positions_dict[graph_name1].keys():
                    for chromosome2 in node_positions_dict[graph_name2].keys():
                        positions1 = node_positions_dict[graph_name1][chromosome1]
                        positions2 = node_positions_dict[graph_name2][chromosome2]
                        jaccard_sim = self._jaccard_similarity(positions1, positions2)
                        self.jaccard_similarities[(graph_name1, graph_name2, chromosome1, chromosome2)] = jaccard_sim

    def average_jaccard_node_similarity(self):
        total_similarity = sum(self.jaccard_similarities.values())
        average_similarity = total_similarity / len(self.jaccard_similarities)
        return average_similarity

    @staticmethod
    def _jaccard_similarity(set1, set2):
        intersection = set1.intersection(set2)
        union = set1.union(set2)
        return len(intersection) / len(union)

    def get_jaccard_similarity(self):
        return self.jaccard_similarities


class JensenShannonSimilarity:

    def __init__(self, graph_dict):
        self.graph_dict = graph_dict
        self.jensen_shannon_similarities = {}

    def calculate_similarity(self, inter_similarity=None):
        node_positions_dict = defaultdict(lambda: defaultdict(list))

        for graph_name, graph in self.graph_dict.items():
            # Extract node location and add it to the list
            for node in graph.vs:
                location = node["name"]
                chromosome = location.split(':')[0]
                node_positions_dict[graph_name][chromosome].append(location)

        if not inter_similarity:
            # Calculate Jensen-Shannon similarity between every pair of graphs
            for graph_name1, graph_name2 in combinations(node_positions_dict.keys(), 2):
                for chromosome in node_positions_dict[graph_name1].keys():
                    positions1 = node_positions_dict[graph_name1][chromosome]
                    positions2 = node_positions_dict[graph_name2][chromosome]
                    js_sim = self._jensen_shannon_similarity(positions1, positions2)
                    self.jensen_shannon_similarities[(graph_name1, graph_name2, chromosome)] = js_sim

        if inter_similarity:
            # Calculate Jaccard similarity between every pair of graphs
            for graph_name1, graph_name2 in combinations(node_positions_dict.keys(), 2):
                for chromosome1 in node_positions_dict[graph_name1].keys():
                    for chromosome2 in node_positions_dict[graph_name2].keys():
                        positions1 = node_positions_dict[graph_name1][chromosome1]
                        positions2 = node_positions_dict[graph_name2][chromosome2]
                        js_sim = self._jensen_shannon_similarity(positions1, positions2)
                        self.jensen_shannon_similarities[(graph_name1, graph_name2, chromosome1, chromosome2)] = js_sim

    def average_jensen_shannon_node_similarity(self):
        total_similarity = sum(self.jensen_shannon_similarities.values())
        average_similarity = total_similarity / len(self.jensen_shannon_similarities)
        return average_similarity

    @staticmethod
    def _jensen_shannon_similarity(list1, list2):
        # Make universal set of all possible node locations to avoid scipy complaining about unequal vector lengths
        universal_set = set(list1).union(set(list2))

        # Count the occurrences of each element in the universal set
        counter1 = Counter(list1)
        counter2 = Counter(list2)

        # Create probability distributions
        prob_dist1 = [counter1[loc] / len(list1) if loc in counter1 else 0 for loc in universal_set]
        prob_dist2 = [counter2[loc] / len(list2) if loc in counter2 else 0 for loc in universal_set]

        # Calculate Jensen-Shannon divergence
        js_divergence = jensenshannon(prob_dist1, prob_dist2)

        # Convert divergence to similarity?
        # js_similarity = 1 - np.sqrt(js_divergence)

        return js_divergence

    def get_jensen_shannon_similarity(self):
        return self.jensen_shannon_similarities


class JaccardHeatmap:

    def __init__(self, jaccard_similarities):
        self.jaccard_similarities = jaccard_similarities
        print(f"Type of self: {type(self)}")
        print(f"Type of self.jaccard_similarities: {type(self.jaccard_similarities)}")

    def jaccard_heatmap(self):
        # Convert to DataFrame
        df = pd.DataFrame.from_dict(self.jaccard_similarities, orient='index', columns=['Jaccard Similarity'])
        df.reset_index(inplace=True)
        df[['Graph1', 'Graph2', 'Chromosome']] = pd.DataFrame(df['index'].tolist(), index=df.index)
        df.drop(columns=['index'], inplace=True)

        # Create a DataFrame for each pair of graphs
        dfs = []
        for (graph1, graph2), group_df in df.groupby(['Graph1', 'Graph2']):
            group_df = group_df[['Chromosome', 'Jaccard Similarity']]
            group_df.set_index('Chromosome', inplace=True)
            group_df.columns = [f'{graph1}_{graph2}']
            dfs.append(group_df)

        # Merge the DataFrames
        df_final = pd.concat(dfs, axis=1)

        # Create heatmap
        plt.figure(figsize=(10, 8))
        sns.heatmap(df_final, annot=True, fmt=".2f", cmap='YlGnBu')
        plt.title('Jaccard Similarity between Chromosomes')
        plt.show()

    def inter_jaccard_heatmap(self):
        # Convert to DataFrame
        df = pd.DataFrame.from_dict(self.jaccard_similarities, orient='index', columns=['Jaccard Similarity'])
        df.reset_index(inplace=True)
        df[['Graph1', 'Graph2', 'Chromosome1', 'Chromosome2']] = pd.DataFrame(df['index'].tolist(), index=df.index)
        df.drop(columns=['index'], inplace=True)

        # Pivot the DataFrame to create a matrix suitable for a heatmap
        df_pivot = df.pivot(index='Chromosome1', columns='Chromosome2', values='Jaccard Similarity')

        # Create heatmap
        plt.figure(figsize=(10, 8))
        sns.heatmap(df_pivot, annot=True, fmt=".2f", cmap='GnBu')
        plt.title('Jaccard Similarity between Chromosomes')
        plt.show()

    def jaccard_barplot(self):
        # Convert to DataFrame
        df = pd.DataFrame.from_dict(self.jaccard_similarities, orient='index', columns=['Jaccard Similarity'])
        df.reset_index(inplace=True)
        df[['Graph1', 'Graph2', 'Chromosome']] = pd.DataFrame(df['index'].tolist(), index=df.index)
        df.drop(columns=['index'], inplace=True)

        # Create color map
        color_map = cm.get_cmap('GnBu')
        colors = color_map(df['Jaccard Similarity'])

        # Create bar plot
        plt.figure(figsize=(10, 8))
        sns.barplot(x='Chromosome', y='Jaccard Similarity', palette=colors, data=df)
        plt.title('Jaccard Similarity between Chromosomes')
        plt.show()


class JensenShannonHeatmap:

    def __init__(self, jensen_shannon_similarities):
        self.jensen_shannon_similarities = jensen_shannon_similarities
        print(f"Type of self: {type(self)}")
        print(f"Type of self.jensen_shannon_similarities: {type(self.jensen_shannon_similarities)}")

    def jensen_shannon_heatmap(self):
        # Convert to DataFrame
        df = pd.DataFrame.from_dict(self.jensen_shannon_similarities, orient='index', columns=['Jensen-Shannon Similarity'])
        df.reset_index(inplace=True)
        df[['Graph1', 'Graph2', 'Chromosome']] = pd.DataFrame(df['index'].tolist(), index=df.index)
        df.drop(columns=['index'], inplace=True)

        # Create a DataFrame for each pair of graphs
        dfs = []
        for (graph1, graph2), group_df in df.groupby(['Graph1', 'Graph2']):
            group_df = group_df[['Chromosome', 'Jensen-Shannon Similarity']]
            group_df.set_index('Chromosome', inplace=True)
            group_df.columns = [f'{graph1}_{graph2}']
            dfs.append(group_df)

        # Merge the DataFrames
        df_final = pd.concat(dfs, axis=1)

        # Create heatmap
        plt.figure(figsize=(10, 8))
        sns.heatmap(df_final, annot=True, fmt=".2f", cmap='YlGnBu')
        plt.title('Jensen-Shannon Similarity between Chromosomes')
        plt.show()

    def inter_jensen_shannon_heatmap(self):
        # Pivot the DataFrame to create a matrix suitable for a heatmap
        df = pd.DataFrame.from_dict(self.jensen_shannon_similarities, orient='index', columns=['Jensen-Shannon Similarity'])
        df.reset_index(inplace=True)
        df[['Graph1', 'Graph2', 'Chromosome1', 'Chromosome2']] = pd.DataFrame(df['index'].tolist(), index=df.index)
        df.drop(columns=['index'], inplace=True)
        df_pivot = df.pivot(index='Chromosome1', columns='Chromosome2', values='Jensen-Shannon Similarity')

        # Create heatmap
        plt.figure(figsize=(10, 8))
        sns.heatmap(df_pivot, annot=True, fmt=".2f", cmap='GnBu')
        plt.title('Jensen-Shannon Similarity between Chromosomes')
        plt.show()

    def jensen_shannon_barplot(self):
        # Create color map
        df = pd.DataFrame.from_dict(self.jensen_shannon_similarities, orient='index', columns=['Jensen-Shannon Similarity'])
        df.reset_index(inplace=True)
        df[['Graph1', 'Graph2', 'Chromosome']] = pd.DataFrame(df['index'].tolist(), index=df.index)
        df.drop(columns=['index'], inplace=True)
        color_map = cm.get_cmap('GnBu')
        colors = color_map(df['Jensen-Shannon Similarity'])

        # Create bar plot
        plt.figure(figsize=(10, 8))
        sns.barplot(x='Chromosome', y='Jensen-Shannon Similarity', palette=colors, data=df)
        plt.title('Jensen-Shannon Similarity between Chromosomes')
        plt.show()


class CalculateCentralitySimilarity:
    def __init__(self, graph_dict, metric):
        self.graph_dict = graph_dict
        self.metric = metric
        self.pearson_correlation = None
        self.kolmogorov_smirnov = None
        self.spearman_correlation = None
        self.jensen_shannon_divergence = None

    def _extract_metric(self, graph):
        if self.metric == 'degree':
            return graph.degree()
        elif self.metric == 'betweenness':
            return graph.betweenness()
        elif self.metric == 'closeness':
            return graph.closeness()
        else:
            raise ValueError(f'Unknown metric: {self.metric}')

    # KS test to compare two cell lines centrality metrics to see if they are drawn from the same distribution
    def calculate_kolmogorov_smirnov(self):
        # calculate Kolmogorov-Smirnov test
        cell_line1, cell_line2 = self.graph_dict.keys()
        metric1 = self._extract_metric(self.graph_dict[cell_line1])
        metric2 = self._extract_metric(self.graph_dict[cell_line2])
        self.kolmogorov_smirnov = ks_2samp(metric1, metric2)[0]

    # Other metrics need to have the same number of nodes, can randomly sample or downscale the larger dataset
    def calculate_pearson_correlation(self):
        # calculate Pearson correlation
        cell_line1, cell_line2 = self.graph_dict.keys()
        metric1 = self._extract_metric(self.graph_dict[cell_line1])
        metric2 = self._extract_metric(self.graph_dict[cell_line2])
        self.pearson_correlation = pearsonr(metric1, metric2)[0]

    def calculate_spearman_correlation(self):
        # calculate Spearman correlation
        cell_line1, cell_line2 = self.graph_dict.keys()
        metric1 = self._extract_metric(self.graph_dict[cell_line1])
        metric2 = self._extract_metric(self.graph_dict[cell_line2])
        self.spearman_correlation = spearmanr(metric1, metric2)[0]

    def calculate_jensen_shannon_divergence(self):
        # calculate Jensen-Shannon divergence
        cell_line1, cell_line2 = self.graph_dict.keys()
        metric1 = self._extract_metric(self.graph_dict[cell_line1])
        metric2 = self._extract_metric(self.graph_dict[cell_line2])
        self.jensen_shannon_divergence = jensenshannon(metric1, metric2)

    def get_kolmogorov_smirnov(self):
        return self.kolmogorov_smirnov

    def get_jensen_shannon_divergence(self):
        return self.jensen_shannon_divergence

    def get_pearson_correlation(self):
        return self.pearson_correlation

    def get_spearman_correlation(self):
        return self.spearman_correlation


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


### Centrality metrics

def plot_degree_network():
    graphs = gi.all_graphs()
    filter_instance = gm.FilterGraphs(graphs)
    filter_instance.filter_graphs(cell_lines=["mcf10"], interaction_type="intra", condition="intra-split-raw", resolutions=[1000000])  # , chromosomes=["chr1"])
    graph_dict = filter_instance.graph_dict
    print(graph_dict)
    lcc_instance = gm.LargestComponent(graph_dict)
    lcc = lcc_instance.find_lcc()
    print(lcc)
    degree_instance = PlotDegreeNetwork(graph_dict)  # Use LCC for withtin-chromosome plots
    degree_instance.calculate_degree()
    degree_instance.normalize_degree()
    degree_instance.plot_graph(save_as=None, normalize=True, color_edges=True, layout="circle")

# plot_degree_network()

def plot_closeness_network():
    graphs = gi.intra_1mb_graphs()
    filter_instance = gm.FilterGraphs(graphs)
    filter_instance.filter_graphs(cell_lines=["mcf7"], interaction_type="intra", chromosomes=["chr1"])
    graph_dict = filter_instance.graph_dict
    print(graph_dict)
    lcc_instance = gm.LargestComponent(graph_dict)
    lcc = lcc_instance.find_lcc()
    print(lcc)
    closeness_instance = cm.PlotClosenessNetwork(lcc)  # Use LCC for withtin-chromosome plots
    closeness_instance.plot_closeness(save_as=None, normalize=False, color_edges=True, layout="fr")


# plot_closeness_network()

def plot_betweenness_network():
    graphs = gi.intra_1mb_graphs()
    filter_instance = gm.FilterGraphs(graphs)
    filter_instance.filter_graphs(cell_lines=["gsm2824367"], interaction_type="intra", condition="intra-split-raw")  # , chromosomes=["chr1"])
    graph_dict = filter_instance.graph_dict
    print(graph_dict)
    # lcc_instance = gm.LargestComponent(graph_dict)
    # lcc = lcc_instance.find_lcc()
    # print(lcc)
    betweenness_instance = cm.PlotBetweennessNetwork(graph_dict)  # Use LCC for withtin-chromosome plots
    betweenness_instance.plot_betweenness(save_as=None, normalize=True, color_edges=False, layout="fr")


# plot_betweenness_network()

def combined_betweenness():
    # Instantiate the FilterGraphs class
    all_graphs = gg.GraphDatabaseManager.from_default_path().get_all_graphs()
    graph_filter_intra = gm.FilterGraphs(all_graphs)
    graph_filter_inter = gm.FilterGraphs(all_graphs)

    # Filter the intra_graphs and inter_graphs
    filtered_intra_graphs = graph_filter_intra.filter_graphs(cell_lines=["mcf10"], resolutions=[1000000], chromosomes=["chr1"], condition="inter-nosplit-raw")
    filtered_inter_graphs = graph_filter_inter.filter_graphs(cell_lines=["mcf10"], resolutions=[1000000], chromosomes=["chr1"], condition="intra-nosplit-raw")

    # Combine the filtered graphs
    graph_combiner = gm.GraphCombiner([filtered_intra_graphs, filtered_inter_graphs])
    combined_graphs = graph_combiner.combine_matching_graphs()
    # graph_combiner.print_edges(combined_graphs)

    # Plot betweenness
    betweenness_instance = PlotBetweennessNetwork(combined_graphs)
    betweenness_instance.plot_betweenness(save_as=None, normalize=True, color_edges=False, layout="fr")


# combined_betweenness()

def plot_degree():
    graph_dict = gi.intra_1mb_graphs()
    degree_instance = cm.PlotDegreeDistribution(graph_dict)
    degree_instance.calculate_degree()
    degree_instance.normalize_degree()
    degree_instance.plot_degree_distribution(save_as=cm.Directories.degree_path / "degree_distribution.png", normalize=False)


# plot_degree()

def closeness_plot():
    graph_dict = gi.intra_1mb_graphs()
    closeness_instance = ClosenessDistribution(graph_dict)
    closeness_instance.calculate_closeness()
    closeness_instance.calculate_closeness_distribution()
    closeness_instance.normalize_closeness()
    # closeness_instance.plot_closeness_distribution(save_as=None)
    closeness_instance.plot_closeness_distribution(save_as=None)


# closeness_plot()

def betweenness_plot():
    graph_dict = gi.intra_1mb_graphs()
    betweenness_instance = BetweennessDistribution(graph_dict)
    betweenness_instance.calculate_betweenness()
    betweenness_instance.calculate_betweenness_distribution()
    betweenness_instance.normalize_betweenness()
    betweenness_instance.plot_betweenness_distribution(save_as=None)


# betweenness_plot()

def centrality_correlation_plot():
    graph_dict = gi.intra_1mb_graphs()
    filter_instance = gm.FilterGraphs(graph_dict)
    filtered = filter_instance.filter_graphs(cell_lines=["mcf10"], resolutions=[1000000], condition="intra-split-raw")
    centrality_instance = CentralityCorrelation(filtered)
    centrality_instance.calculate_centrality_correlation()
    centrality_instance.plot_centrality_correlation(save_as=None, metric1="betweenness", metric2="degree")


# centrality_correlation_plot()

def jaccard_plot():
    # Find graphs and filter on them
    graphs = gi.intra_1mb_graphs()
    filter_instance = gm.FilterGraphs(graphs)
    filter_instance.filter_graphs(cell_lines=["mcf10", "mcf7"], interaction_type="intra", resolutions=[1000000], condition="intra-split-raw")
    graph_dict = filter_instance.graph_dict

    # Calculate Jaccard similarity and plot heatmap
    jaccard_instance = JaccardSimilarity(graph_dict)
    jaccard_instance.calculate_similarity(inter_similarity=False)
    jaccard_similarities = jaccard_instance.get_jaccard_similarity()
    print(jaccard_similarities)

    # Create heatmap
    heatmap_instance = JaccardHeatmap(jaccard_similarities)
    heatmap_instance.jaccard_barplot()


# jaccard_plot()


def plot_js():
    # Find graphs and filter on them
    graphs = gi.all_graphs()
    filter_instance = gm.FilterGraphs(graphs)
    filter_instance.filter_graphs(cell_lines=["mcf10", "mcf7"], interaction_type="intra", resolutions=[1000000], condition="intra-split-raw")
    graph_dict = filter_instance.graph_dict

    # Calculate JS similarity and plot heatmap
    js_inst = JensenShannonSimilarity(graph_dict)
    js_inst.calculate_similarity(inter_similarity=False)
    js_sim = js_inst.get_jensen_shannon_similarity()
    print(js_sim)

    # Create heatmap
    heatmap_instance = JensenShannonHeatmap(js_sim)
    heatmap_instance.jensen_shannon_barplot()


# plot_js()


def compare_cents():
    # Find graphs and filter on them
    graphs = gi.intra_1mb_graphs()
    filter_instance = gm.FilterGraphs(graphs)
    filter_instance.filter_graphs(cell_lines=["mcf10", "mcf7"], interaction_type="intra", resolutions=[1000000], condition="intra-split-raw")
    graph_dict = filter_instance.graph_dict

    # Calculate centrality metrics
    degree_instance = DegreeCentrality(graph_dict)
    degree_instance.calculate_degree()
    degree_instance.normalize_degree()

    closeness_instance = ClosenessCentrality(graph_dict)
    closeness_instance.calculate_closeness()
    closeness_instance.normalize_closeness()

    betweenness_instance = BetweennessCentrality(graph_dict)
    betweenness_instance.calculate_betweenness()
    betweenness_instance.normalize_betweenness()

    # Calculate similarity between degree in two cell lines
    similarity_instance = CalculateCentralitySimilarity(graph_dict, "degree")

    similarity_instance.calculate_kolmogorov_smirnov()
    similarity_instance.get_kolmogorov_smirnov()
    print(similarity_instance.get_kolmogorov_smirnov())


# compare_cents()

def sorter_sort():
    graphs = gi.all_graphs()
    filter_instance = gm.FilterGraphs(graphs)
    filter_instance.filter_graphs(cell_lines=["mcf10"], resolutions=[1000000], interaction_type="intra", condition="intra-split-raw", chromosomes=["chr18"])
    graph_dict = filter_instance.graph_dict

    sorter_instance = gg.Sorter(graph_dict)
    sorted_graph_dict = sorter_instance.sort_graph()
    sorter_instance.print_sorted_edgelist()
    return sorted_graph_dict

# sorted_graph = sorter_sort()

def exporter():
    graphs = gi.all_graphs()
    filter_instance = gm.FilterGraphs(graphs)
    filter_instance.filter_graphs(cell_lines=["mcf10"], resolutions=[50000], interaction_type="intra", condition="intra-split-raw")
    graph_dict = filter_instance.graph_dict

    degree_instance = DegreeCentrality(graph_dict)
    degree_instance.calculate_degree()

    closeness_instance = ClosenessCentrality(graph_dict)
    closeness_instance.calculate_closeness()

    betweenness_instance = BetweennessCentrality(graph_dict)
    betweenness_instance.calculate_betweenness()

    exporter_instance = gg.GraphExporter(graph_dict, path.Path(Directories.bed_path))
    exporter_instance.export_gff3(["degree", "betweenness", "closeness"])

exporter()

if __name__ == "__main__":
    print("main running")
    # plot_betweenness_network()
    # plot_degree()
    # closeness_plot()
    # betweenness_plot()
    # centrality_correlation_plot()
    # plot_jaccard()
    # compare_cents()
