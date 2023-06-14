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
import random as rd
import cairo as cr
import cairocffi as crffi
import pathlib as path
import cdlib as cd
from cdlib import algorithms as cdalgs
from cdlib import evaluation as cdeval



# cd.NodeClustering. contains lots of community detection algs




rd.seed = 42

class Directories:
    base_path = path.Path("/Users/GBS/Master/Figures")
    comms_path = base_path / "comms"

    if not comms_path.exists():
        comms_path.mkdir(parents=True)


# TODO:
#  Community Detection:
#     Calculate community structure using Fastgreedy, Louvain, and Infomap methods.
#     Compare the resulting communities across cell lines using the compare_to method (if applicable).
#     ...
#  Community Visualization:
#     Create dendrogram plots of the community clusters.
#     Color nodes in the network plot by the community they belong to.
#     Scale nodes in the community level plot by number of nodes in the community.
#     Label communities in the network plot and dendrograms.
#     ...
#  Community Analysis:
#     Make a histogram of cluster sizes for each cell line.
#     Calculate the partition density for each community.
#     Calculate the mixing parameter between communities.
#     Determine the modularity of each community.
#     Determine the clustering coefficient for each community.
#     Determine the assortativity for each community.
#     ...
#  Community Comparison:
#     Calculate Jaccard similarity between communities (treat communities as sets of nodes).
#     Create a heatmap of the Jaccard index between communities.'

# TODO: Compare everything
#   then make isntances and nice graphs for everything and avoid avoid big render times (if possible: use metrics)
#   ...
#   Lacking:
#   Attribute Analysis: Find out how attribues (centrality) relate to/distribute over communities.
#   External Validation: Lacking external "ground truth" validation of communities, like gene ontology, TADs... Or compartments, which we can integrate.
#   Also, isolate the nodes with highest centrality rank, and the hub nodes, and see where they are. Compare between MCF7 and MCF10A.
#   ...


class CommunityDetection:
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

        return self.graph_dict

    @staticmethod
    def _fastgreedy(graph):
        graph.vs["fg"] = graph.community_fastgreedy().as_clustering().membership


    @staticmethod
    def _louvain(graph):
        graph.vs["lv"] = graph.community_multilevel().membership

    @staticmethod
    def _infomap(graph):
        graph.vs["im"] = graph.community_infomap().membership

    @staticmethod
    def _leiden(graph):
        graph.vs["ld"] = graph.community_leiden().membership

    @staticmethod
    def _label_propagation(graph):
        graph.vs["lp"] = graph.community_label_propagation().membership

    def print_communities(self, method=None):
        if method not in self.methods:
            raise ValueError(f"Invalid method. Choose from: {', '.join(self.methods.keys())}")
        for graph_name, graph in self.graph_dict.items():
            if method not in graph.vs.attributes():
                print(f"No communities calculated for {graph_name} with method {method}")
                continue
            print(graph_name)
            print(graph.vs[method])


class PlotTopology:
    def __init__(self, graph_dict, community_detection):
        self.graph_dict = graph_dict
        self.community_detection = community_detection

    def plot_community(self, method=None, save_as=None):
        if method not in self.community_detection.methods:  # access the methods from the stored instance
            raise ValueError(f"Invalid method. Choose from: {', '.join(self.community_detection.methods.keys())}")

        for graph_name, graph in self.graph_dict.items():
            # Create a color map
            cmap = plt.get_cmap('rainbow', max(graph.vs[method])+1)
            colors = [cmap(i) for i in graph.vs[method]]

            # Set visual style
            fig, ax = plt.subplots()
            node_size = 10 / graph.vcount()
            edge_width = 250 / graph.ecount()
            visual_style = {
                "layout": "kk",
                "vertex_size": node_size,
                "edge_width": edge_width,
                "bbox": (6000, 6000),
                "margin": 200,
                "vertex_color": colors,
                "vertex_label_size": 10,
                "vertex_label_dist": 10,
                "vertex_label_angle": 100,
                "vertex_label_color": "black",
                "vertex_frame_color": colors,  # Remove black outline, use in circle layout
                # "vertex_label": graph.vs["name"],  # Add location as label, can be used in small graphs
            }

            # Plot the graph
            ig.plot(graph, **visual_style, target=ax)
            plt.title(f"{graph_name} - {method}")

            if save_as:
                plt.savefig(Directories.comms_path / f"{graph_name}_{method}.png", dpi=300)

            plt.show()

    def plot_community_level(self):
        """
        Plot the community level of the graph, ie nodes are communities.
        """
        pass


class CommunityComparison:
    def __init__(self, graph_dict1, graph_dict2, method=None, detection_method=None):
        self.graph_dict1 = graph_dict1
        self.graph_dict2 = graph_dict2
        self.method = method
        self.detection_method = detection_method

    def compare_communities(self):
        results = {}
        for graph_name1, graph1 in self.graph_dict1.items():
            for graph_name2, graph2 in self.graph_dict2.items():
                community1 = cdalgs.greedy_modularity(graph1)
                community2 = cdalgs.greedy_modularity(graph2)
                score1 = cd.FuzzyNodeClustering.overlapping_normalized_mutual_information_LFK(community1, community2)
                score2 = cdeval.overlapping_normalized_mutual_information_MGH(community1, community2)

                results[(graph_name1, graph_name2)] = score1, score2
        return results




def compare_coms():
    # Filter on graphs
    graph_dict = gi.all_graphs()
    filter_instance_1 = gm.FilterGraphs(graph_dict)
    filter_instance_2 = gm.FilterGraphs(graph_dict)
    filtered1 = filter_instance_1.filter_graphs(resolutions=[1000000], cell_lines=["mcf10"], condition="intra-split-raw", interaction_type="intra", chromosomes=["chr1"])
    filtered2 = filter_instance_2.filter_graphs(resolutions=[1000000], cell_lines=["gsm2824367"], condition="intra-split-raw", interaction_type="intra", chromosomes=["chr1"])


    # Detect communities
    cd1 = CommunityDetection(filtered1)
    cd2 = CommunityDetection(filtered2)

    filtered1 = cd1.calculate_communities(method="fg")
    filtered2 = cd2.calculate_communities(method="fg")

    # Compare communities
    cc = CommunityComparison(filtered1, filtered2, method="fg", detection_method="fg") # Membership veector must be of equal size, ie same number of nodes, which is not the case.
    comparison_results = cc.compare_communities()

    print(comparison_results)

compare_coms()





class PlotTopologyComparison:
    def __init__(self, graph_dict1, graph_dict2):
        self.graph_dict1 = graph_dict1
        self.graph_dict2 = graph_dict2

    def plot_community(self, method=None, save_as=None):

        # Create a new figure
        fig, ax = plt.subplots()

        for idx, graph_dict in enumerate([self.graph_dict1, self.graph_dict2]):
            for graph_name, graph in graph_dict.items():

                # Create a color map
                cmap = plt.get_cmap('rainbow', max(graph.vs[method]) + 1)
                colors = [cmap(i) for i in graph.vs[method]]

                # Adjust the layout of nodes
                layout = graph.layout("kk")
                for coords in layout:
                    # Shift the x-coordinate of nodes
                    coords[0] += idx * 8
                    # Shift the y-coordinate
                    coords[1] += idx * 2

                # Set visual style
                node_size = 20 / graph.vcount()
                edge_width = 250 / graph.ecount()
                visual_style = {
                    "layout": layout,
                    "vertex_size": node_size,
                    "edge_width": edge_width,
                    "bbox": (6000, 6000),
                    "margin": 200,
                    "vertex_color": colors,
                    "vertex_label_size": 10,
                    "vertex_label_dist": 10,
                    "vertex_label_angle": 100,
                    "vertex_label_color": "black",
                    "vertex_frame_color": colors,
                    #  "vertex_label": graph.vs["name"],
                }

                # Plot the graph
                ig.plot(graph, **visual_style, target=ax)
                # ax.text(idx * 10, 0, graph_name, fontsize=12)

        # Save the figure if required
        if save_as:
            plt.savefig(save_as)

        plt.show()


def mcf10_coms():
    # Filter on graphs
    graph_dict = gi.all_graphs()
    filter_instance = gm.FilterGraphs(graph_dict)
    filtered = filter_instance.filter_graphs(resolutions=[1000000], cell_lines=["mcf10"], condition="intra-split-raw", interaction_type="intra", chromosomes=["chr1"])

    # Detect communities
    cd = CommunityDetection(filtered)
    filtered = cd.calculate_communities(method="fg")  # update filtered with the added community information
    return filtered

def mcf7_coms():
    # Filter on graphs
    graph_dict = gi.all_graphs()
    filter_instance = gm.FilterGraphs(graph_dict)
    filtered = filter_instance.filter_graphs(resolutions=[1000000], cell_lines=["mcf7"], condition="intra-split-raw", interaction_type="intra", chromosomes=["chr1"])

    # Detect communities
    cd = CommunityDetection(filtered)
    filtered = cd.calculate_communities(method="fg")  # update filtered with the added community information
    return filtered

# grapa, grapb = ig.Graph.Erdos_Renyi(100, 0.1), ig.Graph.Erdos_Renyi(100, 0.1)
# dendrogram_a, dendrogram_b = grapa.community_fastgreedy(), grapb.community_fastgreedy()
# a, b = dendrogram_a.as_clustering(), dendrogram_b.as_clustering()
c, d = mcf10_coms(), mcf7_coms()
print(c,d)
# print(ig.compare_communities(c, d, method="nmi"))






class DendrogramPlot:
    def __init__(self, graph_dict, community_detection):
        self.graph_dict = graph_dict
        self.community_detection = community_detection

    def plot_dendrogram(self, method=None):
        if method not in self.community_detection.methods:  # access the methods from the stored instance
            raise ValueError(f"Invalid method. Choose from: {', '.join(self.community_detection.methods.keys())}")
        for graph_name, graph in self.graph_dict.items():
            # Calculate the dendrogram
            dendrogram = graph.community_fastgreedy().as_clustering().as_dendrogram()
            # Plot the dendrogram
            fig, ax = plt.subplots()
            dendrogram.plot(axes=ax)
            ax.set_title(f"{graph_name} - {method}")
            plt.show()

class CommunityMetrics:

    def __init__(self, graph_dict, community_detection):
        self.graph_dict = graph_dict
        self.community_detection = community_detection

    def modularity(self, method=None):
        if method not in self.community_detection.methods:
            raise ValueError(f"Invalid method. Choose from: {', '.join(self.community_detection.methods.keys())}")
        modularity_dict = {}
        for graph_name, graph in self.graph_dict.items():
            community = graph.vs[method]
            modularity = graph.modularity(community)
            modularity_dict[graph_name] = modularity
        return modularity_dict

    def clustering_coefficient(self, method=None):
        if method not in self.community_detection.methods:
            raise ValueError(f"Invalid method. Choose from: {', '.join(self.community_detection.methods.keys())}")
        clustering_coefficient_dict = {}
        for graph_name, graph in self.graph_dict.items():
            clustering_coefficient = graph.transitivity_avglocal_undirected()
            clustering_coefficient_dict[graph_name] = clustering_coefficient
        return clustering_coefficient_dict

    def partition_density(self, method=None):
        if method not in self.community_detection.methods:
            raise ValueError(f"Invalid method. Choose from: {', '.join(self.community_detection.methods.keys())}")
        partition_density_dict = {}
        for graph_name, graph in self.graph_dict.items():
            community = graph.vs[method]
            partition_density = self._calculate_partition_density(graph, community)
            partition_density_dict[graph_name] = partition_density
        return partition_density_dict

    @staticmethod
    def _calculate_partition_density(graph, community):
        densities = []
        for comm in set(community):
            subgraph = graph.subgraph([v for v, c in enumerate(community) if c == comm])
            if subgraph.vcount() < 2:
                density = 0
            else:
                density = subgraph.ecount() / (subgraph.vcount() * (subgraph.vcount() - 1) / 2)
            densities.append(density)
        return sum(densities) / len(densities) if densities else 0

    def conductance(self, method=None):
        if method not in self.community_detection.methods:
            raise ValueError(f"Invalid method. Choose from: {', '.join(self.community_detection.methods.keys())}")
        conductance_dict = {}
        for graph_name, graph in self.graph_dict.items():
            community = graph.vs[method]
            conductance = self._calculate_conductance(graph, community)
            conductance_dict[graph_name] = conductance
        return conductance_dict

    @staticmethod
    def _calculate_conductance(graph, community):
        conductances = []
        for comm in set(community):
            subgraph = graph.subgraph([v for v, c in enumerate(community) if c == comm])
            inner_edges = subgraph.ecount()
            outer_edges = sum(1 for v in subgraph.vs for n in v.neighbors() if n not in subgraph.vs)
            conductance = outer_edges / inner_edges if inner_edges != 0 else 0
            conductances.append(conductance)
        return sum(conductances) / len(conductances)

    def intra_inter_community_edges(self, method=None):
        if method not in self.community_detection.methods:
            raise ValueError(f"Invalid method. Choose from: {', '.join(self.community_detection.methods.keys())}")
        intra_inter_edges_dict = {}
        for graph_name, graph in self.graph_dict.items():
            community = graph.vs[method]
            intra_edges, inter_edges = self._calculate_intra_inter_edges(graph, community)
            intra_inter_edges_dict[graph_name] = (intra_edges, inter_edges)
        return intra_inter_edges_dict

    @staticmethod
    def _calculate_intra_inter_edges(graph, community):
        intra_edges = []
        inter_edges = []
        for comm in set(community):
            subgraph = graph.subgraph([v for v, c in enumerate(community) if c == comm])
            inner_edges = subgraph.ecount()
            outer_edges = sum(1 for v in subgraph.vs for n in v.neighbors() if n not in subgraph.vs)
            intra_edges.append(inner_edges)
            inter_edges.append(outer_edges)
        return intra_edges, inter_edges

class Ideaogram:

    def __init__(self):
        graph_dict = self.graph_dict






def metrics_tester():
    # Filter on graphs
    graph_dict = gi.all_graphs()
    filter_instance = gm.FilterGraphs(graph_dict)
    filtered = filter_instance.filter_graphs(resolutions=[50000], cell_lines=["mcf10"], condition="intra-nosplit-raw", interaction_type="intra", chromosomes=["chr1"])
    print(filtered)

    # Detect communities
    cd = CommunityDetection(filtered)
    filtered = cd.calculate_communities(method="fg")  # update filtered with the added community information

    # Find metrics for communities
    com_met = CommunityMetrics(filtered, cd)
    modularity = com_met.modularity(method="fg")
    clustering_coefficient = com_met.clustering_coefficient(method="fg")
    partition_density = com_met.partition_density(method="fg")
    conductance = com_met.conductance(method="fg")
    intra_inter_edges = com_met.intra_inter_community_edges(method="fg")

    # Print results
    print("Modularity:", modularity)
    print("Clustering Coefficient:", clustering_coefficient)
    print("Partition Density:", partition_density)
    print("Conductance:", conductance)
    print("Intra- and Inter-Community Edges:", intra_inter_edges)

# metrics_tester()








