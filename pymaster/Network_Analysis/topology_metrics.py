# Import modules
import igraph as ig
import matplotlib.pyplot as plt
import numpy as np
from Graph_Processing import graph_metrics as gm
from Graph_Processing import graph_instances as gi
import random as rd
import pathlib as path
import cdlib as cd
from cdlib import algorithms as cdalgs
from cdlib import evaluation as cdeval

# cd.NodeClustering. contains lots of community detection algs




rd.seed = 42

class Directories:
    base_path = path.Path("/Users/GBS/Master/Figures")
    comms_path = path.Path("/Users/GBS/Master/Figures/comms")

    if not comms_path.exists():
        comms_path.mkdir(parents=True)



# TODO:
#     Create dendrogram plots of the community clusters: Compare MCF10-A and MCF-7
#     Compare the largest chromosomes
#     Make similarity NMI plot
#     Annotate communities with compartment A and B


#  Community Analysis???
#     Make a histogram of cluster sizes for each cell line.
#     Calculate the partition density for each community.
#     Calculate the mixing parameter between communities.
#     Determine the modularity of each community.
#     Determine the clustering coefficient for each community.
#     Determine the assortativity for each community.


def chromosome_sort_key(chrom):
    try:
        return int(chrom[3:])
    except ValueError:
        return float("100")


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
            graph_name = graph_name.split("_")[1:2]
            print(graph_name)
            node_size = 30 / graph.vcount()
            edge_width = 60 / graph.ecount()
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
            plt.title(f"{graph_name}")

            if save_as:
                plt.savefig(save_as, dpi=300)

            plt.show()

    def plot_community_level(self):
        """
        Plot the community level of the graph, ie nodes are communities.
        """
        pass

def plot_communities():
    # Filter on graphs
    graph_dict = gi.all_graphs()
    filter_instance = gm.FilterGraphs(graph_dict)
    filtered = filter_instance.filter_graphs(resolutions=[1000000], cell_lines=["mcf7"], condition="intra-split-raw", interaction_type="intra", chromosomes=["chr9"])  # , chromosomes=["chr1"])
    topology_instance = CommunityDetection(filtered)
    topology_instance.calculate_communities(method="fg")
    plot_instance = PlotTopology(filtered, topology_instance)
    plot_instance.plot_community(method="fg", save_as=Directories.comms_path / "chr9_mcf7_fg.png")

# plot_communities()


class CommunitiesPerChromosome:

    def __init__(self, graph_dict):
        self.graph_dict = graph_dict

    @staticmethod
    def _fastgreedy(graph):
        return graph.community_fastgreedy().as_clustering().membership

    def calculate_communities_per_chromosome(self, method=None):
        comm_per_chrom = {}
        avg_comm_size_per_chrom = {}
        for name, graph in self.graph_dict.items():

            # Group vertices by chromosome
            nodes_by_chromosome = {}
            for idx in range(len(graph.vs)):
                chrom = graph.vs[idx]['chromosome']
                if chrom not in nodes_by_chromosome:
                    nodes_by_chromosome[chrom] = []
                nodes_by_chromosome[chrom].append(idx)

            # count unique communities per chrom
            for chrom, nodes in nodes_by_chromosome.items():
                subgraph = graph.subgraph(nodes)
                communities = self._fastgreedy(subgraph)
                for idx, community in enumerate(communities):
                    graph.vs[nodes[idx]][method] = community

                communities = {graph.vs[idx][method] for idx in nodes}
                avg_comm_size = len(nodes) / len(communities)  # average community size
                key = (chrom, graph['cell_line'])
                comm_per_chrom[key] = len(set(communities))
                avg_comm_size_per_chrom[key] = avg_comm_size  # store average community size

        return comm_per_chrom, avg_comm_size_per_chrom


    def plot_communities_per_chromosome(self, method=None, save_as=None):
        comm_per_chrom, _ = self.calculate_communities_per_chromosome(method=method)

        chromosomes = sorted(list(set(chrom for chrom, _ in comm_per_chrom.keys())), key=chromosome_sort_key)
        cell_lines = sorted(list(set(cell_line for _, cell_line in comm_per_chrom.keys())))
        cell_line_colors = ['b', 'r']
        bar_width = 0.35
        for i, cell_line in enumerate(cell_lines):
            values = [comm_per_chrom.get((chrom, cell_line), 0) for chrom in chromosomes]
            plt.bar(np.arange(len(chromosomes)) + i * bar_width, values, color=cell_line_colors[i], width=bar_width, label=cell_line)

        plt.xlabel("Chromosome")
        plt.ylabel("Number of Communities")
        plt.title("Number of Communities per Chromosome")
        plt.xticks(np.arange(len(chromosomes)), [chrom[3:] for chrom in chromosomes])
        plt.legend()

        if save_as:
            plt.savefig(save_as, dpi=300, format='png')
        plt.show()

    def plot_average_community_size_per_chromosome(self, method=None, save_as=None):
        _, avg_comm_size_per_chrom = self.calculate_communities_per_chromosome(method=method)

        chromosomes = sorted(list(set(chrom for chrom, _ in avg_comm_size_per_chrom.keys())), key=chromosome_sort_key)
        cell_lines = sorted(list(set(cell_line for _, cell_line in avg_comm_size_per_chrom.keys())))
        cell_line_colors = ['b', 'r']

        for i, cell_line in enumerate(cell_lines):
            values = [avg_comm_size_per_chrom.get((chrom, cell_line), 0) for chrom in chromosomes]
            # Plot scatter
            plt.scatter(np.arange(len(chromosomes)), values, color=cell_line_colors[i], label=cell_line)
            # Plot polyfit curve
            fit = np.polyfit(np.arange(len(chromosomes)), values, 1)
            fit_fn = np.poly1d(fit)
            plt.plot(np.arange(len(chromosomes)), fit_fn(np.arange(len(chromosomes))), '--', color=cell_line_colors[i])

        plt.xlabel("Chromosome")
        plt.ylabel("Average Community Size")
        plt.title("Average Community Size per Chromosome")
        plt.xticks(np.arange(len(chromosomes)), [chrom[3:] for chrom in chromosomes])
        plt.legend()

        if save_as:
            plt.savefig(save_as, dpi=300, format='png')
        plt.show()


def get_communities_per_chromosome():
    # Get the graphs
    graph_dict = gi.all_graphs()

    # Filter on graphs
    filter_instance = gm.FilterGraphs(graph_dict)
    filtered = filter_instance.filter_graphs(resolutions=[1000000], cell_lines=["mcf7", "mcf10"], condition="intra-split-raw", interaction_type="intra")

    # Calculate communities
    topology_instance = CommunityDetection(filtered)
    topology_instance.calculate_communities(method="fg")

    # Calculate and plot communities per chromosome
    comm_per_chrom_instance = CommunitiesPerChromosome(topology_instance.graph_dict)
    comm_per_chrom_instance.calculate_communities_per_chromosome(method="fg")
    comm_per_chrom_instance.plot_communities_per_chromosome(method="fg", save_as=Directories.comms_path / "mcf7-mcf10_comnum.png")

# get_communities_per_chromosome()

def get_communities_size_per_chromosome():
    # Get the graphs
    graph_dict = gi.all_graphs()

    # Filter on graphs
    filter_instance = gm.FilterGraphs(graph_dict)
    filtered = filter_instance.filter_graphs(resolutions=[1000000], cell_lines=["mcf7", "mcf10"], condition="intra-split-raw", interaction_type="intra")

    # Calculate communities
    topology_instance = CommunityDetection(filtered)
    topology_instance.calculate_communities(method="fg")

    # Calculate and plot communities per chromosome
    comm_per_chrom_instance = CommunitiesPerChromosome(topology_instance.graph_dict)
    comm_per_chrom_instance.calculate_communities_per_chromosome(method="fg")
    comm_per_chrom_instance.plot_average_community_size_per_chromosome(method="fg", save_as=Directories.comms_path / "mcf7-mcf10_avgcomsize.png")

# get_communities_size_per_chromosome()


class NMI_comparison:

    def __init__(self, graph_dict):
        self.graph_dict = graph_dict

    def calculate_communities_per_chromosome(self):
        communities_per_chrom = {}

        for name, graph in self.graph_dict.items():
            nodes_by_chromosome = {}
            for idx in range(len(graph.vs)):
                chrom = graph.vs[idx]['chromosome']
                if chrom not in nodes_by_chromosome:
                    nodes_by_chromosome[chrom] = []
                nodes_by_chromosome[chrom].append(idx)

            for chrom, nodes in nodes_by_chromosome.items():
                subgraph = graph.subgraph(nodes)
                community = cdalgs.greedy_modularity(subgraph)
                key = (chrom, graph['cell_line'])
                print(f"Key: {key}")  # Debugging line
                communities_per_chrom[key] = community
                print(communities_per_chrom)
                print(type(communities_per_chrom))

        return communities_per_chrom

    def calculate_communities_per_chromosome_network(self):
        communities_per_chrom = {}

        for name, graph in self.graph_dict.items():
            nodes_by_chromosome = {}
            for idx in range(len(graph.vs)):  # Iterate over the indices of the nodes
                chrom = graph.vs[idx]['chromosome']
                if chrom not in nodes_by_chromosome:
                    nodes_by_chromosome[chrom] = []
                nodes_by_chromosome[chrom].append(idx)

            for chrom, nodes in nodes_by_chromosome.items():
                subgraph = graph.subgraph(nodes)
                community = cdalgs.louvain(subgraph)
                key = (chrom, name)  # Use the graph name instead of cell line name
                communities_per_chrom[key] = community

        return communities_per_chrom

    def calculate_nmi(self):
        communities_per_chrom = self.calculate_communities_per_chromosome_network()
        nmi_per_chrom = {}

        for chrom in set(chrom for chrom, _ in communities_per_chrom.keys() if _ == "mcf7"):
            communities_mcf10 = communities_per_chrom[(chrom, "mcf10")]
            communities_mcf7 = communities_per_chrom[(chrom, "mcf7")]
            nmi = cdeval.overlapping_normalized_mutual_information_MGH(communities_mcf10, communities_mcf7)
            print(type(communities_mcf7))
            print(communities_mcf7)
            nmi_per_chrom[chrom] = nmi.score

        return nmi_per_chrom

    def calculate_nmi_network(self, communities_per_chrom):
        nmi_per_chrom = {}

        for chrom in set(chrom for chrom, _ in communities_per_chrom.keys() if _ == "intra-split-raw_mcf7_1000000"):
            communities_mcf10 = communities_per_chrom[(chrom, "intra-split-raw_mcf10_1000000")].communities
            communities_mcf7 = communities_per_chrom[(chrom, "intra-split-raw_mcf7_1000000")].communities
            nmi = cdeval.overlapping_normalized_mutual_information_MGH(communities_mcf10, communities_mcf7)
            nmi_per_chrom[chrom] = nmi.score

        return nmi_per_chrom

    def plot_nmi_per_chromosome(self, save_as=None):
        nmi_per_chrom = self.calculate_nmi()

        chromosomes = sorted(list(nmi_per_chrom.keys()), key=chromosome_sort_key)
        nmi_values = [nmi_per_chrom.get(chrom, 0) for chrom in chromosomes]

        plt.bar(np.arange(len(chromosomes)), nmi_values, color='b', width=0.35)

        plt.xlabel("Chromosome")
        plt.ylabel("NMI")
        plt.title("NMI of Communities per Chromosome")
        plt.xticks(np.arange(len(chromosomes)), [chrom[3:] for chrom in chromosomes])

        if save_as:
            plt.savefig(save_as, dpi=300, format='png')
        plt.show()

    def print_same_and_different_nodes(self):
        communities_per_chrom = self.calculate_communities_per_chromosome_network()

        for chrom in set(chrom for chrom, _ in communities_per_chrom.keys() if _ == "intra-split-raw_mcf7_1000000"):
            communities_mcf10 = communities_per_chrom[(chrom, "intra-split-raw_mcf10_1000000")].communities
            communities_mcf7 = communities_per_chrom[(chrom, "intra-split-raw_mcf7_1000000")].communities

            for i, (community_mcf10, community_mcf7) in enumerate(zip(communities_mcf10, communities_mcf7)):
                community_mcf10_set = set(community_mcf10)
                community_mcf7_set = set(community_mcf7)

                same_nodes = community_mcf10_set.intersection(community_mcf7_set)
                different_nodes_mcf10 = community_mcf10_set.difference(community_mcf7_set)
                different_nodes_mcf7 = community_mcf7_set.difference(community_mcf10_set)

                print(f"Community {i} in chromosome {chrom}:")
                print(f"Same nodes: {same_nodes}")
                print(f"Different nodes in mcf10: {different_nodes_mcf10}")
                print(f"Different nodes in mcf7: {different_nodes_mcf7}")
                print("\n")

def get_nmi_per_chromosome():
    # Get the graphs
    graph_dict = gi.all_graphs()

    # Filter on graphs
    filter_instance = gm.FilterGraphs(graph_dict)
    filtered = filter_instance.filter_graphs(resolutions=[1000000], cell_lines=["mcf10", "mcf7"], condition="intra-split-raw", interaction_type="intra", chromosomes=["chr1"])

    # Calculate communities
    topology_instance = NMI_comparison(filtered)
    topology_instance.calculate_communities_per_chromosome()

    # Calculate and plot communities per chromosome
    comm_per_chrom_instance = NMI_comparison(topology_instance.graph_dict)
    comm_per_chrom_instance.calculate_nmi()
    # comm_per_chrom_instance.plot_nmi_per_chromosome(save_as=Directories.comms_path / "mcf7-mcf10_nmi.png")
    comm_per_chrom_instance.print_same_and_different_nodes()

# get_nmi_per_chromosome()

class NMItopology:
    def __init__(self, original_graph_dict, community_detection):
        self.original_graph_dict = original_graph_dict
        self.community_detection = community_detection
        print(self.community_detection)

    def compute_overlap(self):
        same_nodes = set()
        different_nodes_mcf10 = set()
        different_nodes_mcf7 = set()

        for graph_name, graph in self.original_graph_dict.items():
            key = ('chr1', graph_name)  # Construct the key using the chromosome and graph name
            communities = self.community_detection[key].communities
            if isinstance(communities, list):
                # Handle the case where communities is a list of communities
                igraph_communities = []
                for community in communities:
                    igraph_community = [graph.vs[node_index]['cell_line'] for node_index in community]
                    igraph_communities.append(igraph_community)
            elif isinstance(communities, dict):
                # Handle the case where communities is a dictionary mapping nodes to community labels
                igraph_communities = []
                node_community_map = communities
                for community_label in set(node_community_map.values()):
                    igraph_community = [node['name'] for node in graph.vs if node_community_map[node.index] == community_label]
                    igraph_communities.append(igraph_community)

            # Update the community labels to be consistent with igraph's community structure
            igraph_communities = []
            for community in communities:
                igraph_community = [graph.vs[node_index]['name'] for node_index in community]
                igraph_communities.append(igraph_community)

            if graph_name == 'intra-split-raw_mcf10_1000000':
                different_nodes_mcf10.update(set(graph.vs['name']).difference(*igraph_communities))
            elif graph_name == 'intra-split-raw_mcf7_1000000':
                different_nodes_mcf7.update(set(graph.vs['name']).difference(*igraph_communities))

            same_nodes = set(graph.vs['name']).difference(different_nodes_mcf10, different_nodes_mcf7)

        return same_nodes, different_nodes_mcf10, different_nodes_mcf7

    def plot_community(self, same_nodes=None, different_nodes_mcf10=None, different_nodes_mcf7=None, save_as=None):
        # if method not in self.community_detection.methods:  # access the methods from the stored instance
        #     raise ValueError(f"Invalid method. Choose from: {', '.join(self.community_detection.methods.keys())}")

        for graph_name, graph in self.original_graph_dict.items():
            # Create color map for each set of nodes
            node_colors = []
            for node in graph.vs:
                if node['name'] in same_nodes:
                    node_colors.append('green')  # green for overlapping nodes
                elif node['name'] in different_nodes_mcf10:
                    node_colors.append('blue')  # blue for unique nodes in mcf10
                elif node['name'] in different_nodes_mcf7:
                    node_colors.append('red')  # red for unique nodes in mcf7

            # Set visual style
            fig, ax = plt.subplots()
            graph_name = graph_name.split("_")[1:2]
            node_size = 30 / graph.vcount()
            edge_width = 60 / graph.ecount()
            visual_style = {
                "layout": "kk",
                "vertex_size": node_size,
                "edge_width": edge_width,
                "bbox": (2000, 2000),
                "margin": 200,
                "vertex_color": node_colors,
                "vertex_label_size": 10,
                "vertex_label_dist": 10,
                "vertex_label_angle": 100,
                "vertex_label_color": "black",
                "vertex_frame_color": node_colors,
            }

            # Plot the graph
            ig.plot(graph, **visual_style, target=ax)
            plt.title(f"{graph_name}")

            if save_as:
                plt.savefig(save_as, dpi=300)

            plt.show()

def nmi_network():
    # Get the graphs
    graph_dict = gi.all_graphs()

    # Filter on graphs
    filter_instance = gm.FilterGraphs(graph_dict)
    filtered = filter_instance.filter_graphs(resolutions=[1000000], cell_lines=["mcf10", "mcf7"], condition="intra-split-raw", interaction_type="intra", chromosomes=["chr1"])

    # Calculate communities
    topology_instance = NMI_comparison(filtered)
    communities_per_chrom = topology_instance.calculate_communities_per_chromosome_network()

    # Calculate and plot communities per chromosome
    comm_per_chrom_instance = NMI_comparison(communities_per_chrom)
    nmi_per_chrom = comm_per_chrom_instance.calculate_nmi_network(communities_per_chrom)
    # comm_per_chrom_instance.print_same_and_different_nodes()

    # Plot the communities
    nmi_topology_instance = NMItopology(filtered, communities_per_chrom)
    same_nodes, different_nodes_mcf10, different_nodes_mcf7 = nmi_topology_instance.compute_overlap()

    # Iterate over each graph and plot
    for graph_name, graph in topology_instance.graph_dict.items():
        nmi_topology_instance.plot_community(graph, same_nodes, different_nodes_mcf10, different_nodes_mcf7, save_as=None)  # Directories.comms_path / f"{graph_name}_nmi.png")

# nmi_network()

def nmi_net2():
    graph_dict = gi.all_graphs()

    filter_instance = gm.FilterGraphs(graph_dict)
    filtered = filter_instance.filter_graphs(resolutions=[1000000], cell_lines=["mcf10", "mcf7"], condition="intra-split-raw", interaction_type="intra", chromosomes=["chr1"])

    nmi_comparison_instance = NMI_comparison(filtered)
    communities_per_chrom = nmi_comparison_instance.calculate_communities_per_chromosome_network()

    nmi_topology_instance = NMItopology(filtered, communities_per_chrom)  # Provide the community_detection argument
    same_nodes, different_nodes_mcf10, different_nodes_mcf7 = nmi_topology_instance.compute_overlap()
    nmi_topology_instance.plot_community(same_nodes, different_nodes_mcf10, different_nodes_mcf7, save_as=None)

# nmi_net2()
# WElL none of this worked ...



class Dendogram:

    def __init__(self, graph_dict):
        self.graph_dict = graph_dict

    def plot_dendogram(self):

        for graph_name, graph in self.graph_dict.items():
            dendogram = ig.Graph.community_fastgreedy(graph)
            clusters = dendogram.as_clustering()
            fig, ax = plt.subplots()
            ig.plot(clusters, target=ax)
            plt.show()
def dendogram():
    graph_dict = gi.all_graphs()

    filter_instance = gm.FilterGraphs(graph_dict)
    filtered = filter_instance.filter_graphs(resolutions=[1000000], cell_lines=["mcf10", "mcf7"], condition="intra-split-raw", interaction_type="intra", chromosomes=["chr1"])

    dendogram_instance = Dendogram(filtered)
    dendogram_instance.plot_dendogram()

# dendogram()




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
                score3 = cdeval.partition_closeness_simple(community1, community2)
                score4 = cd.FuzzyNodeClustering.f1(community1, community2)

                results[(graph_name1, graph_name2)] = score1, score2, score3, score4
        return results


    def compare_node_communities(self):
        # Make method to use the fg community membership mapped to nodes, and find differences between cell lines
        pass


def compare_coms():
    # Filter on graphs
    graph_dict = gi.all_graphs()
    filter_instance_1 = gm.FilterGraphs(graph_dict)
    filter_instance_2 = gm.FilterGraphs(graph_dict)
    filtered1 = filter_instance_1.filter_graphs(resolutions=[1000000], cell_lines=["mcf10"], condition="intra-split-raw", interaction_type="intra", chromosomes=["chr2"])
    filtered2 = filter_instance_2.filter_graphs(resolutions=[1000000], cell_lines=["mcf7"], condition="intra-split-raw", interaction_type="intra", chromosomes=["chr2"])

    # {('intra-split-raw_mcf10_50000', 'intra-split-raw_mcf7_50000'): (MatchingResult(score=0.09531429795246016, std=None), MatchingResult(score=0.06330616073371602, std=None),
    # MatchingResult(score=0.18415529905561384, std=None), MatchingResult(score=0.38147117296222666, std=0.23803499472082973))}

    # {('intra-split-raw_mcf10_1000000', 'intra-split-raw_mcf7_1000000'): (MatchingResult(score=0.23712347613448248, std=None), MatchingResult(score=0.1679255378736726, std=None),
    # MatchingResult(score=0.0784313725490196, std=None), MatchingResult(score=0.49666666666666676, std=0.15684387141358122))}

    # Detect communities
    cd1 = CommunityDetection(filtered1)
    cd2 = CommunityDetection(filtered2)

    filtered1 = cd1.calculate_communities(method="fg")
    filtered2 = cd2.calculate_communities(method="fg")

    # Compare communities
    cc = CommunityComparison(filtered1, filtered2, method="fg", detection_method="fg") # Membership veector must be of equal size, ie same number of nodes, which is not the case.
    comparison_results = cc.compare_communities()

    print(comparison_results)

# compare_coms()





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
                    "bbox": (2000, 2000),
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
    filtered = filter_instance.filter_graphs(resolutions=[1000000], cell_lines=["mcf10"], condition="intra-split-raw", interaction_type="intra")  # , chromosomes=["chr1"])

    # Detect communities
    cd = CommunityDetection(filtered)
    filtered = cd.calculate_communities(method="fg")  # update filtered with the added community information
    return filtered

# mcf10_coms()

def mcf7_coms():
    # Filter on graphs
    graph_dict = gi.all_graphs()
    filter_instance = gm.FilterGraphs(graph_dict)
    filtered = filter_instance.filter_graphs(resolutions=[1000000], cell_lines=["mcf7"], condition="intra-split-raw", interaction_type="intra", chromosomes=["chr1"])

    # Detect communities
    com = CommunityDetection(filtered)
    filtered = com.calculate_communities(method="fg")  # update filtered with the added community information
    return filtered

# grapa, grapb = ig.Graph.Erdos_Renyi(100, 0.1), ig.Graph.Erdos_Renyi(100, 0.1)
# dendrogram_a, dendrogram_b = grapa.community_fastgreedy(), grapb.community_fastgreedy()
# a, b = dendrogram_a.as_clustering(), dendrogram_b.as_clustering()
# c, d = mcf10_coms(), mcf7_coms()
# print(c,d)
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



def metrics_tester():
    # Filter on graphs
    graph_dict = gi.all_graphs()
    filter_instance = gm.FilterGraphs(graph_dict)
    filtered = filter_instance.filter_graphs(resolutions=[1000000], cell_lines=["mcf7"], condition="intra-nosplit-raw", interaction_type="intra")  # , chromosomes=["chr1"])
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







if __name__ == "__main__":
    print("DONE: TOPOLGGY METRICS")






