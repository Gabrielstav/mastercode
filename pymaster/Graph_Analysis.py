# Import modules
from collections import defaultdict
from pathlib import Path
from typing import Union
import Graph_Processing as Gp


class NetworkMetrics:
    """
    Class to calculate network metrics for a given graph object.
    """

    def __init__(self, graph_dict):
        self.graph_dict = graph_dict

    available_metrics = {
        "size": lambda g: g.vcount(),
        "edges": lambda g: g.ecount(),
        "density": lambda g: g.density(),
        "fg_communities": lambda g: g.community_fastgreedy(),
        "community": lambda g: g.community_multilevel(),
        "assortativity": lambda g: g.assortativity_degree(),
        "betweenness": lambda g: g.betweenness(),
        "degree": lambda g: g.degree(),
        "closeness": lambda g: g.closeness(),
        "eigen_centrality": lambda g: g.eigenvector_centrality(),
        "shortest_path_length": lambda g: g.distances(),
        "radius": lambda g: g.radius(),
        "diameter": lambda g: g.diameter(),
        "average_path_length": lambda g: g.average_path_length(),
        "clustering_coefficient": lambda g: g.transitivity_undirected(),
        "jaccard_coefficient": lambda g: g.similarity_jaccard(),
        "pagerank": lambda g: g.pagerank(),
        "fg_modularity": lambda g: g.community_fastgreedy().as_clustering(),
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
    def get_metrics(cls, graph_dict_or_function=None, metrics: Union[None, dict, list] = None):
        """
        Filter metrics from dict, function returning dict or from root directory containing edgelists
        :param cls: instance of class
        :param graph_dict_or_function: input as dictionary or function returning dictionary
        :param metrics: network metrics to calculate
        :return: dict containing graph objects and metrics
        """

        graph_dict = {}

        # If graph_dict_or_function is a function, call it to get the graph_dict
        if callable(graph_dict_or_function):
            graph_dict = graph_dict_or_function()
        # If graph_dict_or_function is a dictionary, use it directly
        elif isinstance(graph_dict_or_function, dict):
            graph_dict = graph_dict_or_function

        # Use metrics
        if metrics is None:
            metrics = cls.available_metrics
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
    def print_metrics(cls, graph_dict, metric_names=None):
        if metric_names is None:
            metric_names = []

        for graph_name, metrics in graph_dict.items():
            print(f"Metrics for {graph_name}:")
            for metric_name, metric_value in metrics.items():
                print(f"  {metric_name}: {metric_value}")
            print()


class LCC_Ratio:

    def __init__(self, graph_dict):
        self.graph_dict = graph_dict

    def calculate_lcc_ratio(self):
        lcc_ratio_dict = {}
        for graph_name, graph in self.graph_dict.items():
            lcc_graph_size = graph.connected_components().giant().vcount()
            parent_graph_size = graph.vcount()
            lcc_ratio_dict[graph_name] = parent_graph_size, lcc_graph_size
        return lcc_ratio_dict

    def calculate_lcc_ratio_per_chromosome(self, chromosomes=None):
        lcc_ratio_dict = defaultdict(list)

        # Filter on chromosome and store unique chromosomes in set
        if chromosomes is None:
            # Filter the graphs by chromosome and store the unique chromosome names in a set
            chromosomes = set()
            for graph_name, graph in self.graph_dict.items():
                for edge in graph.es:
                    source_chromosome = edge.source_vertex['name'].split(':')[0]
                    target_chromosome = edge.target_vertex['name'].split(':')[0]
                    chromosomes.add(source_chromosome)
                    chromosomes.add(target_chromosome)

        for chromosome in chromosomes:
            filtered_graphs = Gp.FilterGraphs(self.graph_dict).filter_graphs(chromosomes=chromosome)
            for graph_name, graph in filtered_graphs.items():
                lcc_graph_size = graph.connected_components().giant().vcount()
                parent_graph_size = graph.vcount()
                lcc_ratio = lcc_graph_size / parent_graph_size
                lcc_ratio_dict[chromosome].append(lcc_ratio)

        return lcc_ratio_dict

    def print_lcc_ratio(self):
        for graph_name, graph in self.calculate_lcc_ratio().items():
            print(f"LCC ratio for: {graph_name} \n size: {graph[0]} \n size: {graph[1]}")


class Filter_Graphs_on_Metrics:

    def __init__(self, graph_dict, degree):
        self.graph_dict = graph_dict
        self.degree = degree

    def filter_degree(self):
        filtered_graph_dict = {}
        for graph_name, graph in self.graph_dict.items():
            nodes_to_keep = [v.index for v in graph.vs if v.degree() >= self.degree]
            subgraph = graph.subgraph(nodes_to_keep)
            filtered_graph_dict[graph_name] = subgraph
        return filtered_graph_dict

    def filter_on_size(self, size):
        filtered_graph_dict = {}
        for graph_name, graph in self.graph_dict.items():
            components = graph.components()
            large_components = [comp for comp in components if len(comp) >= size]
            if large_components:
                subgraph = graph.subgraph(large_components[0])
                filtered_graph_dict[graph_name] = subgraph
        return filtered_graph_dict

# TODO: Make clustering and hub detection class?

# Just for quick testing:
if __name__ == "__main__":

    # # All (intra, raw) graphs
    # def intrachromosomal_graphs():
    #     root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/intra/raw")
    #     graph_creator = Gp.CreateGraphsFromDirectory(root_dir)
    #     graph_creator.from_edgelists()
    #     all_intra_graphss = graph_creator.graph_dict
    #     return all_intra_graphss

    # IMR90
    def imr90_graphs():
        root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/intra/imr90")
        graph_creator = Gp.CreateGraphsFromDirectory(root_dir)
        graph_creator.from_edgelists()
        imr90_graphss = graph_creator.graph_dict
        return imr90_graphss
    print(imr90_graphs())

    def mcf10_graphs():
        root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/intra/mcf10/raw")
        graph_creator = Gp.CreateGraphsFromDirectory(root_dir)
        graph_creator.from_edgelists()
        mcf10_graphss = graph_creator.graph_dict
        return mcf10_graphss
    print(mcf10_graphs())

    def imr90_chr18():
        graph_filter = Gp.FilterGraphs(imr90_graphs())
        filtered_graph = graph_filter.filter_graphs(chromosomes="chr18", resolutions="1000000")
        graph_filter.print_filtered_edges()
        return filtered_graph

    # HUVEC
    def huvec_graphs():
        root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/intra/huvec/edgelists")
        graph_creator = Gp.CreateGraphsFromDirectory(root_dir)
        graph_creator.from_edgelists()
        huvec_graphss = graph_creator.graph_dict
        return huvec_graphss

    def huvec_chr18():
        graph_filter = Gp.FilterGraphs(huvec_graphs())
        filtered_graph = graph_filter.filter_graphs(chromosomes="chr18", resolutions="1000000")
        graph_filter.print_filtered_edges()
        return filtered_graph

    def imr90_degree():
        graphs = Gp.FilterGraphs(imr90_graphs())  # instance
        filtered_graph = graphs.filter_graphs(chromosomes="chr1", resolutions="1000000")  # instance to filter on
        network_metrics = NetworkMetrics(filtered_graph)  # instance to calculate metrics on
        metrics = NetworkMetrics.calculate_metrics(network_metrics, selected_metrics=["degree"])  # instance to calculate metrics on
        return metrics

    print(imr90_degree())


    # # All raw graphs:
    # def mcf710_raw_graphs():
    #     root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/lowres_mcf7_mcf10/raw")
    #     graph_creator = CreateGraphsFromDirectory(root_dir)
    #     graph_creator.from_edgelists()
    #     mcf710_graphs = graph_creator.graph_dict
    #     return mcf710_graphs
    # print(mcf710_raw_graphs())

    # # All norm graphs:
    # def mcf710_norm_graphs():
    #     root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/lowres_mcf7_mcf10/norm")
    #     graph_creator = CreateGraphsFromDirectory(root_dir)
    #     graph_creator.from_edgelists()
    #     mcf7_10_graphs = graph_creator.graph_dict
    #     return mcf7_10_graphs
    #
    # print(mcf710_norm_graphs())

    # # MCF7 raw filtering
    # def mcf7_raw_filtered():
    #     graph_filter = FilterGraphs(mcf710_raw_graphs())
    #     filtered_graphs = graph_filter.filter_graphs(cell_lines=["mcf7"], chromosomes=["chr18"], resolutions=["1000000", "500000"])
    #     # graph_filter.print_filtered_edges()
    #     return filtered_graphs

    # mcf7_raw_filtered()
    #
    # # MCF7 norm filtering
    # def mcf7_norm_filtered():
    #     graph_filter = FilterGraphs(mcf710_norm_graphs())
    #     filtered_graphs = graph_filter.filter_graphs(cell_lines=["mcf7"], chromosomes=["chr18"], resolutions=["1000000", "500000"])
    #     return filtered_graphs
    #
    # mcf7_norm_filtered()
    #
    # # MCF10 raw filtering
    # def mcf10_raw_filtered():
    #     graph_filter = FilterGraphs(mcf710_raw_graphs())
    #     filtered_graphs = graph_filter.filter_graphs(cell_lines=["mcf10"], chromosomes=["chr18"], resolutions=["1000000", "500000"])
    #     return filtered_graphs
    #
    # # MCF10 norm filtering
    # def mcf10_norm_filtered():
    #     graph_filter = FilterGraphs(mcf710_norm_graphs())
    #     filtered_graphs = graph_filter.filter_graphs(cell_lines=["mcf10"], chromosomes=["chr18"], resolutions=["1000000", "500000"])
    #     return filtered_graphs
    #
    # # LCC for MCF10 raw vs norm
    # LargestComponent(mcf10_norm_filtered()).print_lcc()
    # LargestComponent(mcf10_raw_filtered()).print_lcc()
    #
    # # LCC for MCF7 raw vs norm
    # LargestComponent(mcf7_norm_filtered()).print_lcc()
    # LargestComponent(mcf7_raw_filtered()).print_lcc()
    #
    # # Metrics for MCF7 raw filtered:
    # def mcf7_raw_filtered_metrics():
    #     metrics = NetworkMetrics.get_metrics(
    #         graph_dict_or_function=mcf7_raw_filtered(),
    #         metrics=[
    #             "size", "edges", "fg_communities", "betweenness", "closeness", "degree", "average_path_length"
    #         ])
    #     return NetworkMetrics.print_metrics(metrics)
    #
    # mcf7_raw_filtered_metrics()
    #
    #
    # # Metrics for MCF7 norm filtered LCC:
    # def mcf7_norm_filtered_metrics():
    #     metrics = NetworkMetrics.get_metrics(
    #         graph_dict_or_function=LargestComponent(mcf7_norm_filtered()).find_lcc(),
    #         metrics=[
    #             "size", "edges", "fg_communities"
    #         ])
    #     return NetworkMetrics.print_metrics(metrics)
    #
    #
    # mcf7_norm_filtered_metrics()
    #
    # # Metrics for MCF10 raw filtered:
    # def mcf10_raw_filtered_metrics():
    #     metrics = NetworkMetrics.get_metrics(
    #         graph_dict_or_function=mcf10_raw_filtered(),
    #         metrics=[
    #             "size", "edges", "fg_communities", "betweenness", "closeness", "degree", "average_path_length"
    #         ])
    #     return NetworkMetrics.print_metrics(metrics)
    #
    # # Metrics for MCF10 norm filtered:
    # def mcf10_norm_filtered_metrics():
    #     metrics = NetworkMetrics.get_metrics(
    #         graph_dict_or_function=LargestComponent(mcf10_norm_filtered()).find_lcc(),
    #         metrics=[
    #             "size", "edges", "fg_communities"
    #         ])
    #     return NetworkMetrics.print_metrics(metrics)
    #
    # # LCC ratio for MCF7 raw vs norm
    # LCC_ratio(mcf7_raw_filtered()).print_lcc_ratio()
    # LCC_ratio(mcf7_norm_filtered()).print_lcc_ratio()
    # LCC_ratio(mcf10_raw_filtered()).print_lcc_ratio()
    # LCC_ratio(mcf10_norm_filtered()).print_lcc_ratio()
