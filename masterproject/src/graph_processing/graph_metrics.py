# Import modules
from collections import defaultdict
from typing import Union
import itertools as it
from graph_processing import graph_generator as gg


class ConnectedComponents:
    """
    Class to find all connected components from a given graph object.
    """

    def __init__(self, graph_dict):
        self.graph_dict = graph_dict

    def find_components(self):
        components_dict = {}
        for graph_name, graph in self.graph_dict.items():
            components = graph.components()
            for i, component in enumerate(components):
                if len(component) == 1:
                    print(f"Isolated node in {graph_name}: {graph.vs[component[0]]['name']}")
            components_dict[graph_name] = components
        return components_dict

    def print_components(self):
        for graph_name, components in self.find_components().items():
            print(f"Components for: {graph_name}")
            for i, component in enumerate(components):
                # Get the node names for this component from the original graph
                component_nodes = [v['name'] for v in self.graph_dict[graph_name].vs if v.index in component]
                print(f"Component {i} nodes: {component_nodes}")


class LargestComponent:
    """
    Class to find the LCC from a given graph object.
    Pass any graph obj dict and return largest conected component.
    """

    def __init__(self, graph_dict):
        self.graph_dict = graph_dict
        self.largest_component_dict = self.find_lcc()

    def find_lcc(self):
        largest_component_dict = {}
        for graph_name, graph in self.graph_dict.items():
            components = graph.components()
            largest_component = components.giant()
            largest_component_dict[graph_name] = largest_component
        return largest_component_dict

    def find_lcc_per_chromosome(self):
        largest_component_dict = {}
        for graph_name, graph in self.graph_dict.items():
            nodes_by_chromosome = {}
            for idx in range(len(graph.vs)):
                chrom = graph.vs[idx]['name'].split(':')[0]
                if chrom not in nodes_by_chromosome:
                    nodes_by_chromosome[chrom] = []
                nodes_by_chromosome[chrom].append(idx)

            for chrom, nodes in nodes_by_chromosome.items():
                subgraph = graph.subgraph(nodes)
                components = subgraph.components()
                largest_component = components.giant()
                largest_component_dict[(graph_name, chrom)] = largest_component

        return largest_component_dict

    def lcc_membership(self):
        lcc_membership_dict = {}
        for graph_name, graph in self.graph_dict.items():
            components = graph.components()
            sizes = components.sizes()
            largest_component_index = max(range(len(sizes)), key=sizes.__getitem__)

            lcc_membership_dict[graph_name] = [
                1 if membership == largest_component_index else 0
                for membership in components.membership
            ]
        return lcc_membership_dict

    def print_lcc(self):
        for graph_name, graph in self.largest_component_dict.items():
            print(f"LCC for: {graph_name}\nsize: {graph.vcount()}\nedges: {graph.ecount()}")

    def print_edges(self):
        for graph_name, graph in self.graph_dict.items():
            print(f"Edges in filtered graph {graph_name}:")
            for edge in graph.es:
                source = graph.vs[edge.source]['name']
                target = graph.vs[edge.target]['name']
                print(f"{source} -- {target}")


class FilterGraphs:
    """
    Class that filters the graphs based on cell line, chromosome, resolution and interaction type.
    """

    def __init__(self, graph_dict):
        self.graph_dict = graph_dict

    def filter_graphs(self, cell_lines=None, chromosomes=None, resolutions=None, interaction_type=None, split_statuses=None, norm_statuses=None, pipeline_conds=None, condition=None):
        filtered_graph_dict = {}

        for graph_name, graph in self.graph_dict.items():
            if "cell_line" not in graph.attributes():
                raise KeyError(f"Graph '{graph_name}' does not have a 'cell_line' attribute")
            cell_line = graph["cell_line"]
            resolution = graph["resolution"]
            split_status = graph["split_status"]
            norm_status = graph["norm_status"]
            pipeline_cond = graph["pipeline_condition"]
            conditions = graph["condition"]

            if cell_lines is not None and cell_line not in cell_lines:
                continue

            if resolutions is not None:
                if resolution is None:
                    print(f"Warning: Graph '{graph_name}' does not have a resolution attribute.")
                    continue
                if resolution not in resolutions:
                    continue

            if split_statuses is not None and split_status not in split_statuses:
                continue

            if norm_statuses is not None and norm_status not in norm_statuses:
                continue

            if pipeline_conds is not None and pipeline_cond not in pipeline_conds:
                continue

            if condition is not None and conditions is not None and condition not in conditions:
                continue

            if chromosomes is not None:
                selected_edges = [edge for edge in graph.es if any(graph.vs[node_index]['chromosome'] in chromosomes for node_index in edge.tuple)]

                if interaction_type is not None and 'interaction_type' in graph.es.attributes():
                    selected_edges = [edge for edge in selected_edges if edge['interaction_type'] == interaction_type]
                    if not selected_edges:
                        print(f"No edges with interaction type '{interaction_type}' found in graph '{graph_name}'.")

                if not selected_edges:
                    print(f"No edges fulfilling the filtering criteria found in graph '{graph_name}'.")
                    continue
                else:
                    graph = graph.subgraph_edges(selected_edges, delete_vertices=True)

            if interaction_type is not None:
                if 'interaction_type' in graph.es.attributes():
                    selected_edges = [edge for edge in graph.es if edge['interaction_type'] == interaction_type]
                    if selected_edges:
                        graph = graph.subgraph_edges(selected_edges, delete_vertices=True)
                    else:
                        print(f"No edges with interaction type '{interaction_type}' found in graph '{graph_name}'.")
                else:
                    print(f"Warning: Graph '{graph_name}' does not have an 'interaction_type' attribute.")

            filtered_graph_dict[graph_name] = graph

        self.graph_dict = filtered_graph_dict

        return filtered_graph_dict

    def print_graph_attributes(self):
        for graph_name, graph in self.graph_dict.items():
            print(graph.attributes())
            print(graph.es.attributes())
            print(graph.vs.attributes())

    def print_filtered_edges(self):
        for graph_name, graph in self.graph_dict.items():
            print(f"Edges in filtered graph {graph_name}:")
            for edge in graph.es:
                source = graph.vs[edge.source]['name']
                target = graph.vs[edge.target]['name']
                print(f"{source} -- {target}")


class GraphCombiner:

    def __init__(self, graph_dicts, filter_parameters=None):
        if filter_parameters is None:
            self.graph_dicts = graph_dicts
        else:
            self.graph_dicts = [self.filter_graph_dict(graph_dict, filter_params)
                                for graph_dict, filter_params in zip(graph_dicts, filter_parameters)]

    @staticmethod
    def filter_graph_dict(graph_dict, filter_params):
        graph_filter = FilterGraphs(graph_dict)
        return graph_filter.filter_graphs(chromosomes=filter_params.get("chromosomes"),
                                          resolutions=filter_params.get("resolutions"),
                                          interaction_type=filter_params.get("interaction_type"))

    @staticmethod
    def combine_graphs(graphs, graph_names):
        combined_graph = graphs[0].copy()
        vertex_names_to_indices = {v["name"]: v.index for v in combined_graph.vs}

        # Add graph_origin attribute to vertices and edges in the first graph
        for v in combined_graph.vs:
            v["graph_origin"] = {graph_names[0]}
        for e in combined_graph.es:
            e["graph_origin"] = {graph_names[0]}

        for graph, graph_name in zip(graphs[1:], graph_names[1:]):
            for v in graph.vs:
                v_index = vertex_names_to_indices.get(v["name"])
                if v_index is None:
                    new_vertex = combined_graph.add_vertex(name=v["name"])
                    new_vertex["graph_origin"] = {graph_name}
                    vertex_names_to_indices[v["name"]] = new_vertex.index
                else:
                    combined_graph.vs[v_index]["graph_origin"].add(graph_name)

            # Add edges from the filtered graphs
            for e in graph.es:
                source_vertex = graph.vs[e.source]
                target_vertex = graph.vs[e.target]

                # Get edge index if the edge already exists
                edge_index = combined_graph.get_eid(source_vertex["name"], target_vertex["name"], error=False)

                if edge_index == -1:
                    # If edge does not exist, add it with graph_name attribute
                    combined_graph.add_edge(source_vertex["name"], target_vertex["name"], graph_origin={graph_name})
                else:
                    # If edge already exists, add graph_name to its graph_origin attribute
                    combined_graph.es[edge_index]["graph_origin"].add(graph_name)

        return combined_graph

    @staticmethod
    def matching_graphs():
        return True
        # return graph1["cell_line"] == graph2["cell_line"] and graph1["resolution"] != graph2["resolution"]
        # make different combinators later based on attributes, like nested nodes using position etc

    def combine_matching_graphs(self):
        combined_graph_dict = {}

        # Flatten the list of dicts into list of tuples with graph name + graph
        graph_list = [(graph_name, graph) for graph_dict in self.graph_dicts for graph_name, graph in graph_dict.items()]

        # Create set to store pairs of combined graph indices
        combined_pairs = set()

        # Compare each graph with every other graph
        for i, (graph_name1, graph1) in enumerate(graph_list):
            for j, (graph_name2, graph2) in enumerate(graph_list[i + 1:], start=i + 1):  # start parameter adjusts the index
                # Check if the pair has already been combined
                if (i, j) in combined_pairs or (j, i) in combined_pairs:
                    continue

                if self.matching_graphs:
                    combined_graph_name = self.create_combined_graph_name([graph_name1, graph_name2])
                    combined_graph = self.combine_graphs([graph1, graph2], [graph_name1, graph_name2])
                    combined_graph_dict[combined_graph_name] = combined_graph

                # Add the pair to the combined_pairs set
                combined_pairs.add((i, j))

        return combined_graph_dict

    @staticmethod
    def create_combined_graph_name(graph_names):
        combined_graph_name_parts = [graph_name.split("_") for graph_name in graph_names]
        combined_graph_name_parts = list(set(it.chain.from_iterable(combined_graph_name_parts)))
        combined_graph_name = "_".join(sorted(combined_graph_name_parts, key=lambda x: x.strip("chr")))
        combined_graph_name = "combined_" + combined_graph_name  # Add 'combined_' prefix
        return combined_graph_name

    @staticmethod
    def print_edges(combined_graphs):
        for graph_name, graph in combined_graphs.items():
            if graph_name.startswith("combined_"):
                print(f"Edges in graph {graph_name}:")
                if not graph.es:
                    print("No edges in this graph")
                else:
                    for edge in graph.es:
                        source = graph.vs[edge.source]['name']
                        target = graph.vs[edge.target]['name']
                        print(f"{source} -- {target}")

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
        "infomap_communities": lambda g: g.community_infomap(),
        "loouvain_communities": lambda g: g.community_leading_eigenvector(),
        "community": lambda g: g.community_multilevel(),
        "assortativity": lambda g: g.assortativity_degree(),
        "betweenness": lambda g: g.betweenness(),
        "edge_betweenness": lambda g: g.edge_betweenness(),
        "closeness": lambda g: g.closeness(),
        "degree": lambda g: g.degree(),
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
    def print_metrics(cls, graph_dict):
        for graph_name, metrics in graph_dict.items():
            print(f"Metrics for {graph_name}:")
            for metric_name, metric_value in metrics.items():
                print(f"{metric_name}: {metric_value}")
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
            filtered_graphs = FilterGraphs(self.graph_dict).filter_graphs(chromosomes=chromosome)
            for graph_name, graph in filtered_graphs.items():
                lcc_graph = graph.connected_components().giant()
                lcc_graph_size = lcc_graph.vcount()
                parent_graph_size = graph.vcount()
                lcc_ratio = lcc_graph_size / parent_graph_size
                lcc_ratio_dict[chromosome].append(lcc_ratio)

        # Normalize the ratios to sum up to 1
        for chromosome in lcc_ratio_dict:
            ratios = lcc_ratio_dict[chromosome]
            total_ratio = sum(ratios)
            lcc_ratio_dict[chromosome] = [ratio / total_ratio for ratio in ratios]

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

class GraphEdgelistPrint:
    def __init__(self, graph):
        self.graph = graph

    def print_as_edgelist(self):
        print(f"Edge list for graph:")
        for edge in self.graph.es:
            source = self.graph.vs[edge.source]['name']
            target = self.graph.vs[edge.target]['name']
            print(f"{source} -- {target}")
class GetGraph:
    def __init__(self, graph_name):
        self.graph_name = graph_name
        self.filter_conditions = {}

    def Filtered_on(self, **filter_conditions):
        self.filter_conditions = filter_conditions
        return self

    def __call__(self):
        graph_db_manager = gg.GraphDatabaseManager.from_default_path()
        graph_names = graph_db_manager.get_graph_names()

        if self.graph_name not in graph_names:
            raise ValueError(f"Graph with name '{self.graph_name}' does not exist.")

        graph = graph_db_manager.get_graph_by_name(self.graph_name)

        graph_filter = FilterGraphs({self.graph_name: graph})
        filtered_graph_dict = graph_filter.filter_graphs(**self.filter_conditions)
        filtered_graph = filtered_graph_dict.get(self.graph_name)

        if filtered_graph is None:
            print(f"No edges fulfilling the filtering criteria found in graph '{self.graph_name}'.")

        return filtered_graph

    @staticmethod
    def combine_graphs(graphs, graph_names):
        graph_combiner = GraphCombiner()
        combined_graph = graph_combiner.combine_graphs(graphs, graph_names)
        return combined_graph

    def print_as_edgelist(self):
        graph_db_manager = gg.GraphDatabaseManager.from_default_path()
        graph_names = graph_db_manager.get_graph_names()

        if self.graph_name not in graph_names:
            raise ValueError(f"Graph with name '{self.graph_name}' does not exist.")

        graph_filter = FilterGraphs({self.graph_name: None})
        graph_filter.print_filtered_edges()
