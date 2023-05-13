# Import modules
import igraph as ig
from pathlib import Path
import types as types
import networkx as nx
import pandas as pd
import os as os
import re
import itertools as it

# Set backend of igraph to matplotlib:
ig.config["plotting.backend"] = "matplotlib"


class CreateGraphsFromDirectory:
    """
    Class to create igraph objects from edgelists in a directory.
    """

    def __init__(self, input_directory):
        self.input_directory = Path(input_directory).absolute()
        self.graph_dict = {}

    def get_files(self, pattern):
        files = [file_path for file_path in self.input_directory.rglob(pattern) if file_path.is_file() and not file_path.name.startswith('.')]
        files_list = [str(file_path.resolve()) for file_path in files]
        return files_list

    def create_graphs(self):
        for file in os.listdir(self.input_directory):
            # Create the graph object from the file
            graph = ig.Graph.Read_Ncol(os.path.join(self.input_directory, file))

            # Extract cell line, resolution, and chromosomes
            cell_line = re.search(r"(\w+)(?=_\d)", file).group(1)
            resolution = re.search(r"(\d+)(?=\D*$)", file).group(1)
            chromosomes = sorted(set(v["name"].split(":")[0] for v in graph.vs))

            # Add graph-level attributes
            graph["cell_line"] = cell_line
            graph["resolution"] = int(resolution)
            graph["chromosomes"] = chromosomes

            # Add the graph object to the dictionary
            self.graph_dict[file] = graph

    def from_edgelists(self, cell_lines=None, chromosomes=None, resolutions=None):
        all_files = self.get_files("*")

        for file_path in all_files:
            file_name = os.path.basename(file_path)
            graph_name = re.sub(r"_edgelist\.txt$", "", file_name)

            if cell_lines is not None:
                if not any(cell_line in graph_name for cell_line in cell_lines):
                    print(f"File name '{file_name}' does not contain a valid cell line name.")
                    continue

            if resolutions is not None:
                graph_name_resolution_match = re.search(r"_(\d+)(?=_)", graph_name)
                if graph_name_resolution_match is None:
                    print(f"File name '{file_name}' does not contain a valid resolution.")
                    continue
                graph_name_resolution = int(graph_name_resolution_match.group(1))
            else:
                graph_name_resolution = None

            edges = []
            with open(file_path, "rb") as f:
                for line_number, line in enumerate(f, start=1):
                    try:
                        line = line.decode("utf-8").strip()
                        edge = tuple(filter(None, line.split()))
                        edges.append(edge)
                    except UnicodeDecodeError as e:
                        print(f"Error decoding line {line_number} in file {file_path}: {e}")

            df = pd.DataFrame(edges, columns=['source', 'target'])

            # Filter the edges based on chromosomes
            if chromosomes is not None:
                df = df[df['source'].str.startswith(tuple(chromosomes)) & df['target'].str.startswith(tuple(chromosomes))]

            graph = ig.Graph.TupleList(df.itertuples(index=False), directed=False)

            # Assign name attribute to edges
            graph.vs['name'] = [v['name'] for v in graph.vs]

            # Assign graph-level attributes
            graph["cell_line"] = graph_name.split("_")[0]
            graph["resolution"] = graph_name_resolution
            graph["chromosomes"] = sorted(set(v["name"].split(":")[0] for v in graph.vs))

            self.graph_dict[graph_name] = graph

        return self.graph_dict

    # Old implementation: Works, but does not assign attributes to edges/nodes on creation of graphs
    # def from_edgelists(self, cell_lines=None, chromosomes=None, resolutions=None):
    #     all_files = self.get_files("*")
    #
    #     for file_path in all_files:
    #         file_name = os.path.basename(file_path)
    #         graph_name = re.sub(r"_edgelist\.txt$", "", file_name)
    #
    #         if cell_lines is not None:
    #             if not any(cell_line in graph_name for cell_line in cell_lines):
    #                 continue
    #
    #         if resolutions is not None:
    #             graph_name_resolution = re.search(r"_(\d+)$", graph_name).group(1)
    #             if graph_name_resolution not in resolutions:
    #                 continue
    #
    #         edges = []
    #         with open(file_path, "rb") as f:
    #             for line_number, line in enumerate(f, start=1):
    #                 try:
    #                     line = line.decode("utf-8").strip()
    #                     edge = tuple(filter(None, line.split()))
    #                     edges.append(edge)
    #                 except UnicodeDecodeError as e:
    #                     print(f"Error decoding line {line_number} in file {file_path}: {e}")
    #
    #         df = pd.DataFrame(edges, columns=['source', 'target'])
    #
    #         # Filter the edges based on chromosomes
    #         if chromosomes is not None:
    #             df = df[df['source'].str.startswith(tuple(chromosomes)) & df['target'].str.startswith(tuple(chromosomes))]
    #
    #         graph = ig.Graph.TupleList(df.itertuples(index=False), directed=False)
    #         graph.vs['name'] = [v['name'] for v in graph.vs]  # Assign name attribute to edges
    #         self.graph_dict[graph_name] = graph
    #
    #     return self.graph_dict


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
            subgraphs = components.subgraphs()
            components_dict[graph_name] = subgraphs
        return components_dict

    def print_components(self):
        for graph_name, components in self.find_components().items():
            print(f"Components for: {graph_name}")
            for i, component in enumerate(components):
                print(f"Component {i} nodes: {component.vs['name']}")


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



class FilterGraphs:
    """
    Class that filters the graphs based on cell line, chromosome, and resolution.
    """

    def __init__(self, graph_dict):
        self.graph_dict = graph_dict
        self.filtered_chromosomes = set()
        self.filtered_graph_dict = {}

    def filter_graphs(self, cell_lines=None, chromosomes=None, resolutions=None, interaction_type=None, graph_dict=None, combined=False):

        if graph_dict is None:
            graph_dict = self.graph_dict

        if cell_lines:
            if isinstance(cell_lines, str):
                cell_lines = [cell_lines]
            cell_lines = set(cell_lines)

        if resolutions is not None:
            if isinstance(resolutions, str):
                resolutions = [resolutions]
            resolutions = set(resolutions)

        if chromosomes:
            if isinstance(chromosomes, str):
                chromosomes = [chromosomes]
            chromosomes = set(chromosomes)
            self.filtered_chromosomes = chromosomes

        if interaction_type is not None:
            if interaction_type not in ["intra", "inter"]:
                raise ValueError("Invalid interaction_type: choose 'intra' or 'inter'")
            graph_dict = self.filter_interactions(interaction_type)

        for graph_name, graph in graph_dict.items():

            if cell_lines is not None and not any(cell_line in graph_name for cell_line in cell_lines):
                continue

            if resolutions is not None and not combined:
                print(f"Processing graph name: {graph_name}")
                graph_name_resolution = re.search(r'(\d+)(?=\D*$)', graph_name).group(1)
                if graph_name_resolution not in resolutions:
                    continue

            if chromosomes is not None:
                for chromosome in chromosomes:
                    # Create empty graph with same nodes as input
                    new_graph = ig.Graph()
                    new_graph.add_vertices(graph.vs['name'])
                    new_edges = []
                    for e in graph.es:
                        source_name = graph.vs[e.source]['name']
                        target_name = graph.vs[e.target]['name']
                        source_chromosome = source_name.split(':')[0]
                        target_chromosome = target_name.split(':')[0]

                        if source_chromosome == chromosome or target_chromosome == chromosome:
                            if target_chromosome == chromosome:
                                source_name, target_name = target_name, source_name
                            new_edges.append((source_name, target_name))

                    new_graph.add_edges(new_edges)
                    new_graph.vs['name'] = [v['name'] for v in new_graph.vs]
                    new_graph.es['source'] = [e.source for e in new_graph.es]
                    new_graph.es['target'] = [e.target for e in new_graph.es]

                    # Add filtered chromosome(s) to graph name
                    new_graph_name = f"{graph_name}, {chromosome}"
                    self.filtered_graph_dict[new_graph_name] = new_graph
            else:
                self.filtered_graph_dict[graph_name] = graph

        return self.filtered_graph_dict

    def filter_interactions(self, interaction_type):
        filtered_graph_dict = {}

        for graph_name, graph in self.graph_dict.items():
            # Create an empty graph with the same vertices as the input graph
            new_graph = ig.Graph()
            new_graph.add_vertices(graph.vs["name"])

            if interaction_type == "intra":
                new_edges = [
                    (graph.vs[e.source]["name"], graph.vs[e.target]["name"])
                    for e in graph.es
                    if graph.vs[e.source]["name"].split(":")[0]
                    == graph.vs[e.target]["name"].split(":")[0]
                ]
            elif interaction_type == "inter":
                new_edges = [
                    (graph.vs[e.source]["name"], graph.vs[e.target]["name"])
                    for e in graph.es
                    if graph.vs[e.source]["name"].split(":")[0]
                    != graph.vs[e.target]["name"].split(":")[0]
                ]
            else:
                raise ValueError("Invalid interaction_type: choose 'intra' or 'inter'")

            new_graph.add_edges(new_edges)
            filtered_graph_dict[graph_name] = new_graph

        return filtered_graph_dict

    def print_filtered_edges(self):
        for graph_name, graph in self.filtered_graph_dict.items():
            print(f"Edges in filtered graph {graph_name}:")
            for edge in graph.es:
                source = graph.vs[edge.source]['name']
                target = graph.vs[edge.target]['name']
                print(f"{source} -- {target}")

    def print_graph_attributes(self):
        for graph_dict in self.graph_dict:
            for graph_name, graph in graph_dict.items():

                # Graph-level attributes
                print(f"Graph: {graph_name}")
                graph_attributes = graph.attributes()
                print("Graph attributes:")
                print(graph_attributes)

                # Node-level attributes
                node_attributes = graph.vs.attributes()
                print("Node attributes:")
                print(node_attributes)

                # Edge-level attributes
                edge_attributes = graph.es.attributes()
                print("Edge attributes:")
                print(edge_attributes)
                print("-" * 30)



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
                                          resolutions=filter_params.get("resolutions"))

    @staticmethod
    def combine_graphs(graphs, graph_names):
        combined_graph = graphs[0].copy()

        # Add graph_origin attribute to vertices in the first graph
        for v in combined_graph.vs:
            v["graph_origin"] = {graph_names[0]}

        for graph, graph_name in zip(graphs[1:], graph_names[1:]):
            for v in graph.vs:
                if not combined_graph.vs.select(name_eq=v["name"]):
                    new_vertex = combined_graph.add_vertex(name=v["name"])
                    new_vertex["graph_origin"] = {graph_name}
                else:
                    new_vertex = combined_graph.vs.select(name_eq=v["name"])[0]
                    new_vertex["graph_origin"].add(graph_name)

            # Add edges from the filtered graphs
            for e in graph.es:
                source_vertex = graph.vs[e.source]
                target_vertex = graph.vs[e.target]

                # Add graph_name attribute to the edge
                combined_graph.add_edge(source_vertex["name"], target_vertex["name"], graph_name=graph_name)

        return combined_graph

    @staticmethod
    def matching_graphs(graph_name1, graph_name2):
        graph_name_parts1 = graph_name1.split("_")
        graph_name_parts2 = graph_name2.split("_")
        return graph_name_parts1[0] == graph_name_parts2[0] and graph_name_parts1[1] != graph_name_parts2[1]

    def combine_matching_graphs(self, interaction_type=None):
        combined_graph_dict = {}

        for graph_name1, graph1 in self.graph_dicts[0].items():
            for graph_name2, graph2 in self.graph_dicts[1].items():
                if self.matching_graphs(graph_name1, graph_name2):
                    combined_graph_name = self.create_combined_graph_name([graph_name1, graph_name2])
                    combined_graph = self.combine_graphs([graph1, graph2], [graph_name1, graph_name2])  # , interaction_type=interaction_type)
                    combined_graph_dict[combined_graph_name] = combined_graph

        return combined_graph_dict

    @staticmethod
    def create_combined_graph_name(graph_names):
        combined_graph_name_parts = [graph_name.split("_") for graph_name in graph_names]
        combined_graph_name_parts = list(set(it.chain.from_iterable(combined_graph_name_parts)))
        combined_graph_name = "_".join(sorted(combined_graph_name_parts, key=lambda x: x.strip("chr")))
        return combined_graph_name


    @staticmethod
    def print_edges(graph):
        for graph_name, graph in graph.items():
            print(f"Edges in graph {graph_name}:")
            for edge in graph.es:
                source = graph.vs[edge.source]['name']
                target = graph.vs[edge.target]['name']
                print(f"{source} -- {target}")


class ExportGraphToEdgelist:

    def __init__(self, graph_dict, output_dir):
        self.graph_dict = graph_dict
        self.output_dir = output_dir
        self.export()

    def export(self):
        for graph_name, graph in self.graph_dict.items():
            self.export_edgelist(graph, os.path.join(self.output_dir, f"{graph_name}_edgelist.txt"))

    @staticmethod
    def export_edgelist(graph, file_path):
        with open(file_path, 'w') as f:
            for edge in graph.es:
                source = graph.vs[edge.source]['name']
                target = graph.vs[edge.target]['name']
                f.write(f"{source} {target}\n")


class ConvertIgraphToNetworkx:

    def __init__(self, graph_dict):
        self.graph_dict = graph_dict
        self.nx_graph_dict = {}

    def convert(self):
        for graph_key in self.graph_dict:
            graph = self.graph_dict[graph_key]

            # Create empty NetworkX graph, add nodes and edges
            nx_graph = nx.Graph()
            nx_graph.add_nodes_from(range(graph.vcount()))
            nx_graph.add_edges_from(graph.get_edgelist())

            # Copy node attributes
            for v in graph.vs:
                for attr in v.attributes():
                    nx_graph.nodes[v.index][attr] = v[attr]

            # Copy edge attributes
            for e in graph.es:
                u, v = e.tuple
                for attr in e.attributes():
                    nx_graph.edges[u, v][attr] = e[attr]

            self.nx_graph_dict[graph_key] = nx_graph

        return self.nx_graph_dict

    def __str__(self):
        if not self.nx_graph_dict:
            return "No NetworkX graph generated yet. Run the convert() method first."

        output_str = ""
        for graph_key, nx_graph in self.nx_graph_dict.items():
            nodes = nx_graph.nodes()
            edges = nx_graph.edges()
            output_str += f"Graph: {graph_key}\nNodes: {nodes}\nEdges: {edges}\n\n"

        return output_str
