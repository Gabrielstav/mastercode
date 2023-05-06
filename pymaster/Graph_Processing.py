# Import modules
import igraph as ig
from pathlib import Path
import re as re
import types as types
import networkx as nx
import pandas as pd
import os as os
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
        files = (file_path for file_path in self.input_directory.rglob(pattern) if file_path.is_file())
        files_list = [str(file_path.resolve()) for file_path in files]
        return files_list

    def from_edgelists(self, cell_lines=None, chromosomes=None, resolutions=None):
        all_files = self.get_files("*")

        for file_path in all_files:
            file_name = os.path.basename(file_path)
            graph_name = re.sub(r"_edgelist\.txt$", "", file_name)

            if cell_lines is not None:
                if not any(cell_line in graph_name for cell_line in cell_lines):
                    continue

            if resolutions is not None:
                graph_name_resolution = re.search(r"_([\d]+)$", graph_name).group(1)
                if graph_name_resolution not in resolutions:
                    continue

            # TODO: Maybe the error is due to this?:

            edges = []
            with open(file_path, "r") as f:
                for line_number, line in enumerate(f, start=1):
                    try:
                        edge = tuple(filter(None, line.strip().split()))
                        edges.append(edge)
                    except UnicodeDecodeError:
                        print(f"Error decoding line {line_number} in file {file_path}")

            # with open(file_path, "r") as f:
            #     edges = [tuple(filter(None, line.strip().split())) for line in f]

            df = pd.DataFrame(edges, columns=['source', 'target'])

            # Filter the edges based on chromosomes
            if chromosomes is not None:
                df = df[df['source'].str.startswith(tuple(chromosomes)) & df['target'].str.startswith(tuple(chromosomes))]

            graph = ig.Graph.TupleList(df.itertuples(index=False), directed=False)
            graph.vs['name'] = [v['name'] for v in graph.vs]  # Assign name attribute to edges
            self.graph_dict[graph_name] = graph

class LargestComponent:
    """
    Class to find the LCC from a given graph object.
    Pass any graph obj dict and return largest conected component.
    """

    def __init__(self, graph_dict_or_function):
        if isinstance(graph_dict_or_function, types.FunctionType):
            self.graph_dict = graph_dict_or_function()
        else:
            self.graph_dict = graph_dict_or_function

    def find_lcc(self):
        largest_component_dict = {}
        for graph_name, graph in self.graph_dict.items():
            largest_component = graph.connected_components().giant()  # or g.clusters().giant() or graph.components().giant()?
            largest_component_dict[graph_name] = largest_component
        return largest_component_dict

    def lcc_membership(self):
        lcc_membership_dict = {}
        for graph_name, graph in self.graph_dict.items():
            components = graph.components()
            sizes = components.sizes()
            largest_component_index = sizes.index(max(sizes))

            lcc_membership_dict[graph_name] = [
                1 if membership == largest_component_index else 0
                for membership in components.membership
            ]
        return lcc_membership_dict

    def print_lcc(self):
        for graph_name, graph in self.find_lcc().items():
            print(f"LCC for: {graph_name} \n size: {graph.vcount()} \n edges: {graph.ecount()}")


# TODO: add filtering of inter graphs, where if inter > only return inter edges

class FilterGraphs:
    """
    Class that filters the graphs based on cell line, chromosome, and resolution.
    """

    def __init__(self, graph_dict):
        self.graph_dict = graph_dict
        self.filtered_chromosomes = set()
        self.filtered_graph_dict = {}

    def filter_graphs(self, cell_lines=None, chromosomes=None, resolutions=None, graph_dict=None):
        if graph_dict is None:
            graph_dict = self.graph_dict

        if cell_lines:
            if isinstance(cell_lines, str):
                cell_lines = [cell_lines]
            cell_lines = set(cell_lines)

        if resolutions:
            if isinstance(resolutions, str):
                resolutions = [resolutions]
            resolutions = set(resolutions)

        if chromosomes:
            if isinstance(chromosomes, str):
                chromosomes = [chromosomes]
            chromosomes = set(chromosomes)
            self.filtered_chromosomes = chromosomes

        for graph_name, graph in graph_dict.items():

            if cell_lines is not None and not any(cell_line in graph_name for cell_line in cell_lines):
                continue

            if resolutions is not None:
                # graph_name_resolution = re.search(r"_(\d+)_", graph_name).group(1)
                print(f"Processing graph name: {graph_name}")
                graph_name_resolution = re.search(r'(\d+)(?=\D*$)', graph_name).group(1)
                if graph_name_resolution not in resolutions:
                    continue

            if chromosomes is not None:
                for chromosome in chromosomes:
                    # Create empty graph with same nodes as input
                    new_graph = ig.Graph()
                    new_graph.add_vertices(graph.vs['name'])
                    new_edges = [(graph.vs[e.source]['name'], graph.vs[e.target]['name']) for e in graph.es
                                 if (graph.vs[e.source]['name'].split(':')[0] == chromosome
                                     or graph.vs[e.target]['name'].split(':')[0] == chromosome)]
                    new_graph.add_edges(new_edges)

                    # Add filtered chromosome(s) to graph name
                    new_graph_name = f"{graph_name}, {chromosome}"
                    self.filtered_graph_dict[new_graph_name] = new_graph
            else:
                self.filtered_graph_dict[graph_name] = graph

        return self.filtered_graph_dict


    def print_filtered_edges(self):
        for graph_name, graph in self.filtered_graph_dict.items():
            print(f"Edges in filtered graph {graph_name}:")
            for edge in graph.es:
                source = graph.vs[edge.source]['name']
                target = graph.vs[edge.target]['name']
                print(f"{source} -- {target}")

# If I need to use networkx for some reason later:
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
