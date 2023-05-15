# Import modules
import igraph as ig
from pathlib import Path
import networkx as nx
import pandas as pd
import os as os
import re
import itertools as it


class CreateGraphsFromDirectory:
    """
    Class to create igraph objects from edgelists in a directory.
    """

    def __init__(self, input_directory):
        self.input_directory = Path(input_directory).absolute()
        self.graph_dict = {}

    def get_files(self, pattern):
        files = [file_path for file_path in self.input_directory.rglob(pattern) if file_path.is_file() and not file_path.name.startswith('.')]
        return [str(file_path.resolve()) for file_path in files]

    def from_edgelists(self):
        all_files = self.get_files("*")

        for file_path in all_files:
            file_name = os.path.basename(file_path)
            graph_name = re.sub(r"_edgelist\.txt$", "", file_name)
            cell_line = graph_name.split("_")[0]
            graph_name_resolution_match = re.search(r"(\d+)(?=\D*$)", graph_name)
            graph_name_resolution = int(graph_name_resolution_match.group(1)) if graph_name_resolution_match else None

            edges = []
            with open(file_path, "r") as f:
                for line_number, line in enumerate(f, start=1):
                    try:
                        edge = list(filter(None, line.strip().split()))
                        source_chromosome, target_chromosome = edge[0].split(":")[0], edge[1].split(":")[0]
                        interaction_type = "inter" if source_chromosome != target_chromosome else "intra"
                        edges.append((*edge, interaction_type))
                    except Exception as e:
                        print(f"Error processing line {line_number} in file {file_path}: {e}")

            df = pd.DataFrame(edges, columns=['source', 'target', 'interaction_type'])
            df['source_chromosome'] = df['source'].str.split(":").str[0]
            df['target_chromosome'] = df['target'].str.split(":").str[0]

            graph = ig.Graph.TupleList(df.itertuples(index=False), directed=False, edge_attrs=['interaction_type', 'source_chromosome', 'target_chromosome'])

            # Check for isolated nodes
            for node in graph.vs:
                if graph.degree(node) == 0:
                    print(f"Node: {node['name']}, Degree: {graph.degree(node)} has degree 0.")

            # Assign attributes to nodes
            graph.vs['location'] = [name.split(":")[1] for name in graph.vs['name']]
            graph.vs['chromosome'] = [name.split(":")[0] for name in graph.vs['name']]

            # Assign graph-level attributes
            graph["cell_line"] = cell_line
            graph["resolution"] = graph_name_resolution

            self.graph_dict[graph_name] = graph

        return self.graph_dict

    def print_edgelist(self):
        for graph_name, graph in self.graph_dict.items():
            print(f"Graph: {graph_name}")
            for edge in graph.es:
                source = graph.vs[edge.source]['name']
                target = graph.vs[edge.target]['name']
                print(f"{source} -- {target}")


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
    Class that filters the graphs based on cell line, chromosome, resolution and interaction type.
    """

    def __init__(self, graph_dict):
        self.graph_dict = graph_dict

    def filter_graphs(self, cell_lines=None, chromosomes=None, resolutions=None, interaction_type=None):
        filtered_graph_dict = {}

        for graph_name, graph in self.graph_dict.items():
            cell_line = graph["cell_line"]
            resolution = graph["resolution"]

            if cell_lines is not None and cell_line not in cell_lines:
                continue

            if resolutions is not None:
                if resolution is None:
                    print(f"Warning: Graph '{graph_name}' does not have a resolution attribute.")
                    continue
                if resolution not in resolutions:
                    # print(f"No graph with resolution(s) {resolutions} found. The resolution for '{graph_name}' is {resolution}.")
                    continue

            if chromosomes is not None:
                selected_edges = [edge for edge in graph.es if any(graph.vs[node_index]['chromosome'] in chromosomes for node_index in edge.tuple)]

                if interaction_type is not None and 'interaction_type' in graph.es.attributes():
                    selected_edges = [edge for edge in selected_edges if edge['interaction_type'] == interaction_type]

                if not selected_edges:
                    print(f"No edges fulfilling the filtering criteria found in graph '{graph_name}'.")
                    continue
                else:
                    graph = graph.subgraph_edges(selected_edges, delete_vertices=False)

            if interaction_type is not None:
                if 'interaction_type' in graph.es.attributes():
                    selected_edges = [edge for edge in graph.es if edge['interaction_type'] == interaction_type]
                    if selected_edges:
                        graph = graph.subgraph_edges(selected_edges, delete_vertices=False)
                    else:
                        print(f"No edges with interaction type '{interaction_type}' found in graph '{graph_name}'.")

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
        # return graph1["cell_line"] == graph2["cell_line"] and graph1["resolution"] != graph2["resolution"] # make different combinators later based on attributes

    def combine_matching_graphs(self):
        combined_graph_dict = {}

        # Flatten the list of dicts into a list of tuples with graph name and graph
        graph_list = [(graph_name, graph) for graph_dict in self.graph_dicts for graph_name, graph in graph_dict.items()]

        # Create a set to store pairs of graph indices that have been combined
        combined_pairs = set()

        # Compare each graph with every other graph only once
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



class ExportGraphToEdgelist:

    def __init__(self, graph_dict, output_dir):
        self.graph_dict = graph_dict
        self.output_dir = Path(output_dir)

    def export(self):
        for graph_name, graph in self.graph_dict.items():
            self._export_edgelist(graph, self.output_dir / f"{graph_name}_edgelist.txt")

    @staticmethod
    def _export_edgelist(graph, file_path):
        try:
            with open(file_path, 'w') as f:
                for edge in graph.es:
                    source = graph.vs[edge.source]['name']
                    target = graph.vs[edge.target]['name']
                    f.write(f"{source} {target}\n")
        except IOError as e:
            print(f"Failed to write to file: {file_path}. Error: {e}")

    @staticmethod
    def _export_graphml(graph, file_path):
        try:
            graph.write_graphml(file_path)
        except IOError as e:
            print(f"Failed to write to file: {file_path}. Error: {e}")

    @staticmethod
    def _export_gml(graph, file_path):
        try:
            graph.write_gml(file_path)
        except IOError as e:
            print(f"Failed to write to file: {file_path}. Error: {e}")

    @staticmethod
    def _export_pickle(graph, file_path):
        try:
            graph.write_pickle(file_path)
        except IOError as e:
            print(f"Failed to write to file: {file_path}. Error: {e}")


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
