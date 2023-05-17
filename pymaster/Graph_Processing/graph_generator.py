# Import modules
import igraph as ig
from pathlib import Path
import networkx as nx
import pandas as pd
import os as os
import re
import itertools as it
import sqlite3 as sql
import pickle as pkl

class GraphDatabaseManager:

    @classmethod
    def from_default_path(cls: type = Path):
        database_path = Path("/Users/GBS/Master/HiC-Data/network_db/graphs.db")
        return cls(database_path)

    def __init__(self, db_name=Path("graphs.db")):
        self.conn = sql.connect(db_name.as_posix())
        self.c = self.conn.cursor()
        self.c.execute('''
            CREATE TABLE IF NOT EXISTS graphs (
                name TEXT PRIMARY KEY,
                graph BLOB
            )
        ''')

    def save_graph(self, name, graph):
        graph_blob = pkl.dumps(graph)
        try:
            self.c.execute('INSERT INTO graphs VALUES (?, ?)', (name, graph_blob))
            self.conn.commit()
        except sql.IntegrityError:
            print(f"A graph with the name {name} already exists in the database.")
        print(f"Saving graph {name} with attributes {graph.attributes()}")


    def load_graph(self, name):
        self.c.execute('SELECT graph FROM graphs WHERE name = ?', (name,))
        result = self.c.fetchone()
        if result is not None:
            graph = pkl.loads(result[0])
            print(f"Loaded graph {name} with attributes {graph.attributes()}")
            return graph
        else:
            return None

    def get_graph_names(self):
        self.c.execute("SELECT name FROM graphs")
        result = self.c.fetchall()
        graphnames = [row[0] for row in result]
        print("Fetching all graph names.")  # Debug
        for name in graphnames:
            print(name)
        return graphnames

    def graph_exists(self, name):
        self.c.execute('SELECT name FROM graphs WHERE name = ?', (name,))
        result = self.c.fetchone()
        print(f"Checking if graph {name} exists.")  # Debug
        return result is not None

    def delete_graph(self, graph_name):
        self.c.execute('DELETE FROM graphs WHERE name = ?', (graph_name,))
        self.conn.commit()
        print(f"Graph {graph_name} has been deleted.")

    def get_all_graphs(self):
        graph_dict = {}
        graph_names = self.get_graph_names()
        for graph_name in graph_names:
            graph = self.load_graph(graph_name)
            graph_dict[graph_name] = graph
        return graph_dict

class CreateGraphsFromDirectory:
    """
    Class to create igraph objects from edgelists in a directory.
    """

    def __init__(self, input_directory, database):
        self.input_directory = Path(input_directory).absolute()
        self.graph_dict = {}
        self.db = database

    def get_files(self, pattern):
        files = [file_path for file_path in self.input_directory.rglob(pattern) if file_path.is_file() and not file_path.name.startswith('.')]
        return [str(file_path.resolve()) for file_path in files]

    def from_edgelists(self, file_path):

        file_path = Path(file_path)
        parent_folder = file_path.parent.name
        parent_folder_components = parent_folder.split("-")

        file_name = os.path.basename(file_path)
        graph_name = re.sub(r"_edgelist\.txt$", "", file_name)

        # Extract metadata from the folder name (this should be done in a more robust way, using metadata from pipeline)
        # pipeline_interaction_type = parent_folder_components[0]
        split_status = parent_folder_components[1]
        norm_status = parent_folder_components[2].split("_")[0]

        # Add the parent directory name as a prefix to the graph name
        graph_name = f"{parent_folder}_{graph_name}" if parent_folder else graph_name

        cell_line = graph_name.split("_")[1]
        print(f"Cell line: {cell_line}")  # Debug
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
        graph["split_status"] = split_status
        graph["norm_status"] = norm_status
        # graph["pipeline_interaction_type"] = pipeline_interaction_type

        self.graph_dict[graph_name] = graph

        print(f"Graph created: {graph_name}, Nodes: {graph.vcount()}, Edges: {graph.ecount()}")

        if graph["cell_line"] != cell_line:
            raise ValueError(f"Expected graph['cell_line'] to be {cell_line}, but got {graph['cell_line']}")

        if graph["resolution"] != graph_name_resolution:
            raise ValueError(f"Expected graph['resolution'] to be {graph_name_resolution}, but got {graph['resolution']}")

        return self.graph_dict

    def generate_and_store_graphs(self):
        for graph_file in self.input_directory.rglob("*_edgelist.txt"):
            print(f"Processing file: {graph_file}")
            graph_dir = graph_file.parent
            graph_name = f"{graph_dir.name.lstrip('_')}_{graph_file.stem}"
            graph_name = re.sub(r"_edgelist$", "", graph_name)
            print(f"Graph name: {graph_name}")
            if not self.db.graph_exists(graph_name):
                self.from_edgelists(graph_file)
                self.db.save_graph(graph_name, self.graph_dict[graph_name])

    def print_edgelist(self):
        for graph_name, graph in self.graph_dict.items():
            print(f"Graph: {graph_name}")
            for edge in graph.es:
                source = graph.vs[edge.source]['name']
                target = graph.vs[edge.target]['name']
                print(f"{source} -- {target}")

def make_graphs():
    graph_db_manager = GraphDatabaseManager.from_default_path()
    graph_creator = CreateGraphsFromDirectory("/Users/GBS/Master/HiC-Data/edgelists", graph_db_manager)
    graph_creator.generate_and_store_graphs()
    graph_names = graph_db_manager.get_graph_names()
    print(graph_names)
    return graph_names

make_graphs()

def get_graphs():
    graph_db_manager = GraphDatabaseManager.from_default_path()
    graph_dict = graph_db_manager.get_all_graphs()
    return graph_dict



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

    def filter_graphs(self, cell_lines=None, chromosomes=None, resolutions=None, interaction_type=None, split_statuses=None, norm_statuses=None):
        filtered_graph_dict = {}

        for graph_name, graph in self.graph_dict.items():
            if "cell_line" not in graph.attributes():
                raise KeyError(f"Graph '{graph_name}' does not have a 'cell_line' attribute")
            cell_line = graph["cell_line"]
            resolution = graph["resolution"]
            split_status = graph["split_status"]
            norm_status = graph["norm_status"]

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

def mcf10_1Mb_intra_split():
    graph_db_manager = GraphDatabaseManager.from_default_path()
    graph_dict = graph_db_manager.get_all_graphs()
    graph_filter = FilterGraphs(graph_dict)  # Instantiate the FilterGraphs class
    filtered_graph = graph_filter.filter_graphs(chromosomes=["chr1"], resolutions=[1000000], cell_lines=["mcf10"], interaction_type="intra", split_statuses="split", norm_statuses="raw")
    graph_filter.print_filtered_edges()
    return filtered_graph



mcf10_1Mb_intra_split()

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

class GraphAttributeInspector:
    def __init__(self, *graph_dicts):
        # Merge dirs into one
        self.graph_dict = {k: v for d in graph_dicts for k, v in d.items()}

    def print_graph_attributes(self):
        for graph_name, graph in self.graph_dict.items():
            print(f"\nGraph: {graph_name}")
            print("Graph attributes:", graph.attributes())

    def print_node_attributes(self):
        for graph_name, graph in self.graph_dict.items():
            print(f"\nGraph: {graph_name}")
            print("Node attributes:", graph.vs.attributes())

    def print_edge_attributes(self):
        for graph_name, graph in self.graph_dict.items():
            print(f"\nGraph: {graph_name}")
            print("Edge attributes:", graph.es.attributes())

    def print_edge_list_with_attributes(self):
        for graph_name, graph in self.graph_dict.items():
            print(f"\nGraph: {graph_name}")
            for edge in graph.es:
                source = graph.vs[edge.source]
                target = graph.vs[edge.target]
                print(f"Edge: {source['name']} -- {target['name']}", end='. ')
                print("Source node attributes -->", end=' ')
                for index, attribute in enumerate(graph.vs.attributes()):
                    if index != len(graph.vs.attributes()) - 1:
                        print(f"{attribute}: {source[attribute]}", end=', ')
                    else:
                        print(f"{attribute}: {source[attribute]}", end='. ')
                print("Target node attributes -->", end=' ')
                for index, attribute in enumerate(graph.vs.attributes()):
                    if index != len(graph.vs.attributes()) - 1:
                        print(f"{attribute}: {target[attribute]}", end=', ')
                    else:
                        print(f"{attribute}: {target[attribute]}", end='. ')
                print("Edge attributes -->", end=' ')
                for index, attribute in enumerate(edge.attributes()):
                    if index != len(edge.attributes()) - 1:
                        print(f"{attribute}: {edge[attribute]}", end=', ')
                    else:
                        print(f"{attribute}: {edge[attribute]}", end='. ')
                print()


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
