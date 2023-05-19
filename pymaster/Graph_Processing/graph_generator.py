# Import modules
import igraph as ig
from pathlib import Path
import networkx as nx
import pandas as pd
import os as os
import re
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
            return graph
        else:
            return None

    def get_graph_names(self):
        self.c.execute("SELECT name FROM graphs")
        result = self.c.fetchall()
        graphnames = [row[0] for row in result]
        return graphnames

    def get_graph_by_name(self, name):
        self.c.execute("SELECT graph FROM graphs WHERE name = ?", (name,))
        result = self.c.fetchone()
        if result is not None:
            graph = pkl.loads(result[0])
            return graph
        else:
            return None

    def graph_exists(self, name):
        self.c.execute('SELECT name FROM graphs WHERE name = ?', (name,))
        result = self.c.fetchone()
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
        pipeline_condition = parent_folder_components[0].strip(")")
        split_status = parent_folder_components[1]
        norm_status = parent_folder_components[2].split("_")[0]
        condition_slice = parent_folder_components[:3]
        condition = "-".join(condition_slice)
        print(f"Pipeline condition: {pipeline_condition} \nSplit status: {split_status} \nNormalization status: {norm_status} \nCondition: {condition}")

        # Add the parent directory name as a prefix to the graph name
        graph_name = f"{parent_folder}_{graph_name}" if parent_folder else graph_name

        cell_line = graph_name.split("_")[1]
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
        graph["pipeline_condition"] = pipeline_condition
        graph["condition"] = condition
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

def make_all_graphs():
    graph_db_manager = GraphDatabaseManager.from_default_path()
    graph_creator = CreateGraphsFromDirectory("/Users/GBS/Master/HiC-Data/edgelists", graph_db_manager)
    graph_creator.generate_and_store_graphs()
    graph_names = graph_db_manager.get_graph_names()
    print(graph_names)
    return graph_names

def get_graphs():
    graph_db_manager = GraphDatabaseManager.from_default_path()
    graph_dict = graph_db_manager.get_all_graphs()
    return graph_dict



class GraphAttributeInspector:
    def __init__(self, *graph_dicts):
        # Merge dicts into one
        if len(graph_dicts) == 1:
            self.graph_dict = graph_dicts[0]
        else:
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
