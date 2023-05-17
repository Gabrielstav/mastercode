# Import modules
from Graph_Processing import graph_generator as gg
# from Network_Analysis import Graph_Analysis as ga
from pathlib import Path
import igraph as ig
import sqlite3 as sql
import pickle as pkl

# Module to store instances of graphs from all cell lines across all resolutions.

# Make DB after adding metadata from Pipeline and automating graph creation process, so we can filter on metadata in DB and write class
# that wraps the Filtering and Combination of graphs, so that this class retrives the correct graphs from the DB automatically and applies the filters and combinations,
# then returns the correct graphs, so that only one line of code is needed to get the correct graphs (but not needed now).



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
        self.c.execute('INSERT INTO graphs VALUES (?, ?)', (name, graph_blob))
        self.conn.commit()

    def load_graph(self, name):
        self.c.execute('SELECT graph FROM graphs WHERE name = ?', (name,))
        result = self.c.fetchone()
        if result is not None:
            graph = pkl.loads(result[0])
            return graph
        else:
            return None

class GraphGenerator:

    def __init__(self, root_dir, database):
        self.root_dir = Path(root_dir)
        self.db = database

    def generate_and_store_graphs(self):
        for graph_dir in self.root_dir.iterdir():
            if graph_dir.is_dir():
                graph_files = [file for file in graph_dir.iterdir() if file.is_file() and file.suffix == ".edgelist"]
                for graph_file in graph_files:
                    graph = self.create_graph_from_edgelist(graph_file)
                    graph_name = f"{graph_dir.name}_{graph_file.stem}"
                    self.db.save_graph(graph_name, graph)

    @staticmethod
    def create_graph_from_edgelist(file_path):
        # Instantiate CreateGraphsFromDirectory class
        graph_creator = gg.CreateGraphsFromDirectory(file_path.parent)

        # Use the class method to generate the graph from the edgelist
        graph_creator.from_edgelists()

        # Get the graph from the graph_dict
        # Assuming the keys of graph_dict are graph names
        graph_name = file_path.stem.replace("_edgelist.txt", "")
        graph = graph_creator.graph_dict[graph_name]

        return graph

# Instantiate GraphGenerator and GraphDatabaseManager
graph_generator = gg.CreateGraphsFromDirectory(Path("/Users/GBS/Master/HiC-Data/edgelists"))
graph_db_manager = GraphDatabaseManager.from_default_path()

# Generate and store graphs
graph_generator.generate_and_store_graphs()




def mcf10_inter_chr1():
    all_mcf10 = mcf10_intra_raw_graphs()  # generate and collect graphs
    graph_filter = gg.FilterGraphs(all_mcf10)  # initialize filter with the generated graphs
    filtered_graph = graph_filter.filter_graphs(chromosomes=["chr36"], resolutions=[100], interaction_type="intra")  # Filter graphs but with wrong resolution and chromosome
    # graph_filter.print_graph_attributes()
    # graph_filter.print_filtered_edges()
    return filtered_graph


# MCF-10A

def mcf10_intra_raw_graphs():
    root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/intra/mcf10/raw")
    graph_creator = gg.CreateGraphsFromDirectory(root_dir)
    graph_creator.from_edgelists()
    mcf10_graphs = graph_creator.graph_dict
    # graph_creator.print_edgelist()
    return mcf10_graphs


# mcf10_intra_raw_graphs()

def mcf10_LCCs_intra():
    graphs = mcf10_intra_raw_graphs()
    lcc_instance = gg.LargestComponent(graphs)  # instantiate class
    lccs = lcc_instance.find_lcc()  # get LCCs
    lcc_instance.print_lcc()  # print LCCs
    return lccs




# New stuff

def norm_intra_split_graphs():
    root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/norm_intra_split")
    graph_creator = gg.CreateGraphsFromDirectory(root_dir)
    graph_creator.from_edgelists()
    graphs = graph_creator.graph_dict
    return graphs
# print(norm_intra_split_graphs())

def raw_intra_nosplit_graphs():
    root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/intra_nosplit")
    graph_creator = gg.CreateGraphsFromDirectory(root_dir)
    graph_creator.from_edgelists()
    graphs = graph_creator.graph_dict
    return graphs
# print(raw_intra_nosplit_graphs())

def raw_intra_split_graphs():
    root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/intra_split")
    graph_creator = gg.CreateGraphsFromDirectory(root_dir)
    graph_creator.from_edgelists()
    graphs = graph_creator.graph_dict
    return graphs
# print(raw_intra_split_graphs())

def raw_inter_nosplit_graphs():
    root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/inter_nosplit")
    graph_creator = gg.CreateGraphsFromDirectory(root_dir)
    graph_creator.from_edgelists()
    graphs = graph_creator.graph_dict
    return graphs
# print(raw_inter_nosplit_graphs())

def raw_inter_split_graphs():
    root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/inter_split")
    graph_creator = gg.CreateGraphsFromDirectory(root_dir)
    graph_creator.from_edgelists()
    graphs = graph_creator.graph_dict
    return graphs
# print(raw_inter_split_graphs())







def mcf10_inter_graphs():
    root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/inter/mcf10/nosplit")
    graph_creator = gg.CreateGraphsFromDirectory(root_dir)
    graph_creator.from_edgelists()
    mcf10_graphs = graph_creator.graph_dict
    # graph_creator.print_edgelist()
    return mcf10_graphs


def mcf10_inter_split_graphs():
    root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/inter/mcf10/split")
    graph_creator = gg.CreateGraphsFromDirectory(root_dir)
    graph_creator.from_edgelists()
    mcf10_graphs = graph_creator.graph_dict
    return mcf10_graphs


# MCF-7


def mcf7_intra_raw_graphs():
    root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/intra/mcf7/raw")
    graph_creator = gg.CreateGraphsFromDirectory(root_dir)
    graph_creator.from_edgelists()
    mcf7_graphs = graph_creator.graph_dict
    return mcf7_graphs


def mcf7_intra_norm_graphs():
    root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/intra/mcf7/norm")
    graph_creator = gg.CreateGraphsFromDirectory(root_dir)
    graph_creator.from_edgelists()
    mcf7_graphs = graph_creator.graph_dict
    return mcf7_graphs


def mcf7_inter_graphs():
    root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/inter/mcf7/nosplit")
    graph_creator = gg.CreateGraphsFromDirectory(root_dir)
    graph_creator.from_edgelists()
    mcf7_graphs = graph_creator.graph_dict
    return mcf7_graphs


def mcf7_inter_split_graphs():
    root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/inter/mcf7/split")
    graph_creator = gg.CreateGraphsFromDirectory(root_dir)
    graph_creator.from_edgelists()
    mcf7_graphs = graph_creator.graph_dict
    return mcf7_graphs


# IMR90


def imr90_intra_graphs():
    root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/intra/imr90")
    graph_creator = gg.CreateGraphsFromDirectory(root_dir)
    graph_creator.from_edgelists()
    imr90_graphs = graph_creator.graph_dict
    return imr90_graphs


def imr90_inter_graphs():
    root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/inter/imr90/nosplit")
    graph_creator = gg.CreateGraphsFromDirectory(root_dir)
    graph_creator.from_edgelists()
    imr90_graphs = graph_creator.graph_dict
    return imr90_graphs


def imr90_inter_split_graphs():
    root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/inter/imr90/split")
    graph_creator = gg.CreateGraphsFromDirectory(root_dir)
    graph_creator.from_edgelists()
    imr90_graphs = graph_creator.graph_dict
    return imr90_graphs


# HUVEC


def huvec_intra_graphs():
    root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/intra/huvec")
    graph_creator = gg.CreateGraphsFromDirectory(root_dir)
    graph_creator.from_edgelists()
    huvec_graphs = graph_creator.graph_dict
    return huvec_graphs


def huvec_inter_graphs():
    root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/inter/huvec/nosplit")
    graph_creator = gg.CreateGraphsFromDirectory(root_dir)
    graph_creator.from_edgelists()
    huvec_graphs = graph_creator.graph_dict
    return huvec_graphs


def huvec_inter_split_graphs():
    root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/inter/huvec/split")
    graph_creator = gg.CreateGraphsFromDirectory(root_dir)
    graph_creator.from_edgelists()
    huvec_graphs = graph_creator.graph_dict
    return huvec_graphs


def gsm2824367_intra_graphs():
    root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/intra/gsm2824367")
    graph_creator = gg.CreateGraphsFromDirectory(root_dir)
    graph_creator.from_edgelists()
    gsm2824367_graphs = graph_creator.graph_dict
    return gsm2824367_graphs


# print(gsm2824367_intra_graphs())

def gsm2824367_inter_graphs():
    root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/inter/gsm2824367/nosplit")
    graph_creator = gg.CreateGraphsFromDirectory(root_dir)
    graph_creator.from_edgelists()
    gsm2824367_graphs = graph_creator.graph_dict
    return gsm2824367_graphs


def gsm2824367_inter_split_graphs():
    root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/inter/gsm2824367/split")
    graph_creator = gg.CreateGraphsFromDirectory(root_dir)
    graph_creator.from_edgelists()
    gsm2824367_graphs = graph_creator.graph_dict
    return gsm2824367_graphs


# Filtered graphs here ?


# Intra filtering works
def imr90_chr18():
    all_imr90 = imr90_intra_graphs()
    graph_filter = gg.FilterGraphs(all_imr90)
    filtered_graph = graph_filter.filter_graphs(chromosomes=["chr18"], resolutions=[1000000], interaction_type="intra")
    # graph_filter.print_filtered_edges()
    return filtered_graph


print(imr90_chr18())


def imr90_chr18_components():
    graphs = imr90_chr18()
    components_instance = gg.ConnectedComponents(graphs)
    components = components_instance.find_components()
    # components_instance.print_components()
    return components


# imr90_chr18_components()

# Do centrality measures on connected components,
# DO other stuff on cliques and communities


def imr90_chr1_inter():
    all_imr90 = imr90_inter_graphs()
    graph_filter = gg.FilterGraphs(all_imr90)
    filtered_graph = graph_filter.filter_graphs(chromosomes=["chr1"], resolutions=["1000000"], interaction_type="intra")
    # graph_filter.print_filtered_edges()
    return filtered_graph


# imr90_chr1_inter()

def mcf10_chr1_inter():
    all_mcf10 = mcf10_inter_graphs()
    graph_filter = gg.FilterGraphs(all_mcf10)
    filtered_graph = graph_filter.filter_graphs(chromosomes="chr10", resolutions="1000000", interaction_type="inter")
    # graph_filter.print_filtered_edges()
    return filtered_graph


# Take graphs, filter on them and combine them into one
def mcf10_combined():
    intra_graphs = mcf10_intra_raw_graphs()
    inter_graphs = mcf10_inter_graphs()

    # Filter the intra_graphs
    filter_intra = gg.FilterGraphs(intra_graphs)
    filtered_intra_graphs = filter_intra.filter_graphs(resolutions=[500000], chromosomes=["chr18"], interaction_type="intra")
    # filter_intra.print_filtered_edges()

    # Filter the inter_graphs
    filter_inter = gg.FilterGraphs(inter_graphs)
    filtered_inter_graphs = filter_inter.filter_graphs(resolutions=[1000000], chromosomes=["chr10"], interaction_type="inter")
    # filter_inter.print_filtered_edges()

    # Combine the filtered graphs
    graph_combiner = gg.GraphCombiner([filtered_intra_graphs, filtered_inter_graphs])
    combined_graphs = graph_combiner.combine_matching_graphs()
    # graph_combiner.print_edges(combined_graphs)
    return combined_graphs


# mcf10_combined()

# ig.crea("/Users/GBS/Master/HiC-Data/edgelists/intra/mcf10/raw/mcf10_1000000_edgelist.txt", format="ncol")
# mcf10_tester = ig.Graph.Load("/Users/GBS/Master/HiC-Data/edgelists/intra/mcf10/raw/mcf10_1000000_edgelist.txt", format="ncol", directed=False)

def mcf10_intra_chr3():
    all_mcf10 = mcf10_intra_raw_graphs()
    graph_filter = gg.FilterGraphs(all_mcf10)
    filtered_graph = graph_filter.filter_graphs(chromosomes=["chr3"], resolutions=[1000000], interaction_type="intra")
    # graph_filter.print_filtered_edges()
    return filtered_graph


# mcf10_intra_chr3()

def imr90_new_intra():
    root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/intra/imr90")
    graph_creator = gg.CreateGraphsFromDirectory(root_dir)
    graph_creator.from_edgelists()  # generate graphs from edge lists
    imr90_graphs = graph_creator.graph_dict  # collect generated graphs
    return imr90_graphs


# print(imr90_new_intra())

def imr90_chr1():
    all_imr90 = mcf10_inter_graphs()  # generate and collect graphs
    graph_filter = gg.FilterGraphs(all_imr90)  # initialize filter with the generated graphs
    filtered_graph = graph_filter.filter_graphs(chromosomes=["chr1"], resolutions=[1000000], interaction_type="intra")  # Filter graphs
    # graph_filter.print_graph_attributes()
    graph_filter.print_filtered_edges()
    return filtered_graph


#
# imr90_chr1()


def mcf10_inter_chr1():
    all_mcf10 = mcf10_intra_raw_graphs()  # generate and collect graphs
    graph_filter = gg.FilterGraphs(all_mcf10)  # initialize filter with the generated graphs
    filtered_graph = graph_filter.filter_graphs(chromosomes=["chr36"], resolutions=[100], interaction_type="intra")  # Filter graphs but with wrong resolution and chromosome
    # graph_filter.print_graph_attributes()
    # graph_filter.print_filtered_edges()
    return filtered_graph


# mcf10_inter_chr1()


# Something wrong with components after filtering ?

def imr90_chr18_components():
    graphs = imr90_chr18()
    components_instance = gg.ConnectedComponents(graphs)
    components = components_instance.find_components()
    # components_instance.print_components()
    return components


# imr90_chr18_components()


def mcf10_chr18_components():
    graphs = mcf10_intra_raw_graphs()
    filter_instance = gg.FilterGraphs(graphs)
    filtered_graphs = filter_instance.filter_graphs(chromosomes=["chr18"], resolutions=[1000000], interaction_type="intra")
    components_instance = gg.ConnectedComponents(filtered_graphs)
    components_instance.print_components()
    return components_instance.find_components()


# mcf10_chr18_components()

def mcf10_components():
    graphs = mcf10_intra_raw_graphs()
    components_instance = gg.ConnectedComponents(graphs)
    components = components_instance.find_components()
    components_instance.print_components()
    return components


# mcf10_components()


if __name__ == "__main__":
    pass
