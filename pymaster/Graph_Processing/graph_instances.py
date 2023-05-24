# Import modules
from Graph_Processing import graph_generator as gg
from Graph_Processing import graph_metrics as gm


def all_graphs():
    graph_db_manager = gg.GraphDatabaseManager.from_default_path()
    return graph_db_manager.get_all_graphs()


# Intra graphs 1 MB
def intra_1mb_graphs():
    graph_db_manager = gg.GraphDatabaseManager.from_default_path()
    graph_filter = gm.FilterGraphs(graph_db_manager.get_all_graphs())
    graphs = graph_filter.filter_graphs(cell_lines=["mcf10", "mcf7", "imr90", "gsm2824367", "huvec"], resolutions=[1000000], interaction_type="intra", condition="intra-split-raw")
    # instance_lcc = gm.LargestComponent(mcf10_graphs)
    # lcc = instance_lcc.find_lcc()
    # return lcc
    return graphs

# print(intra_1mb_graph())

def inter_1mb_graphs():
    graph_db_manager = gg.GraphDatabaseManager.from_default_path()
    graph_filter = gm.FilterGraphs(graph_db_manager.get_all_graphs())
    graphs = graph_filter.filter_graphs(cell_lines=["mcf10", "mcf7", "imr90", "gsm2824367", "huvec"], resolutions=[1000000], interaction_type="inter", condition="inter-nosplit-raw")
    return graphs
# print(inter_1mb_graphs())

def intra_500kb_graphs():
    graph_db_manager = gg.GraphDatabaseManager.from_default_path()
    graph_filter = gm.FilterGraphs(graph_db_manager.get_all_graphs())
    graphs = graph_filter.filter_graphs(cell_lines=["mcf10", "mcf7", "imr90", "gsm2824367", "huvec"], resolutions=[500000], interaction_type="intra", condition="intra-split-raw")
    return graphs

def intra_250kb_graphs():
    graph_db_manager = gg.GraphDatabaseManager.from_default_path()
    graph_filter = gm.FilterGraphs(graph_db_manager.get_all_graphs())
    graphs = graph_filter.filter_graphs(cell_lines=["mcf10", "mcf7", "imr90", "gsm2824367", "huvec"], resolutions=[250000], interaction_type="intra", condition="intra-split-raw")
    return graphs

def intra_120kb_graphs():
    graph_db_manager = gg.GraphDatabaseManager.from_default_path()
    graph_filter = gm.FilterGraphs(graph_db_manager.get_all_graphs())
    graphs = graph_filter.filter_graphs(cell_lines=["mcf10", "mcf7", "imr90", "gsm2824367", "huvec"], resolutions=[120000], interaction_type="intra", condition="intra-split-raw")
    return graphs

def intra_80kb_graphs():
    graph_db_manager = gg.GraphDatabaseManager.from_default_path()
    graph_filter = gm.FilterGraphs(graph_db_manager.get_all_graphs())
    graphs = graph_filter.filter_graphs(cell_lines=["mcf10", "mcf7", "imr90", "gsm2824367", "huvec"], resolutions=[80000], interaction_type="intra", condition="intra-split-raw")
    return graphs

def intra_50kb_graphs():
    graph_db_manager = gg.GraphDatabaseManager.from_default_path()
    graph_filter = gm.FilterGraphs(graph_db_manager.get_all_graphs())
    graphs = graph_filter.filter_graphs(cell_lines=["mcf10", "mcf7", "imr90", "gsm2824367", "huvec"], resolutions=[50000], interaction_type="intra", condition="intra-split-raw")
    return graphs

def norm_mcf10_graphs():
    graph_db_manager = gg.GraphDatabaseManager.from_default_path()
    graph_filter = gm.FilterGraphs(graph_db_manager.get_all_graphs())
    mcf10_graphs = graph_filter.filter_graphs(cell_lines=["mcf10"], interaction_type="intra", condition="intra-split-norm")
    return mcf10_graphs

def norm_mcf7_graphs():
    graph_db_manager = gg.GraphDatabaseManager.from_default_path()
    graph_filter = gm.FilterGraphs(graph_db_manager.get_all_graphs())
    mcf7_graphs = graph_filter.filter_graphs(cell_lines=["mcf7"], interaction_type="intra", condition="intra-split-norm")
    return mcf7_graphs








# Make all graphs
def make_all_graphs():
    graph_db_manager = gg.GraphDatabaseManager.from_default_path()
    graph_creator = gg.CreateGraphsFromDirectory("/Users/GBS/Master/HiC-Data/edgelists", graph_db_manager)
    graph_creator.generate_and_store_graphs()
    graph_names = graph_db_manager.get_graph_names()
    print(graph_names)
    return graph_names

# Make isntances and filter graphs


def mcf7_mcf10_intra_1Mb():
    graph_db_manager = gg.GraphDatabaseManager.from_default_path()
    graph_filter = gm.FilterGraphs(graph_db_manager.get_all_graphs())
    graphs = graph_filter.filter_graphs(cell_lines=["mcf10", "mcf7"], resolutions=[1000000], interaction_type="intra", condition="intra-split-raw")
    return graphs


# Or use facade to filter graphs (bad)
def get_mcf10_intra_graphs_facade():
    mcf10_graph = gm.GetGraph("inter-nosplit-raw_mcf10_1000000").Filtered_on(cell_line="mcf10", interaction_type="intra")()
    printer = gm.GraphEdgelistPrint(mcf10_graph)
    printer.print_as_edgelist()
    return mcf10_graph
# mcf7_mcf10_intra_1Mb().degree_distribution()

# Combine graphs
def combined():
    # Instantiate the FilterGraphs class
    all_graphs = gg.GraphDatabaseManager.from_default_path().get_all_graphs()
    graph_filter_intra = gm.FilterGraphs(all_graphs)
    graph_filter_inter = gm.FilterGraphs(all_graphs)

    # Filter the intra_graphs and inter_graphs
    filtered_intra_graphs = graph_filter_intra.filter_graphs(cell_lines=["mcf10"], resolutions=[1000000], chromosomes=["chr1"], condition="inter-nosplit-raw")
    filtered_inter_graphs = graph_filter_inter.filter_graphs(cell_lines=["mcf7"], resolutions=[1000000], chromosomes=["chr1"], condition="inter-nosplit-raw")
    print(filtered_intra_graphs)
    print(filtered_inter_graphs)

    # Combine the filtered graphs
    graph_combiner = gm.GraphCombiner([filtered_intra_graphs, filtered_inter_graphs])
    combined_graphs = graph_combiner.combine_matching_graphs()
    # graph_combiner.print_edges(combined_graphs)
    return combined_graphs



