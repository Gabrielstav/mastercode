# Import modules
from Graph_Processing import graph_generator as gg
# from Network_Analysis import Graph_Analysis as ga
from pathlib import Path

# Module to store instances of graphs from all cell lines across all resolutions.


# MCF-10A

def mcf10_intra_raw_graphs():
    root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/intra/mcf10/raw")
    graph_creator = gg.CreateGraphsFromDirectory(root_dir)
    graph_creator.from_edgelists()
    mcf10_graphs = graph_creator.graph_dict
    return mcf10_graphs
# print(mcf10_intra_raw_graphs())

def mcf10_LCCs_intra():
    graphs = mcf10_intra_raw_graphs()
    lcc_instance = gg.LargestComponent(graphs)  # instantiate class
    lccs = lcc_instance.find_lcc()  # get LCCs
    return lccs

print(mcf10_LCCs_intra())

def mcf10_components():
    graphs = mcf10_intra_raw_graphs()
    components_instance = gg.ConnectedComponents(graphs)
    components = components_instance.find_components()
    components_instance.print_components()
    return components

print(mcf10_components())

def mcf10_intra_norm_graphs():
    root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/intra/mcf10/norm")
    graph_creator = gg.CreateGraphsFromDirectory(root_dir)
    graph_creator.from_edgelists()
    mcf10_graphs = graph_creator.graph_dict
    return mcf10_graphs

def mcf10_inter_graphs():
    root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/inter/mcf10/nosplit")
    graph_creator = gg.CreateGraphsFromDirectory(root_dir)
    graph_creator.from_edgelists()
    mcf10_graphs = graph_creator.graph_dict
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
print(gsm2824367_intra_graphs())

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
    filtered_graph = graph_filter.filter_graphs(chromosomes="chr18", resolutions="1000000", interaction_type="intra")
    graph_filter.print_filtered_edges()
    return filtered_graph



def imr90_chr1_inter():
    all_imr90 = imr90_inter_graphs()
    graph_filter = gg.FilterGraphs(all_imr90)
    filtered_graph = graph_filter.filter_graphs(chromosomes="chr1", resolutions="1000000")  # , interaction_type="intra")
    graph_filter.print_filtered_edges()
    return filtered_graph

def mcf10_chr1_inter():
    all_mcf10 = mcf10_inter_graphs()
    graph_filter = gg.FilterGraphs(all_mcf10)
    filtered_graph = graph_filter.filter_graphs(chromosomes="chr10", resolutions="1000000", interaction_type="inter")
    graph_filter.print_filtered_edges()
    return filtered_graph


# YESSS IT WORKS !!!!!

# Take graphs, filter on them and combine them into one
def mcf10_combined():
    intra_graphs = mcf10_intra_raw_graphs()
    inter_graphs = mcf10_inter_graphs()

    # Filter the intra_graphs
    filter_intra = gg.FilterGraphs(intra_graphs)
    filtered_intra_graphs = filter_intra.filter_graphs(resolutions="1000000", interaction_type="intra")
    # filter_intra.print_filtered_edges()

    # Filter the inter_graphs
    filter_inter = gg.FilterGraphs(inter_graphs)
    filtered_inter_graphs = filter_inter.filter_graphs(resolutions="1000000", interaction_type="inter")
    # filter_inter.print_filtered_edges()

    # Combine the filtered graphs
    graph_combiner = gg.GraphCombiner([filtered_intra_graphs, filtered_inter_graphs])
    combined_graphs = graph_combiner.combine_matching_graphs()
    # graph_combiner.print_edges(combined_graphs)
    return combined_graphs




if __name__ == "__main__":

    for graph in mcf10_combined():
        print(graph.edges())

    # Export as edge-list
    gg.ExportGraphToEdgelist(mcf10_combined(), "/Users/GBS/Master/HiC-Data/edgelists/others/")






