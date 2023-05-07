# Import modules
from Network_Analysis import graph_processing as gp
# from Network_Analysis import Graph_Analysis as Ga
from pathlib import Path

# Module to store instances of graphs from all cell lines across all resolutions.


# MCF-10A

# def mcf10_combined():
#     graph_combiner = gp.GraphCombiner([mcf10_intra_raw_graphs(), mcf10_inter_graphs()])
#     combined_graphs = graph_combiner.combine_intra_inter_graphs()
#     return combined_graphs

def mcf10_intra_raw_graphs():
    root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/intra/mcf10/raw")
    graph_creator = gp.CreateGraphsFromDirectory(root_dir)
    graph_creator.from_edgelists()
    mcf10_graphs = graph_creator.graph_dict
    return mcf10_graphs
# print(mcf10_intra_raw_graphs())

def mcf10_intra_norm_graphs():
    root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/intra/mcf10/norm")
    graph_creator = gp.CreateGraphsFromDirectory(root_dir)
    graph_creator.from_edgelists()
    mcf10_graphs = graph_creator.graph_dict
    return mcf10_graphs

def mcf10_inter_graphs():
    root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/inter/mcf10/nosplit")
    graph_creator = gp.CreateGraphsFromDirectory(root_dir)
    graph_creator.from_edgelists()
    mcf10_graphs = graph_creator.graph_dict
    return mcf10_graphs
# print(mcf10_inter_graphs())

def mcf10_inter_split_graphs():
    root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/inter/mcf10/split")
    graph_creator = gp.CreateGraphsFromDirectory(root_dir)
    graph_creator.from_edgelists()
    mcf10_graphs = graph_creator.graph_dict
    return mcf10_graphs

# MCF-7

# def mc7_combined():
#     graph_combiner = gp.GraphCombiner([mcf7_intra_raw_graphs(), mcf7_inter_graphs()])
#     combined_graphs = graph_combiner.combine_intra_inter_graphs()
#     return combined_graphs

def mcf7_intra_raw_graphs():
    root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/intra/mcf7/raw")
    graph_creator = gp.CreateGraphsFromDirectory(root_dir)
    graph_creator.from_edgelists()
    mcf7_graphs = graph_creator.graph_dict
    return mcf7_graphs

def mcf7_intra_norm_graphs():
    root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/intra/mcf7/norm")
    graph_creator = gp.CreateGraphsFromDirectory(root_dir)
    graph_creator.from_edgelists()
    mcf7_graphs = graph_creator.graph_dict
    return mcf7_graphs

def mcf7_inter_graphs():
    root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/inter/mcf7/nosplit")
    graph_creator = gp.CreateGraphsFromDirectory(root_dir)
    graph_creator.from_edgelists()
    mcf7_graphs = graph_creator.graph_dict
    return mcf7_graphs

def mcf7_inter_split_graphs():
    root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/inter/mcf7/split")
    graph_creator = gp.CreateGraphsFromDirectory(root_dir)
    graph_creator.from_edgelists()
    mcf7_graphs = graph_creator.graph_dict
    return mcf7_graphs


# IMR90

# def imr90_combined():
#     graph_combiner = gp.GraphCombiner([imr90_intra_graphs(), imr90_inter_graphs()])
#     combined_graphs = graph_combiner.combine_intra_inter_graphs()
#     return combined_graphs

def imr90_intra_graphs():
    root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/intra/imr90")
    graph_creator = gp.CreateGraphsFromDirectory(root_dir)
    graph_creator.from_edgelists()
    imr90_graphs = graph_creator.graph_dict
    return imr90_graphs

def imr90_inter_graphs():
    root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/inter/imr90/nosplit")
    graph_creator = gp.CreateGraphsFromDirectory(root_dir)
    graph_creator.from_edgelists()
    imr90_graphs = graph_creator.graph_dict
    return imr90_graphs

def imr90_inter_split_graphs():
    root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/inter/imr90/split")
    graph_creator = gp.CreateGraphsFromDirectory(root_dir)
    graph_creator.from_edgelists()
    imr90_graphs = graph_creator.graph_dict
    return imr90_graphs


# HUVEC

# def huvec_combined():
#     graph_combiner = gp.GraphCombiner([huvec_intra_graphs(), huvec_inter_graphs()])
#     combined_graphs = graph_combiner.combine_intra_inter_graphs()
#     return combined_graphs

def huvec_intra_graphs():
    root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/intra/huvec")
    graph_creator = gp.CreateGraphsFromDirectory(root_dir)
    graph_creator.from_edgelists()
    huvec_graphs = graph_creator.graph_dict
    return huvec_graphs
# print(huvec_intra_graphs())

def huvec_inter_graphs():
    root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/inter/huvec/nosplit")
    graph_creator = gp.CreateGraphsFromDirectory(root_dir)
    graph_creator.from_edgelists()
    huvec_graphs = graph_creator.graph_dict
    return huvec_graphs
# print(huvec_inter_graphs())

def huvec_inter_split_graphs():
    root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/inter/huvec/split")
    graph_creator = gp.CreateGraphsFromDirectory(root_dir)
    graph_creator.from_edgelists()
    huvec_graphs = graph_creator.graph_dict
    return huvec_graphs

# GSM2824367

# def gsm2824367_combined():
#     graph_combiner = gp.GraphCombiner([gsm2824367_intra_graphs(), gsm2824367_inter_graphs()])
#     combined_graphs = graph_combiner.combine_intra_inter_graphs()
#     return combined_graphs

def gsm2824367_intra_graphs():
    root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/intra/gsm2824367")
    graph_creator = gp.CreateGraphsFromDirectory(root_dir)
    graph_creator.from_edgelists()
    gsm2824367_graphs = graph_creator.graph_dict
    return gsm2824367_graphs

def gsm2824367_inter_graphs():
    root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/inter/gsm2824367/nosplit")
    graph_creator = gp.CreateGraphsFromDirectory(root_dir)
    graph_creator.from_edgelists()
    gsm2824367_graphs = graph_creator.graph_dict
    return gsm2824367_graphs

def gsm2824367_inter_split_graphs():
    root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/inter/gsm2824367/split")
    graph_creator = gp.CreateGraphsFromDirectory(root_dir)
    graph_creator.from_edgelists()
    gsm2824367_graphs = graph_creator.graph_dict
    return gsm2824367_graphs


# Filtered graphs here ?


# Intra filtering works
def imr90_chr18():
    all_imr90 = imr90_intra_graphs()
    graph_filter = gp.FilterGraphs(all_imr90)
    filtered_graph = graph_filter.filter_graphs(chromosomes="chr18", resolutions="1000000", interaction_type="intra")
    graph_filter.print_filtered_edges()
    return filtered_graph

# imr90_chr18()


def imr90_chr1_inter():
    all_imr90 = imr90_inter_graphs()
    graph_filter = gp.FilterGraphs(all_imr90)
    filtered_graph = graph_filter.filter_graphs(chromosomes="chr1", resolutions="1000000")  # , interaction_type="intra")
    graph_filter.print_filtered_edges()
    return filtered_graph

# imr90_chr1_inter()

def mcf10_chr1_inter():
    all_mcf10 = mcf10_inter_graphs()
    graph_filter = gp.FilterGraphs(all_mcf10)
    filtered_graph = graph_filter.filter_graphs(chromosomes="chr10", resolutions="1000000", interaction_type="inter")
    graph_filter.print_filtered_edges()
    return filtered_graph
# mcf10_chr1_inter()

# def mcf10_chr1_combined():
#     graph_combiner = gp.GraphCombiner([mcf10_intra_raw_graphs(), mcf10_inter_graphs()])
#     filtered_combined_graph = graph_combiner.filter_combined_graph(chromosomes="chr1", resolutions="1000000")
#     return filtered_combined_graph
#
# print(mcf10_chr1_combined())


# YESSS IT WORKS !!!!!
def mcf10_chr18_combined_single_function():
    intra_graphs = mcf10_intra_raw_graphs()
    inter_graphs = mcf10_inter_graphs()

    # Filter the intra_graphs
    filter_intra = gp.FilterGraphs(intra_graphs)
    filtered_intra_graphs = filter_intra.filter_graphs(resolutions="500000", interaction_type="intra")
    filter_intra.print_filtered_edges()

    # Filter the inter_graphs
    filter_inter = gp.FilterGraphs(inter_graphs)
    filtered_inter_graphs = filter_inter.filter_graphs(interaction_type="inter")
    filter_inter.print_filtered_edges()

    # Combine the filtered graphs
    graph_combiner = gp.GraphCombiner([filtered_intra_graphs, filtered_inter_graphs])
    combined_graphs = graph_combiner.combine_matching_graphs()
    graph_combiner.print_edges(combined_graphs)
    return combined_graphs

print(mcf10_chr18_combined_single_function())




