# Import modules
from Network_Analysis import Graph_Processing as Gp
from Network_Analysis import Graph_Analysis as Ga
from pathlib import Path

# Module to store instances of graphs from all cell lines across all resolutions.


# MCF-10A

def mcf10_intra_raw_graphs():
    root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/intra/mcf10/raw")
    graph_creator = Gp.CreateGraphsFromDirectory(root_dir)
    graph_creator.from_edgelists()
    mcf10_graphs = graph_creator.graph_dict
    return mcf10_graphs

def mcf10_intra_norm_graphs():
    root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/intra/mcf10/norm")
    graph_creator = Gp.CreateGraphsFromDirectory(root_dir)
    graph_creator.from_edgelists()
    mcf10_graphs = graph_creator.graph_dict
    return mcf10_graphs

def mcf10_inter_graphs():
    root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/inter/mcf10/nosplit")
    graph_creator = Gp.CreateGraphsFromDirectory(root_dir)
    graph_creator.from_edgelists()
    mcf10_graphs = graph_creator.graph_dict
    return mcf10_graphs

def mcf10_inter_split_graphs():
    root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/inter/mcf10/split")
    graph_creator = Gp.CreateGraphsFromDirectory(root_dir)
    graph_creator.from_edgelists()
    mcf10_graphs = graph_creator.graph_dict
    return mcf10_graphs

# MCF-7

def mcf7_intra_raw_graphs():
    root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/intra/mcf7/raw")
    graph_creator = Gp.CreateGraphsFromDirectory(root_dir)
    graph_creator.from_edgelists()
    mcf7_graphs = graph_creator.graph_dict
    return mcf7_graphs

def mcf7_intra_norm_graphs():
    root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/intra/mcf7/norm")
    graph_creator = Gp.CreateGraphsFromDirectory(root_dir)
    graph_creator.from_edgelists()
    mcf7_graphs = graph_creator.graph_dict
    return mcf7_graphs

def mcf7_inter_graphs():
    root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/inter/mcf7/nosplit")
    graph_creator = Gp.CreateGraphsFromDirectory(root_dir)
    graph_creator.from_edgelists()
    mcf7_graphs = graph_creator.graph_dict
    return mcf7_graphs

def mcf7_inter_split_graphs():
    root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/inter/mcf7/split")
    graph_creator = Gp.CreateGraphsFromDirectory(root_dir)
    graph_creator.from_edgelists()
    mcf7_graphs = graph_creator.graph_dict
    return mcf7_graphs




# IMR90

def imr90_intra_graphs():
    root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/intra/imr90")
    graph_creator = Gp.CreateGraphsFromDirectory(root_dir)
    graph_creator.from_edgelists()
    imr90_graphs = graph_creator.graph_dict
    return imr90_graphs

def imr90_inter_graphs():
    root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/inter/imr90/nosplit")
    graph_creator = Gp.CreateGraphsFromDirectory(root_dir)
    graph_creator.from_edgelists()
    imr90_graphs = graph_creator.graph_dict
    return imr90_graphs

def imr90_inter_split_graphs():
    root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/inter/imr90/split")
    graph_creator = Gp.CreateGraphsFromDirectory(root_dir)
    graph_creator.from_edgelists()
    imr90_graphs = graph_creator.graph_dict
    return imr90_graphs


# HUVEC

def huvec_intra_graphs():
    root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/intra/huvec")
    graph_creator = Gp.CreateGraphsFromDirectory(root_dir)
    graph_creator.from_edgelists()
    huvec_graphs = graph_creator.graph_dict
    return huvec_graphs

def huvec_inter_graphs():
    root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/inter/huvec/nosplit")
    graph_creator = Gp.CreateGraphsFromDirectory(root_dir)
    graph_creator.from_edgelists()
    huvec_graphs = graph_creator.graph_dict
    return huvec_graphs

def huvec_inter_split_graphs():
    root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/inter/huvec/split")
    graph_creator = Gp.CreateGraphsFromDirectory(root_dir)
    graph_creator.from_edgelists()
    huvec_graphs = graph_creator.graph_dict
    return huvec_graphs

# GSM2824367

def gsm2824367_intra_graphs():
    root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/intra/gsm2824367")
    graph_creator = Gp.CreateGraphsFromDirectory(root_dir)
    graph_creator.from_edgelists()
    gsm2824367_graphs = graph_creator.graph_dict
    return gsm2824367_graphs

def gsm2824367_inter_graphs():
    root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/inter/gsm2824367/nosplit")
    graph_creator = Gp.CreateGraphsFromDirectory(root_dir)
    graph_creator.from_edgelists()
    gsm2824367_graphs = graph_creator.graph_dict
    return gsm2824367_graphs

def gsm2824367_inter_split_graphs():
    root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/inter/gsm2824367/split")
    graph_creator = Gp.CreateGraphsFromDirectory(root_dir)
    graph_creator.from_edgelists()
    gsm2824367_graphs = graph_creator.graph_dict
    return gsm2824367_graphs


# Filtered graphs here ?


def imr90_chr18():
    all_imr90 = imr90_intra_graphs()
    graph_filter = Gp.FilterGraphs(all_imr90)
    filtered_graph = graph_filter.filter_graphs(chromosomes="chr15", resolutions="1000000")
    graph_filter.print_filtered_edges()
    return filtered_graph

# imr90_chr18()

def imr90_chr1_inter():
    all_imr90 = imr90_inter_graphs()
    graph_filter = Gp.FilterGraphs(all_imr90)
    filtered_graph = graph_filter.filter_graphs(chromosomes="chr1", resolutions="1000000", interaction_type="inter")
    graph_filter.print_filtered_edges()
    return filtered_graph

imr90_chr1_inter()