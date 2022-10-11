""" Importing packages"""

from dataclasses import dataclass, astuple, fields
import matplotlib.pyplot as plt
import numpy as np
import pybedtools as pb
import os
from os import path

import pandas as pd
from igraph import *
import re

from cairo import *

# import ast
# import csv
# import pandas as pd
# import io


""" Creating classes """


@dataclass
class Node:
    id: str
    edges: list[str]
    chr: str

    def __iter__(self):
        return iter(astuple(self))

    def print_class(self):
        print("Node '" + self.id + "', chromosome '" + self.chr + "'")




@dataclass
class NodeList:
    nodes: list

    def __iter__(self):
        return iter(astuple(self))

    def as_list(self):
        return self.nodes

    def group_cellines(self):
        pass

    def group_chromosomes(self):
        pass

    def group_connected_nodes(self):
        pass

    def group_isolated_nodes(self):
        pass

    def isolated_degree(self):
        pass

    def node_overlap(self):
        pass

    def standardize_nodes(self):
        pass


@dataclass
class Celline:
    strain: str
    nodes: NodeList

    def __iter__(self):
        return self

    def print_self(self):
        print(f"Strain: {self.strain} with nodelist: {self.NodeList}")

Celline.print_self(NodeList)


""" Pre-processing data from gtrack to edgelist """


# Selecting files to pre-process and instantiating Node class
def process_file_to_node(*args):
    nodes = []
    for arg in args:
        with open(arg, encoding="latin-1") as file:
            file_content: list[str] = file.readlines()
    for index, line in enumerate(file_content):
        if line.startswith("chr"):
            columns = line.strip("\n").split("\t")
            if len(columns) != 7:
                print(f"Bad line format: {columns} with index: {index}")
            else:
                if columns[6] == ".":
                    nodes.append(Node(columns[3], columns[6], columns[0]))
                else:
                    nodes.append(Node(columns[3], columns[6].split(";"), columns[0]))
    return NodeList(nodes)

nodes = process_file_to_node("/Users/GBS/Master/HiC-Data/Hi-C_data_fra_Jonas/4linescopy/IMR90_50kb.domain.RAW.no_cen.NCHG_fdr.o_by_e5_to_plot.gtrack")

# Secelting directory containing files to process and instantiating Celline class
def process_directory_to_celline(directory):
    paths = []
    for file in os.listdir(directory):
        paths.append(os.path.join(directory, file))
    return process_files_to_celline(paths)

def process_files_to_celline(files):
    cellines = []
    for arg in files:
        with open(arg):
            strain = arg.split("/")[7].split("_")[0]
            if not strain.startswith("."):
                cellines.append(Celline(strain, process_file_to_node(arg)))
    return cellines

cellines = process_directory_to_celline("/Users/GBS/Master/HiC-Data/Hi-C_data_fra_Jonas/4cell_lines_Hi-C")




""" Statistics on empty nodes and connectedness ratio (empty vs connected))"""

# Lists for isolated and connected nodes across all chromosomes called int he get_nodes_from_chromosomes function
chromlist_isolated_nodes = []
chromlist_connected_nodes = []
# Isolated nodes are the first value for every key, connected nodes the econd value for every key
node_ratio_dict = {}
node_ratio_iso_dict = {}
node_ratio_con_dict = {}



def get_nodes_from_chromosome(*args):

    """
    creates lists containing empty and connected nodes from selected chromosomes, e.g: ("chr1", "chrX")
    and calculates some statistics on empty vs connected nodes for the selected chromosome(s)
    """

    # empty and connected node lists
    for node in nodes.as_list():
        if node.edges == "." and node.chr in args:
            chromlist_isolated_nodes.append(node)
        elif node.edges != "." and node.chr in args:
            chromlist_connected_nodes.append(node)

    # stats for all chromosomes combined (for comparing cell lines)
    total_nodes = len(chromlist_isolated_nodes) + len(chromlist_connected_nodes)

    def sum_total_nodes(*args):
        for arg in args:
            # sum total nodes
            print(f"Sum of nodes in {arg} = {total_nodes}")
            # sum empty nodes
            print(f"Isolated nodes in {arg} = {len(chromlist_isolated_nodes)}")
            # sum connected nodes
            print(f"Connected nodes in {arg} = {len(chromlist_connected_nodes)}")
            # sum connected node ratio
            total_connectedness_ratio = len(chromlist_connected_nodes) / total_nodes
            print(f"Connectedness ratio in {arg} = {round(total_connectedness_ratio, 5)}")
            # sum isolated node ratio
            total_isolated_ratio = len(chromlist_isolated_nodes) / total_nodes
            print(f"Isolated ratio in {arg} = {round(total_isolated_ratio, 5)}")

    sum_total_nodes()

    # nodes per chromosome
    def get_nodes_in_chrom(*args, nodelist):
        nodes_in_total = []
        for arg in args:
            get_chrom = list(filter(lambda x: arg == str(x.chr), nodes))
            print(f"Chromosome {arg} has {len(get_chrom)} nodes in total")
            nodes_in_total.append(len(get_chrom))
        # Test (unit test this some day?)
        if sum(nodes_in_total) != total_nodes:
            print(
                f"ERROR: len(nodes_in_total) = {sum(nodes_in_total)} and len(total_nodes) = {total_nodes} are different."
                f" Maybe you duplicated chromosomes when calling the function?")

    # get_nodes_in_chrom(*args, nodelist=nodes)

    # Empty to connected ratio per chromosome
    def empty_nodes_in_chrom(*args, nodes):
        for arg in args:
            get_chroms = list(filter(lambda x: arg == str(x.chr), nodes.as_list()))

            # number of isolated nodes
            get_iso = list(filter(lambda y: y.edges == ".", get_chroms))
            print(f"Chromosome {arg} has {len(get_iso)} isolated nodes")

            # number of connected nodes
            get_con = list(filter(lambda z: z.edges != ".", get_chroms))
            print(f"Chromosome {arg} has {len(get_con)} connected nodes")

            # isolated node ratio
            iso_ratio = len(get_iso) / len(get_con + get_iso)
            print(f"Chromosome {arg} has isolation ratio = {round(iso_ratio, 5)}")

            # connective node ratio
            con_ratio = len(get_con) / len(get_con + get_iso)
            print(f"Chromosome {arg} has connective ratio = {round(con_ratio, 5)}")

            # outputting connective and isolated node ratios to global dict for plotting
            def add_element_node_ratio_dict(dictionary, key, iso_value, con_value):
                if key not in node_ratio_dict:
                    node_ratio_dict[key] = []
                dictionary[key].append(iso_value)
                dictionary[key].append(con_value)

            add_element_node_ratio_dict(node_ratio_dict, arg, iso_ratio, con_ratio)

            # creating isolated node ratio dict for plotting
            def add_element_node_ratio_dict_iso(dictionary, key, value):
                if key not in node_ratio_iso_dict:
                    node_ratio_iso_dict[key] = []
                dictionary[key].append(value)

            add_element_node_ratio_dict_iso(node_ratio_iso_dict, arg, iso_ratio)

            # creating connected ratio dictionary for plotting
            def add_element_node_ratio_dict_con(dictionary, key, value):
                if key not in node_ratio_con_dict:
                    node_ratio_con_dict[key] = []
                dictionary[key].append(value)

            add_element_node_ratio_dict_con(node_ratio_con_dict, arg, con_ratio)

            # tests (because IDK how to unit testing yet)
            if con_ratio + iso_ratio != 1:
                print(f"ERROR: sum of con+iso ratios are not equal to 1. They are = {con_ratio + iso_ratio}")
            if len(get_chroms) != len(get_iso) + len(get_con):
                print(f"ERROR: len(get_chroms) = {len(get_chroms)} is not equal to len of isolated"
                      f" + connected nodes {len(get_iso) + len(get_con)}")

    empty_nodes_in_chrom(*args, nodes=nodes)


# get_nodes_from_chromosome("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
#                           "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")

""" Plotting node stats """

# Dictionaries created in the node stats function (remove once plotting is complete)

# print(node_ratio_dict)
# print(node_ratio_con_dict)
# print(node_ratio_iso_dict)

# Make inputs to function "1, 2, 3-7" or "all" etc instead of "chr1, chr2"
# maybe use a decoratr function
# def arg_changer(*args):
#     def wrapper():
#         for arg in args:
#             if args == int:

# Plot node_ratio_iso_dict

# print(node_ratio_iso_dict)
# def select_chromosomes_for_plotting_iso(*args):
#         # Selecting chromosomes for plotting
#
#     node_ratio_keys = []
#     node_ratio_values = []
#     for arg in args:
#         for key, val in node_ratio_dict.items():
#             if arg == key:
#                 node_ratio_keys.append(key)
#                 node_ratio_values.append(val)
#     node_ratio_values = np.array(node_ratio_values)
#
#     # Plotting node ratio
#     fig, ax = plt.subplots()
#     ax.scatter(range(len(node_ratio_keys)), node_ratio_values[:, 1])
#     # plt.xticks(range(len(node_ratio_keys)), node_ratio_keys)
#     plt.show()
#
# select_chromosomes_for_plotting_iso("chr1", "chrX", "chr18", "chr5")

# Node isolation ratio stacked bar graph

# Vil plotte hele ratioen (altså node_ratio_dict) ikke bare node_ratio_iso_dict
# def select_chromosomes_for_plotting(*args): #og celline?
#
#     # Selecting chromosomes for plotting
#     node_ratio_keys = []
#     node_ratio_values = []
#     for arg in args:
#         print(arg)
#         for key, val in node_ratio_dict.items():
#             if arg == key:
#                 node_ratio_keys.append(key)
#                 node_ratio_values.append(val)
#     node_ratio_values = np.array(node_ratio_values)
#
#     # Plotting node ratio
#     fig, ax = plt.subplots()
#     plt.bar(range(len(node_ratio_keys)), node_ratio_values[:, 0])
#     plt.bar(range(len(node_ratio_keys)), node_ratio_values[:, 1], bottom=node_ratio_values[:, 0])
#     plt.xticks(range(len(node_ratio_keys)), node_ratio_keys)
#     plt.show()
#
#
#
# select_chromosomes_for_plotting("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
#                                "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")


# # Trying to add subtitle to each bar, but nothing showed up?
#     y_offset = -15
#     for bar in ax.patches:
#         ax.text(
#             bar.get_x() + bar.get_width() / 2,
#             bar.get_height() + bar.get_y() + y_offset, round(bar.get_height()),
#             ha="center",
#             color="k",
#             weight="bold",
#             size=8
#         )

# # Forsøk på å sortere etter verdi, men: '<' not supported between instances of 'numpy.ndarray' and 'list'
#     x = df.index
#     indexes = np.argsort(df.values).T
#     heights = np.sort(df.values).T
#     order = -1
#     bottoms = heights[::order].cumsum(axis=0)
#     bttoms = np.insert(bottoms, 0, np.zeros(len(bottoms[0])), axis=0)
#     mpp_colors = dict(zip(df.columns, plt.rcParams["axes.prop_cycle"].by_key()["color"]))
#     for btms, (idxs, vals) in enumerate(list(zip(indexes, heights))[::order]):
#         mps = np.take(np.array(df.columns), idxs)
#         ax.bar(x, height=vals, bottom=bottoms[btms], color=[mpp_colors[m] for m in mps])
#     ax.set_ylim(bottom=0, top=2)
#     plt.legend((np.take(np.array(df.columns), np.argsort(df.values)[0]))[::order], loc="upper right")
#
#     ax.set_title("Node connectedness in IMR90")
#     plt.xlabel("Chromosome")
#     plt.ylabel("Connectedness ratio")
#     plt.show()


""" Node read length statistics """

# We need to create positional information from the reads in the nodes, so we can sort them e.g by length and
# plot them from start to finish or something? Later we want to integrate RNAseq data so we can overlay this on the nodes?

# Use pybadtools or do it organically in python, think we should avoid hte command line if we can.

# Can we just write a function in python that takes each connected node as input from the data class. The nodes are already ordered from start to end of chromosome/genome.
# Then we calculate the length of the first node for the 4 different cell lines.
# Then we quantify their positions, by using the start/end point.
# Maybe the largest node is the anchor node, and the smaller nodes are matched against it?
# After we have length and positions, we can calculate the percentage of overlap.
# If nodes overlap less than say 50%, the node is discarded (not deleted from the data class, but just not included in the function output).

# We want to examine node overlap both on the "node level" (whatever that means) and on the bp level.
# Find overlap using pybadtools intersect? or just write the function to calculate the percentage overlap/ bp overlap.
# Do this for every node, iterating through each instance of the data class for the 4 cell lines.
# Before we do this we need to write functionality for reading in multiple cell lines and store their data as instances of the Node class.
# Then we can iterate over the instances, for each cell line, all the way through the data class.
# Inside this iterator, we need to do the percentage/bp/node overlap calculations
# OR
# We can output the mapping of the nodes (start and end to find length) for each isntance to a dict or something,
# and then use this for further coverage calculations.


""" Writing pre-processed output to new files """

# Use to OS module to check what directory we are in, read the files put in a specific directory (to be pre-processed),
# checks their format (only .gtrack compatibility as of now) and then run the pre-processing parts of the code,
# and outputs pre-processed files in the edge list format to a specific directory with specific file names,
# maybe just the original file names with "processed" or "edge-list" behind the name.


# Pseudocode-function to automatically process files and write to a new dir if we provide the dir:
# count = 0
# directory = "/Users/GBS/Master/HiC-Data/Processed_Data/HiC_from_Jonas/FullGenome"
# for filename in os.listdir(directory):
#     if filename.endswith(".txt"):
#         count += 1
#         with open(directory + filename, "r") as read_file:
#            return_of_your_function = do_something_with_data() # Putt pre-processinga over inn en funksjon som calles her
#         with open(directory + count + filename, "w") as write_file:
#             write_file.write(return_of_your_function)  !!! igjen bare caller funksjonen her !!!

# Selecting full genome and writing edge list to new file
# with open("/Users/GBS/Master/HiC-Data/Processed_Data/HiC_from_Jonas/FullGenome/TESTING", "w") as fp:
#     for node in nodes:
#         for edge in node.edges:
#             fp.write((node.id + " " + edge + "\n"))


# Selecting specific chromsome and writing edge list to new file
# with open("/Users/GBS/Master/HiC-Data/Processed_Data/Processed_J_data/HiC_from_Jonas/Chr18/IMR90_processed_chr18.txt", "w") as fp:
#     for l in nodes:
#         for edge in l.edges:
#              if l.chr == "chr18":  # Make a function that we can call to select id + edges from chromosome(s)?
#                  fp.write((l.id + " " + edge + "\n"))

# This outputs a file with the correct format but with header: id, edges (which should be removed).
# Use with open(file) that points to file opened earlier, should be able to find the path to the file opened when
# pre-processing the data, so that we only need to point to the full path once and the rest of the code knows what
# file(s) we are working with.


""" Network analysis """

# Check if we can create graph objects directly without having to create an edgelist format file each time and
# import it from a different directory.

# # K562 chromosome 18 (just testing iGraph stuff)
# K562_chr18 = Graph.Load("/Users/GBS/Master/HiC-Data/Processed_Data/HiC_from_Jonas/Chr18/K562_processed_chr18.txt",
#                         format="ncol")
# # K562_chr18.save("K562_chr18", format = "ncol") Saving doesn't work?
# # Finding degree of graph
# degree_K562_chr18 = K562_chr18.degree()
# #print(degree_K562_chr18)
#
# # Finding betweenness of edges (same as using .es but our graph is not directed)
# edge_betweenness_K562_chr18 = K562_chr18.edge_betweenness()
# print(edge_betweenness_K562_chr18)
#
# # Finding betweenness of vertices (.vs: nodes)
# nodes_betweenness_K562_chr18 = K562_chr18.vs.betweenness()
# print(nodes_betweenness_K562_chr18)

# # Finding adjacency matrix for graph
# adjacency_K562_chr18 = K562_chr18.get_adjacency()
# print(adjacency_K562_chr18)

# # Is our graph directed:
# print("Our graph is directed:", K562_chr18.is_directed())

# Testing if the Cairo package works (it works)
# g = Graph.Famous("petersen")
# plot(g)
# And it works on my data
# plot(K562_chr18, layout = "fr") # 2D layouts can be: circle, drl, fr, kk, lgl, random, rt


# Frequency distribution test

# Må ha frequency på Y-aksen og Degree på X-aksen
# Må ogsp ha dotplot og ikke histogram, det ser ut som de har log-transformert grafen i paperet.
# Sjekk paper for å se hvordan det skal se ut: https://pubmed.ncbi.nlm.nih.gov/27618581/

# bins = 20
# plt.hist(K562_chr18.degree(), bins)
# plt.show()

# After we standardize nodes, we can run community detection (e.g fast greedy first), then look at
# clustering, degree distritution, betweeness centrality and other network metrics.
