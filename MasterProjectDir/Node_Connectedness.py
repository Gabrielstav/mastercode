from Processing import Node
from Processing import Celline
from Processing import Cellines

# Lists for isolated and connected nodes across all chromosomes called in the get_nodes_from_chromosomes function
chromlist_isolated_nodes = []
chromlist_connected_nodes = []
# Isolated nodes are the first value for every key, connected nodes the econd value for every key
node_ratio_dict = {}
node_ratio_iso_dict = {}
node_ratio_con_dict = {}

# Iterate over celline as list


def get_nodes_from_chromosome(*args):

    """
    creates lists containing empty and connected nodes from selected chromosomes, e.g: ("chr1", "chrX")
    and calculates some statistics on empty vs connected nodes for the selected chromosome(s)
    """

    # empty and connected node lists
    for node in nodelist:
        if nodelist.edges == "." and nodelist.chr in args:
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


get_nodes_from_chromosome("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
                          "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")


""" 
Node read length statistics 
"""

# We need to create positional information from the reads in the nodes, so we can sort them e.g by length and
# plot them from start to finish or something? Later we want to integrate RNAseq data so we can overlay this on the nodes?

# Use pybedtools or do it organically in python, think we should avoid hte command line if we can.

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