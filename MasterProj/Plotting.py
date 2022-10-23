
import igraph as ig
import matplotlib as mpl
import pandas as pd
import

"""
Plotting node stats
"""

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


""" 
Node read length statistics 
"""

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


