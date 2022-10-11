""" Plotting node stats """

# Dictionaries created in the node stats function (remove once plotting is complete)
#
# print(node_ratio_dict)
# print(node_ratio_con_dict)
# print(node_ratio_iso_dict)
#
# Make inputs to function "1, 2, 3-7" or "all" etc instead of "chr1, chr2"
# maybe use a decoratr function
# def arg_changer(*args):
#     def wrapper():
#         for arg in args:
#             if args == int:
#
# Plot node_ratio_iso_dict
#
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
#
# Node isolation ratio stacked bar graph
#
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
#
#
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
#
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
