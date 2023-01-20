"""
Network analysis
"""

import igraph as ig
import numpy as np
from igraph import plot
from matplotlib import pyplot as plt
from collections import Counter
import altair as alt

# Check if we can create graph objects directly without having to create an edgelist format file each time and
# import it from a different directory.

# K562 chromosome 18 (just testing iGraph stuff)
K562_chr18 = ig.Graph.Load("/Users/GBS/Master/HiC-Data/Processed_Data_Edgelist/HiC_from_Jonas/Chr18/K562_processed_chr18.txt",
                           format="ncol")

HUVEC_chr18 = ig.Graph.Load("/Users/GBS/Master/HiC-Data/Processed_Data_Edgelist/HiC_from_Jonas/Chr18/HUVEC_processed_chr18.txt", format="ncol")

IMR90 = ig.Graph.Load("/Users/GBS/Master/HiC-Data/Processed_Data_Edgelist/HiC_from_Jonas/FullGenome/IMR90_processed.txt", format="ncol")

# K562_chr18.save("K562_chr18", format = "ncol") Saving doesn't work?
# Finding degree of graph
degree_K562_chr18 = K562_chr18.degree()
# print(degree_K562_chr18)

# Finding betweenness of edges (same as using .es but our graph is not directed)
edge_betweenness_K562_chr18 = K562_chr18.edge_betweenness()
print(edge_betweenness_K562_chr18)

# Finding betweenness of vertices (.vs: nodes)
nodes_betweenness_K562_chr18 = K562_chr18.vs.betweenness()
print(nodes_betweenness_K562_chr18)

# Finding adjacency matrix for graph
adjacency_K562_chr18 = K562_chr18.get_adjacency()
# print(adjacency_K562_chr18)

# Is our graph directed:
print("Our graph is directed:", K562_chr18.is_directed())

# Testing if the Cairo package works (it works)
# g = ig.Graph.Famous("petersen")
# plot(g)
# And it works on my data
# plot(K562_chr18, layout = "fr", directed="false") # 2D layouts can be: circle, drl, fr, kk, lgl, random, rt

# plot(IMR90, layout = "")
# plot(HUVEC_chr18, layout = "fr")


# Frequency distribution test

# Må ha frequency på Y-aksen og Degree på X-aksen
# Må ogsp ha dotplot og ikke histogram, det ser ut som de har log-transformert grafen i paperet.
# Sjekk paper for å se hvordan det skal se ut: https://pubmed.ncbi.nlm.nih.gov/27618581/

# Plot this using Altair? Looked nice on OMgenomics.

# Something is wrong with iGraph, troubleshoot later: https://github.com/scverse/scanpy/issues/961

ig.add_vertices(5)
ig.add_edges([(0,1), (0,2), (0,3), (1,2), (2,3), (3,4)])

# Compute the degree distribution using the Counter class
degree_values = g.degree()
degree_counts = Counter(degree_values)

# Create a dataframe with the degree and frequency data
df = pd.DataFrame({"degree": list(degree_counts.keys()), "frequency": list(degree_counts.values())})

# Plot the degree distribution using the Chart and Scatter classes
chart = alt.Chart(df).mark_circle(size=60).encode(
    x=alt.X("degree:Q", scale=alt.Scale(type="log")),
    y=alt.Y("frequency:Q", scale=alt.Scale(type="log")),
)
chart.show()

# bins = 20
# # log_data = np.log(K562_chr18)
# plt.scatter(K562_chr18.degree(), bins)
# plt.semilogx()
# plt.xscale("log")
# plt.ylim(0.1, 100)
# plt.xlim(9, max(xscale))
# plt.ylabel("Degree")
# plt.xlabel("Frequency")
# plt.show()

# Something like this for degree distribution?

# plt.scatter(degree_counts.keys(), degree_counts.freq(), linewidth=0)
# plt.semilogx(degree_counts.keys(), degree_counts.values())
# plt.semilogy(degree_counts.keys, degree_counts.values())
# plt.xlabel("Degree")
# plt.ylabel("Frequency")
# plt.show()















# After we standardize nodes, we can run community detection (e.g fast greedy first), then look at
# clustering, degree distritution, betweeness centrality and other network metrics.

# communities = g.communnity_fastgreedy()
# print(communities)
# communitites_list = communitites.as_clustering() # communitites as list of vertices
# The index represents the ID of each edge:
# membership = communitites.membership
# print(membership)

# Plot the graph, coloring the vertices by their community
# plot_options = {
#     "vertex_color": communities.membership,
#     "vertex_label": g.vs["name"],
#     "vertex_size": 30,
#     "edge_color": "lightgray",
#     "margin": 20
# }
# ig.plot(g, **plot_options)




