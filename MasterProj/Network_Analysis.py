"""
Network analysis
"""

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