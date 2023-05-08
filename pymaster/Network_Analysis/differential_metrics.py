# Import modules
import os as os
import igraph as ig
import rpy2.robjects as robjects
import rpy2.robjects.packages as rpackages
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

# TODO:
#    1. Make alpaca class
#    2. Find other algorithms for differential community detection if time ?

# TODO: ALPACA (ALtered partitons Across Community Architectures)

# 1. Takes graphs with community structure assigned (from fast greedy, louvian, infomap)
# 2. Then it aligns the communities between the graphs, assigning the same node to the same community across graphs
# 3. Then the algo does differential community detection, by finding subgraphs in each input graph that are sig diff
#    from other community strucuees, by analyzing the edge distribution withtin/outside of communities?

# Alpaca needs the community structure as inputs, so use as_clustering to get the struct from the graph obj.
# Then I need to pre-process this, adding 1 as weight for each edge, and 0 for each non-edge (but we have only connected edges).
# I could use bedtools to fill in the network to make a segmented genome, containing mostly 0 edges for all non-connected nodes?

# Make class that takes two graphs as input, calls the fast greedy, louvian or infomap algos
# Pre-processes the comm struct into the correct format, and makes adjacency matrices for each graph
# Then integrate the alpaca package into Python or reverse engineer it to make it work in python
# Then call the alpaca algo on the graphs

# So step-by-step:
# 1. Find community structure for each graph
# 2. Pre-process the community structure into the correct format, convertign to R objects (1 weighted for connection)
# 3. Run alpaca on the R objects
# 4. Process the alpaca reults


# Basic implemenataion (test this first, before making class implementation):

# Take two graphs (actual implementation will take graph dict obj as inputs)
baseline_graph = ig.Graph.GRG(50, 0.2)
perturbed_graph = ig.Graph.GRG(50, 0.3)

# Do comm detecton
baseline_communities = baseline_graph.community_fastgreedy().as_clustering() # Baseline graph has changed in docs, find new way to get it
perturbed_communities = perturbed_graph.community_fastgreedy().as_clustering()

# Convert to R objects
with localconverter(robjects.default_converter + pandas2ri.converter):
    r_baseline_graph = robjects.conversion.py2rpy(baseline_graph)
    r_perturbed_graph = robjects.conversion.py2rpy(perturbed_graph)

# Convert comm struct to R objects
r_baseline_communities = robjects.IntVector(baseline_communities.membership)
r_perturbed_communities = robjects.IntVector(perturbed_communities.membership)

# Run Alpaca
alpaca = rpackages.importr("alpaca")
alpaca_results = alpaca.alpaca(r_baseline_graph, r_perturbed_graph, r_baseline_communities, r_perturbed_communities)

# Process results
differential_modularity_matrix = alpaca_results.rx2("D")
differential_communities = alpaca_results.rx2("M")

# Convert back to python objects
with localconverter(robjects.default_converter + pandas2ri.converter):
    differential_modularity_matrix = robjects.conversion.rpy2py(differential_modularity_matrix)
    differential_communities = robjects.conversion.rpy2py(differential_communities)

# Export results to useful format for plotting and analysis, convert back to igraph objects
differential_modularity_matrix = ig.Graph.Adjacency(differential_modularity_matrix.tolist())
differential_communities = ig.Graph.Adjacency(differential_communities.tolist())

# Make grap dict from results ?






