# Import modules
import scipy.stats as stats
import numpy as np
import igraph as ig

# Take in any two graphs and calculate the differences between them
# I think just make use of spearman and pearson and jaccard for centrality metrics
# And for topology metrics just use Jaccard


# TODO:
#    1. Make Jaccard class (does not make that much sense)
#    2. Make Pearson class (compares linear releationship between continuous variables)
#    3. Make Spearman class (compares monotonic releationship between ranked data)
#    4. Make KS test class (compares distributions between two samples, finds if they are drawn from the same distribution)

# This will take too long to implement and is computationally expensive, but could be done after project:
# TODO:
#    5. Make cosine similarity class (angle between n-dimensional vectors, closer to 1 means more similar, can plot two graphs against each other)
#    6. Maake Euclidean distance class (Straight line distance between two n-dimensional vectors, close to 0 means more similar)
#    7. Make graph edit distance class (measure of dissimilariry, transforms one graph into another, comp expensive)
#    8. Make maximum common subgraph class (finds largest common subgraph between two graphs, comp expensive)
#    9. Graph kernels but too time consuming and not that useful and comp expensive


# Example:

graph1 = ig.Graph.GRG(50, 0.2)
graph2 = ig.Graph.GRG(50, 0.3)

degree_seq1 = np.array(graph1.degree())
degree_seq2 = np.array(graph2.degree())

pearson_corr, _ = stats.pearsonr(degree_seq1, degree_seq2)
print("Pearson correlation:", pearson_corr)


class JaccardIndex:

    def __init__(self, graph_dict):
        self.graph_dict = graph_dict

    def calculate_jaccard_index(self):
        jaccard_index_dict = {}
        for graph_name, graph in self.graph_dict.items():
            jaccard_index_dict[graph_name] = graph.similarity_jaccard()
        return jaccard_index_dict

    def print_jaccard_index(self):
        for graph_name, graph in self.calculate_jaccard_index().items():
            print(f"Jaccard index for: {graph_name} \n {graph}")


class PearsonCorrelation:

    def __init__(self, graph_dict):
        self.graph_dict = graph_dict

    def calculate_pearson_correlation(self):
        pearson_correlation_dict = {}
        for graph_name, graph in self.graph_dict.items():
            pearson_correlation_dict[graph_name] = graph.similarity_jaccard()
        return pearson_correlation_dict

    def print_pearson_correlation(self):
        for graph_name, graph in self.calculate_pearson_correlation().items():
            print(f"Pearson correlation for: {graph_name} \n {graph}")


class SpearmanCorrelation:

    def __init__(self, graph_dict):
        self.graph_dict = graph_dict