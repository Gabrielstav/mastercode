# Import modules


# TODO: Find lcc to total size ratio for each graph
#   this class should take a lcc graph and return

# TODO: Calculate betweenness centrality for LCC (on cell line --> resolution level).

# TODO: Calculate closeness centrality for LCC (on cell line --> resolution level).

# TODO: Make a class that does differential community detection for two graphs: Two cell lines on same resolution.

# TODO: Plot the jaccard index for two cell lines on same resolution (for all resolutions).

# TODO: Make Alpaca differential community detection that takes two graphs: Two cell lines on same resolution and calculates the difference in community detection metrics.

# TODO: Plot the differential community detection stuff.

# Move betweenness and closeness to graph analysis? Since this is not necessarily comparative. "But it is comparative in the sense that it is comparing the same graph to itself" -copilot 24/4 23

class BetweennessCentrality:

    def __init__(self, graph_dict):
        self.graph_dict = graph_dict

    def calculate_betweenness_centrality(self):
        betweenness_centrality_dict = {}
        for graph_name, graph in self.graph_dict.items():
            betweenness_centrality_dict[graph_name] = graph.betweenness()
        return betweenness_centrality_dict

    def print_betweenness_centrality(self):
        for graph_name, graph in self.calculate_betweenness_centrality().items():
            print(f"Betweenness centrality for: {graph_name} \n {graph}")

class ClosenessCentrality:

    def __init__(self, graph_dict):
        self.graph_dict = graph_dict

    def calculate_closeness_centrality(self):
        closeness_centrality_dict = {}
        for graph_name, graph in self.graph_dict.items():
            closeness_centrality_dict[graph_name] = graph.closeness()
        return closeness_centrality_dict

    def print_closeness_centrality(self):
        for graph_name, graph in self.calculate_closeness_centrality().items():
            print(f"Closeness centrality for: {graph_name} \n {graph}")

class Degree:
    pass






