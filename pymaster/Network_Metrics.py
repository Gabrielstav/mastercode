# Import modules
import igraph as ig
from matplotlib import pyplot as plt
from pathlib import Path
import re as re
import types as types
from typing import Union
import networkx as nx
import plotly.graph_objs as go
import pandas as pd
import os as os

# Set backend of igraph to matplotlib:
ig.config["plotting.backend"] = "matplotlib"

class CreateGraphsFromDirectory:
    """
    Class to create igraph objects from edgelists in a directory.
    """

    def __init__(self, input_directory):
        self.input_directory = Path(input_directory).absolute()
        self.graph_dict = {}

    def get_files(self, pattern):
        files = (file_path for file_path in self.input_directory.rglob(pattern) if file_path.is_file())
        files_list = [str(file_path.resolve()) for file_path in files]
        return files_list

    def from_edgelists(self, cell_lines=None, chromosomes=None, resolutions=None):
        all_files = self.get_files("*")

        for file_path in all_files:
            file_name = os.path.basename(file_path)
            graph_name = re.sub(r"_edgelist\.txt$", "", file_name)

            if cell_lines is not None:
                if not any(cell_line in graph_name for cell_line in cell_lines):
                    continue

            if resolutions is not None:
                graph_name_resolution = re.search(r"_([\d]+)$", graph_name).group(1)
                if graph_name_resolution not in resolutions:
                    continue

            with open(file_path, "r") as f:
                edges = [tuple(filter(None, line.strip().split())) for line in f]

            df = pd.DataFrame(edges, columns=['source', 'target'])

            # Filter the edges based on chromosomes
            if chromosomes is not None:
                df = df[df['source'].str.startswith(tuple(chromosomes)) & df['target'].str.startswith(tuple(chromosomes))]

            graph = ig.Graph.TupleList(df.itertuples(index=False), directed=False)
            graph.vs['name'] = [v['name'] for v in graph.vs]  # Assign name attribute to edges
            self.graph_dict[graph_name] = graph


################################
# Creating graphs from edgelists
################################

# Creating graph objects from raw and normalized edge lists, as well as filtered on resolution and chromosome:
# Move this to another module later.

def mcf7_10_lowres_graphs():
    root_dir = Path("/Users/GBS/Master/HiC-Data/edgelists/lowres_mcf7_mcf10")
    graph_creator = CreateGraphsFromDirectory(root_dir)
    graph_creator.from_edgelists()
    mcf7_10_graphs = graph_creator.graph_dict
    return mcf7_10_graphs


class FilterGraphs:
    """
    Class that filters the graphs based on cell line, chromosome, and resolution.
    """

    def __init__(self, graph_dict):
        self.graph_dict = graph_dict


    def filter_graphs(self, cell_lines=None, chromosomes=None, resolutions=None, graph_dict=None):
        if graph_dict is None:
            graph_dict = self.graph_dict

        if cell_lines:
            if isinstance(cell_lines, str):
                cell_lines = [cell_lines]
            cell_lines = set(cell_lines)
        if chromosomes:
            if isinstance(chromosomes, str):
                chromosomes = [chromosomes]
            chromosomes = set(chromosomes)
        if resolutions:
            if isinstance(resolutions, str):
                resolutions = [resolutions]
            resolutions = set(resolutions)

        filtered_graph_dict = {}

        for graph_name, graph in graph_dict.items():

            if cell_lines is not None and not any(cell_line in graph_name for cell_line in cell_lines):
                continue

            if resolutions is not None:
                graph_name_resolution = re.search(r"_([^_]+)$", graph_name).group(1)
                if graph_name_resolution not in resolutions:
                    continue

            if chromosomes is not None:
                # Create a new empty graph with the same vertices as the original graph
                new_graph = ig.Graph()
                new_graph.add_vertices(graph.vs['name'])
                new_edges = [(graph.vs[e.source]['name'], graph.vs[e.target]['name']) for e in graph.es
                             if (graph.vs[e.source]['name'].split(':')[0] in chromosomes
                                 or graph.vs[e.target]['name'].split(':')[0] in chromosomes)]
                new_graph.add_edges(new_edges)
                # new_graph = new_graph.simplify()  # Remove self-loops and multiple edges if any
                graph = new_graph

            filtered_graph_dict[graph_name] = graph

        return filtered_graph_dict

    @staticmethod
    def print_filtered_edges(filtered_graph_dict):
        for graph_name, graph in filtered_graph_dict.items():
            print(f"Edges in filtered graph {graph_name}:")
            for edge in graph.es:
                source = graph.vs[edge.source]['name']
                target = graph.vs[edge.target]['name']
                print(f"{source} -- {target}")

# MCF10 chr18 1Mbp norm:
def mcf7_chr18_1mb():
    graph_filter = FilterGraphs(mcf7_10_lowres_graphs())
    filtered_graphs = graph_filter.filter_graphs(cell_lines=["mcf7"], chromosomes=["chr18", "chrX"], resolutions=["1000000"])
    graph_filter.print_filtered_edges(filtered_graphs)
    return filtered_graphs
mcf7_chr18_1mb()


def mcf10_chr18_1mb():
    graph_filter = FilterGraphs(mcf7_10_lowres_graphs())
    filtered_graphs = graph_filter.filter_graphs(cell_lines=["mcf10"], chromosomes=["chr18"], resolutions=["1000000"])
    graph_filter.print_filtered_edges(filtered_graphs)
# mcf10_chr18_1mb()

class NetworkMetrics:
    """
    Class to calculate network metrics for a given graph object.
    """

    def __init__(self, graph_dict_or_function, metrics=None):
        if isinstance(graph_dict_or_function, types.FunctionType):
            self.graph_dict = graph_dict_or_function()
        else:
            self.graph_dict = graph_dict_or_function
        self.metrics_dict = metrics if metrics is not None else {}

    available_metrics = {
        "size": lambda g: g.vcount(),
        "edges": lambda g: g.ecount(),
        "density": lambda g: g.density(),
        "fg_communities": lambda g: g.community_fastgreedy(),
        "community": lambda g: g.community_multilevel(),
        "assortativity": lambda g: g.assortativity_degree(),
        "betweenness": lambda g: g.betweenness(),
        "degree": lambda g: g.degree(),
        "closeness": lambda g: g.closeness(),
        "eigen_centrality": lambda g: g.eigenvector_centrality(),
        "shortest_path_length": lambda g: g.distances(),
        "radius": lambda g: g.radius(),
        "diameter": lambda g: g.diameter(),
        "average_path_length": lambda g: g.average_path_length(),
        "clustering_coefficient": lambda g: g.transitivity_undirected(),
        "jaccard_coefficient": lambda g: g.similarity_jaccard(),
        "pagerank": lambda g: g.pagerank(),
        "fg_modularity": lambda g: g.community_fastgreedy().as_clustering(),
        "average_clustering": lambda g: g.transitivity_avglocal_undirected(),
        "average_local_clustering": lambda g: g.transitivity_local_undirected(),
        "transitivity_undirected": lambda g: g.transitivity_undirected(),
    }

    # Calculate metrics for all graphs passed to class
    def calculate_metrics(self, selected_metrics=None):
        if selected_metrics is None:
            selected_metrics = self.available_metrics.keys()

        calculated_metrics = {}

        for graph_name, graph in self.graph_dict.items():
            metrics_for_graph = {}
            for metric in selected_metrics:
                metric_function = self.available_metrics[metric]
                metrics_for_graph[metric] = metric_function(graph)
            calculated_metrics[graph_name] = metrics_for_graph

        return calculated_metrics

    def get_metrics(self, graph_dict_or_function=None, cell_lines=None, chromosomes=None, resolutions=None, metrics: Union[None, dict, list] = None, root_dir=None):
        """
        Filter metrics from dict, function returning dict or from root directory containing edgelists
        :param cell_lines: cell lines to filter on
        :param graph_dict_or_function: input as dictionary or function returning dictionary
        :param chromosomes: chromosome to filter on
        :param resolutions: specific resolution to filter on
        :param metrics: network metrics to calculate
        :param root_dir: root dir containing edge lists (if calculating metrics from edge lists)
        :return: dict containing graph objects and metrics
        """

        graph_dict = {}

        # If root_dir is provided, create graph_dict from the directory
        if root_dir is not None:
            graph_creator = CreateGraphsFromDirectory(root_dir)
            graph_creator.from_edgelists()
            graph_dict = graph_creator.graph_dict
        else:
            # If graph_dict_or_function is a function, call it to get the graph_dict
            if callable(graph_dict_or_function):
                graph_dict = graph_dict_or_function()
            # If graph_dict_or_function is a dictionary, use it directly
            elif isinstance(graph_dict_or_function, dict):
                graph_dict = graph_dict_or_function

        # Filter the graph_dict based on the filtering criteria
        if cell_lines or chromosomes or resolutions:
            graph_dict = CreateGraphsFromDirectory("").filter_graphs(graph_dict=graph_dict, chromosomes=chromosomes, resolutions=resolutions)

        # Filter the edges within the graph objects based on the chromosome information
        if chromosomes:
            for graph_name, graph in graph_dict.items():
                df = pd.DataFrame([(e.source_vertex['name'], e.target_vertex['name']) for e in graph.es], columns=['source', 'target'])
                df = df[df['source'].str.startswith(tuple(chromosomes)) & df['target'].str.startswith(tuple(chromosomes))]
                graph = ig.Graph.TupleList(df.itertuples(index=False), directed=False)
                graph_dict[graph_name] = graph

        # If metrics is None, use all available metrics
        if metrics is None:
            metrics = cls.available_metrics
        # If metrics is a list, filter the default_metrics dictionary
        elif isinstance(metrics, list):
            metrics = {metric.lower(): cls.available_metrics[metric.lower()] for metric in metrics}

        metrics_data = {}
        for graph_name, graph in graph_dict.items():
            metrics_data[graph_name] = {}
            for metric_name, metric_function in metrics.items():
                metric_value = metric_function(graph)
                metrics_data[graph_name][metric_name] = metric_value

        return metrics_data

    @classmethod
    def print_metrics(cls, metrics_dict, metric_names=None):
        if metric_names is None:
            metric_names = []

        for graph_name, metrics in metrics_dict.items():
            print(f"Metrics for {graph_name}:")
            for metric_name, metric_value in metrics.items():
                if metric_name in metric_names and callable(metric_value):
                    try:
                        metric_value = metric_value()
                    except TypeError:
                        if metric_name == "fg_modularity":
                            membership = metric_value.community_multilevel().membership
                            metric_value = metric_value.modularity(membership)
                        else:
                            raise ValueError(f"Unsupported metric '{metric_name}'")
                print(f"  {metric_name}: {metric_value}")
            print()


graph_dict = mcf7_chr18_1mb()
print("Graphs in the dictionary:")
for name, graph in graph_dict.items():
    print(f"  {name}: {graph.summary()}")

metrics = NetworkMetrics.get_metrics(
    graph_dict_or_function=graph_dict,
    metrics=[
        "size", "edges", "fg_communities"
    ])
NetworkMetrics.print_metrics(metrics)

###################
# Calculate metrics
###################

# From function returning dictionary containing graph objects (or from dictionary):

# MCF10 1MB norm chr18:
# def mcf10_1mb_norm_metrics_chr18():
#     metrics = NetworkMetrics.get_metrics(cell_lines=["mcf10"], chromosomes=["chr18"], resolutions=["1000000"],
#                                          graph_dict_or_function=mcf7_10_lowres_graphs,
#                                          metrics=[
#                                              "size", "edges", "fg_communities"
#                                          ])
#     return NetworkMetrics.print_metrics(metrics)
# mcf10_1mb_norm_metrics_chr18()

# def mcf10_1mb_norm_metrics_chr18():
#     network_metrics = NetworkMetrics(mcf7_10_lowres_graphs)
#     metrics = network_metrics.get_metrics(cell_lines=["mcf7"], chromosomes=["chr18"], resolutions=["1000000"],
#                                           metrics=[
#                                               "size", "edges", "fg_communities"
#                                           ])
#     return network_metrics.print_metrics(metrics)
#
# mcf10_1mb_norm_metrics_chr18()



# # MCF7 1MB norm:
# def mcf7_chr18_1mb_norm_metrics():
#     metrics = NetworkMetrics.get_metrics(
#         graph_dict_or_function=mcf7_1MB_norm_chr18,
#         metrics=[
#             "size", "edges", "fg_communities"
#         ])
#     return NetworkMetrics.print_metrics(metrics)


# mcf7_chr18_1mb_norm_metrics()


# Or calculate metrics from root directory containing edge lists:
# def chr18_size_mod_from_directory():
#     metrics = NetworkMetrics.get_metrics(cell_lines="mcf7", chromosomes="chr18", resolutions="1000000", metrics=["size", "nodes", "fg_communities"], root_dir=Path("/Users/GBS/Master/HiC-Data/edgelists/lowres_mcf7_mcf10"))
#     return NetworkMetrics.print_metrics(metrics)
# chr18_size_mod_from_directory()


class LargestComponent:
    """
    Class to find the LCC from a given graph object.
    Pass any graph obj dict and return largest conected component.
    """

    def __init__(self, graph_dict_or_function):
        if isinstance(graph_dict_or_function, types.FunctionType):
            self.graph_dict = graph_dict_or_function()
        else:
            self.graph_dict = graph_dict_or_function

    def find_lcc(self):
        largest_component_dict = {}
        for graph_name, graph in self.graph_dict.items():
            largest_component = graph.connected_components().giant()  # or g.clusters().giant() or graph.components().giant()?
            largest_component_dict[graph_name] = largest_component
        return largest_component_dict

    def print_lcc(self):
        for graph_name, graph in self.find_lcc().items():
            print(f"LCC for: {graph_name} \n size: {graph.vcount()} \n edges: {graph.ecount()}")


# LargestComponent(chr18_inc_norm_graphs_50kb()).print_lcc()

# def print_50kb_norm():
#     metrics = NetworkMetrics.get_metrics(chr18_inc_norm_graphs_50kb(), metrics=["size", "edges"])
#     return NetworkMetrics.print_metrics(metrics)
# # print_50kb_norm()
#
# LargestComponent(mcf10_1mb_norm_graph()).print_lcc()


# TODO: Find lcc to total size ratio for each graph

# TODO: Plot the lcc to size ratio for each cell line and each resolution as stacked bar plot

# TODO: Plot degree distribution for each cell line and each resolution

# TODO: Calculate betweenness centrality for LCC (on cell line --> resolution level).

# TODO: Calculate closeness centrality for LCC (on cell line --> resolution level).

# TODO:

# TODO: Make a class that does differential community detection for two graphs: Two cell lines on same resolution.

# TODO: Plot the jaccard index for two cell lines on same resolution (for all resolutions).

# TODO: Make Alpaca differential community detection that takes two graphs: Two cell lines on same resolution and calculates the difference in community detection metrics.

# TODO: Plot the differential community detection stuff.


# If I need to use networkx for some reason later:

class ConvertIgraphToNetworkx:

    def __init__(self, graph_dict):
        self.graph_dict = graph_dict
        self.nx_graph_dict = {}

    def convert(self):
        for graph_key in self.graph_dict:
            graph = self.graph_dict[graph_key]

            # Create an empty NetworkX graph
            nx_graph = nx.Graph()

            # Add nodes from the iGraph graph
            nx_graph.add_nodes_from(range(graph.vcount()))

            # Add edges from the iGraph graph
            nx_graph.add_edges_from(graph.get_edgelist())

            # Copy vertex attributes
            for v in graph.vs:
                for attr in v.attributes():
                    nx_graph.nodes[v.index][attr] = v[attr]

            # Copy edge attributes
            for e in graph.es:
                u, v = e.tuple
                for attr in e.attributes():
                    nx_graph.edges[u, v][attr] = e[attr]

            self.nx_graph_dict[graph_key] = nx_graph

        return self.nx_graph_dict

    def __str__(self):
        if not self.nx_graph_dict:
            return "No NetworkX graph generated yet. Please run the convert() method first."

        output_str = ""
        for graph_key, nx_graph in self.nx_graph_dict.items():
            nodes = nx_graph.nodes()
            edges = nx_graph.edges()
            output_str += f"Graph: {graph_key}\nNodes: {nodes}\nEdges: {edges}\n\n"

        return output_str


# ConvertIgraphToNetworkx(chr18_inc_norm_graphs_50kb()).convert()

class PlotNetworkxGraphs:

    def __init__(self, networkx_graph_dict):
        self.networkx_graph_dict = networkx_graph_dict

    def plot(self):
        for graph_key, nx_graph in self.networkx_graph_dict.items():
            plt.figure()
            plt.title(f"{graph_key}")
            nx.draw(nx_graph, with_labels=True, node_color="skyblue", font_weight="bold", node_size=1000)
            plt.show()


# plotter = PlotNetworkxGraphs(ConvertIgraphToNetworkx(chr18_inc_norm_graphs_50kb()).convert())
# plotter.plot()

class InteractivePlotNetworkxGraphs:

    def __init__(self, networkx_graph_dict):
        self.networkx_graph_dict = networkx_graph_dict

    def plot(self):
        for graph_key, nx_graph in self.networkx_graph_dict.items():
            # Use a spring layout to position the nodes
            pos = nx.spring_layout(nx_graph, seed=42)

            # Create node trace
            node_trace = go.Scatter(mode='markers+text', textposition='bottom center',
                                    hoverinfo='text', marker=dict(showscale=False, size=20, colorscale='Viridis', reversescale=True, colorbar=dict(thickness=15, title='Node Connections', xanchor='left', titleside='right')))

            # Create edge trace
            edge_trace = go.Scatter(line=dict(width=0.5, color='#888'), hoverinfo='none', mode='lines')

            node_x, node_y, node_text, node_colors = [], [], [], []
            edge_x, edge_y = [], []

            # Add nodes and edges to the traces
            for node, adjacencies in enumerate(nx_graph.adjacency()):
                x, y = pos[node]
                node_x.append(x)
                node_y.append(y)
                node_colors.append(len(adjacencies[1]))
                node_info = f"{node} - # of connections: {len(adjacencies[1])}"
                node_text.append(node_info)

                for neighbor in adjacencies[1].keys():
                    x0, y0 = pos[node]
                    x1, y1 = pos[neighbor]
                    edge_x.extend([x0, x1, None])
                    edge_y.extend([y0, y1, None])

            node_trace.update(x=node_x, y=node_y, text=node_text, marker=dict(color=node_colors))
            edge_trace.update(x=edge_x, y=edge_y)

            # Customize the plot appearance
            fig = go.Figure(data=[edge_trace, node_trace],
                            layout=go.Layout(title=graph_key, showlegend=False, hovermode='closest',
                                             margin=dict(b=20, l=5, r=5, t=40),
                                             xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                                             yaxis=dict(showgrid=False, zeroline=False, showticklabels=False)))
            fig.show()

# interactive_plotter = InteractivePlotNetworkxGraphs(ConvertIgraphToNetworkx(chr18_inc_norm_graphs_50kb()).convert())
# interactive_plotter.plot()
