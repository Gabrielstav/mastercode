# Import modules
import igraph as ig
import networkx as nx
from Graph_Processing import graph_metrics as Ga
from matplotlib import pyplot as plt
from pathlib import Path

ig.config["plotting.backend"] = "matplotlib"
# Default (?) and deprecated backend:
# ig.config["plotting.backend"] = "cairo"

# Classes for plotting networks


class Settings:

    def __init__(self):
        self.root_dir = Path("/Users/GBS/Master/Figures/iGraph")

    def get_output_dir(self):
        return self.root_dir


# TODO: Make plotting class that takes any graph dict from any class in Network_Metrics and plots it as a network


class plot_graph:

    def __init__(self, graph_dict, output_dir):
        self.graph_dict = graph_dict
        self.output_dir = output_dir

    @staticmethod
    def abbreviate_label(self, label, resolution):
        chrom, region = label.split(':')
        start, end = region.split('-')

        if resolution >= 1000000:
            start_unit = int(start) // 1000000
            end_unit = int(end) // 1000000
            unit = "MB"
        elif resolution >= 1000:
            start_unit = int(start) // 1000
            end_unit = int(end) // 1000
            unit = "KB"
        else:
            start_unit = int(start)
            end_unit = int(end)
            unit = "B"

        # Can be used to abel graphs like this, not useful for large graphs:
        # resolution = int(re.search(r"_([^_]+)_", graph_name).group(1))
        # abbreviated_labels = [self.abbreviate_label(label, resolution) for label in graph.vs["name"]]
        return f"{chrom}:{start_unit}-{end_unit} {unit}"

    def show_graph(self):
        print(self.graph_dict)
        for graph_name, graph in self.graph_dict.items():
            ig.plot(graph)
            plt.show()

    def show_graph_with_lcc(self, largest_component_obj):
        for graph_name, graph in self.graph_dict.items():
            lcc_membership = largest_component_obj.lcc_membership()[graph_name]
            node_colors = ["red" if membership == 1 else "black" for membership in lcc_membership]
            edge_colors = ["red" if lcc_membership[edge.source] == 1 and lcc_membership[edge.target] else "black" for edge in graph.es]
            ig.plot(graph, vertex_color=node_colors, edge_color=edge_colors)
            plt.show()

    def show_only_lcc(self, largest_component_obj):
        for graph_name, graph in self.graph_dict.items():
            lcc_membership = largest_component_obj.lcc_membership()[graph_name]
            lcc_graph = graph.subgraph([node.index for node in graph.vs if lcc_membership[node.index] == 1])
            ig.plot(lcc_graph, vertex_color="red", edge_color="black")
            plt.show()

    def save_graph(self):
        for graph_name, graph in self.graph_dict.items():
            visual_style = {"vertex_size": 1, "vertex_label": graph.vs["name"], "layout": graph.layout("kk")}
            # Remove special characters from graph name to use as file name:
            safe_graph_name = "".join(e for e in graph_name if e.isalnum() or e == "_")
            output_filename = self.output_dir / f"{safe_graph_name}.png"
            ig.plot(graph, output_filename)  # **visual_style)
            print(f"Saved plot to {output_filename}")

# def plot_imr90_chr18_1mb

# def plot_full():
#     dir_manager = Settings()
#     output_dir = dir_manager.get_output_dir()
#     plot = plot_graph(imr90_chr18(), output_dir)
#     plot.show_graph()
# # plot_full()
#
# def plot_lcc():
#     dir_manager = Settings()
#     output_dir = dir_manager.get_output_dir()
#     graph_dict = mcf7_chr18_1mb()
#     plot = plot_graph(graph_dict, output_dir)
#     largest_component_obj = Gp.LargestComponent(graph_dict)
#     plot.show_graph_with_lcc(largest_component_obj)
# # plot_lcc()
#
# def plot_only_lcc():
#     dir_manager = Settings()
#     output_dir = dir_manager.get_output_dir()
#     graph_dict = imr90_chr18()
#     plot = plot_graph(graph_dict, output_dir)
#     largest_component_obj = Gp.LargestComponent(graph_dict)
#     plot.show_only_lcc(largest_component_obj)
# plot_only_lcc()


# Seems redundant, gets almost same results as LCC membership:
# def plot_cc():
#     degree = 1
#     size = 2
#     graph_dict = imr90_chr18()
#
#     filter_components = NM.FilterComponents(graph_dict, degree)
#     filtered_graphs = filter_components.filter_on_size(size)
#
#     dir_manager = SetDirectories()
#     output_dir = dir_manager.get_output_dir()
#     plot = plot_graph(filtered_graphs, output_dir)
#     plot.show_graph()
#
# plot_cc()

# TODO: Make stacked bar plots showing ratio of nodes in LCC to total nodes for each cell line and each resolution, takes graph dicts as input.




class plot_lcc_ratio:
    def __init__(self, graph_dict, output_dir):
        self.graph_dict = graph_dict
        self.output_dir = output_dir

    def plot_lcc_ratio_bar(self):
        lcc_ratio_calculator = Ga.LCC_Ratio(self.graph_dict)
        lcc_ratio_dict = lcc_ratio_calculator.calculate_lcc_ratio()

        for graph_name, graph_sizes in lcc_ratio_dict.items():
            ratio = graph_sizes[1] / graph_sizes[0]
            plt.bar(graph_name, ratio)

        plt.title("LCC ratio")
        plt.xlabel("Graph")
        plt.ylabel("Ratio")
        plt.show()

    def plot_lcc_ratio_per_chromosome(self):
        lcc_ratio_calc = Ga.LCC_Ratio(self)
        lcc_ratio_dict = lcc_ratio_calc.calculate_lcc_ratio_per_chromosome()

        for graph_name, graph_sizes in lcc_ratio_dict.items():
            ratio = graph_sizes[1] / graph_sizes[0]
            plt.bar(graph_name, ratio)

        plt.title("LCC ratio")
        plt.xlabel("Graph")
        plt.ylabel("Ratio")
        plt.show()


# def plot_imr90_lcc_ratios():
#     graph_dict = imr90_graphs()
#     plot_lcc_ratio.plot_lcc_ratio_per_chromosome(graph_dict)

# plot_imr90_lcc_ratios()

# def print_lcc_ratios():
#     graph_dict = imr90_graphs()
#     lcc_ratio_object = Ga.LCC_Ratio(graph_dict)
#     lcc_ratio_dict = lcc_ratio_object.calculate_lcc_ratio_per_chromosome()
#     print(lcc_ratio_dict)
#
# print_lcc_ratios()


# def plot_lcc_ratio_imr90():
#     dir_manager = Settings()
#     output_dir = dir_manager.get_output_dir()
#     graph_dict = imr90_graphs()
#     plot = plot_lcc_ratio(graph_dict, output_dir)
#     plot.plot_lcc_ratio_bar()
# # plot_lcc_ratio_imr90()



# TODO: Make degree distribution plots for each cell lines LCC for each chromosome, per resolution? Takes graph dict as input.

# TODO: Make betweenness plot that colors the nodes by betwenness centrality, takes a graph dict as input. Make

# TODO: Make closeness plot that colors the nodes by closeness, takes graph dict as input.



# Networkx plot

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








