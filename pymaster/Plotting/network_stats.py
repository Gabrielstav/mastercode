# Import modules
import igraph as ig
from Graph_Processing import graph_metrics as gm
from Graph_Processing import graph_instances as gi
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import pandas as pd
import pathlib as path
from Graph_Processing import graph_metrics as gm


class Directories:
    base_path = path.Path("/Users/GBS/Master/Figures")
    networkstats_path = base_path / "networks_stats"

    if not networkstats_path.exists():
        networkstats_path.mkdir(parents=True)

# TODO: Make Interactions per resolution and size per resolution

# TODO: Make connceted components per resolution

# TODO: Make LCC ratio and proportion to all nodes per resolution

# Size per network

# Clean this up later, its ugly and inefficient, make helper functions and combine all plot functions into one

class GetMetrics:
    def __init__(self, graph_dict):
        self.graph_dict = graph_dict

    def get_interactions(self):
        size_interactions = {}
        for name, graph in self.graph_dict.items():
            size_interactions[(graph['cell_line'], graph['norm_status'], graph['resolution'])] = (graph.ecount())
        return size_interactions

    def get_nodes(self):
        node_number = {}
        for name, graph in self.graph_dict.items():
            node_number[(graph['cell_line'], graph['norm_status'], graph['resolution'])] = (graph.vcount())
        return node_number

    def get_components(self):
        component_number = {}
        for name, graph in self.graph_dict.items():
            connected_components = graph.components()
            component_number[(graph['cell_line'], graph['norm_status'], graph['resolution'])] = len(connected_components)
            print(f"(Graph name: {name} number of components: {len(connected_components)}")
        return component_number

    def get_lcc(self):
        lcc_number = {}
        for name, graph in self.graph_dict.items():
            lcc_number[(graph['cell_line'], graph['norm_status'], graph['resolution'])] = (graph.components().giant().vcount())
        return lcc_number

    def get_lcc_ratio(self):
        lcc_ratio = {}
        for name, graph in self.graph_dict.items():
            lcc_ratio[(graph['cell_line'], graph['norm_status'], graph['resolution'])] = (graph.components().giant().vcount()/graph.vcount())
        return lcc_ratio


    # TODO: Maybe
    def get_component_size_distribution(self):
        component_size_distributions = {}
        for name, graph in self.graph_dict.items():
            connected_components = graph.components()
            component_size_distributions[(graph['cell_line'], graph['norm_status'], graph['resolution'])] = [len(component) for component in connected_components]
        return component_size_distributions

    def plot_size(self, save_as=None):
        size_interactions = self.get_nodes()

        unique_keys = set(size_interactions.keys())
        cell_lines = sorted(list(set(key[0] for key in unique_keys)))
        norm_statuses = sorted(list(set(key[1] for key in unique_keys)))  # Collect all normalization statuses
        resolutions = sorted(list(set(key[2] for key in unique_keys)))  # Resolution is now at index 2
        color_map = cm.get_cmap("jet", len(cell_lines)*len(norm_statuses))  # Multiply by the number of norm statuses


        color_index_map = {}  # Map for unique color index
        color_index = 0
        max_color_index = len(cell_lines) * len(norm_statuses) - 1
        for cell_line in cell_lines:
            for norm_status in norm_statuses:
                color_index_map[(cell_line, norm_status)] = color_index
                color_index += 1

        for cell_line in cell_lines:
            for norm_status in norm_statuses:  # Iterate over norm statuses
                y_values = []
                x_values = []
                linestyle = "--" if norm_status == "norm" else "-"  # Choose linestyle based on norm status
                for resolution in resolutions:
                    key = (cell_line, norm_status, resolution)  # Key now includes norm status
                    if key in size_interactions:
                        y_values.append(int(size_interactions[key]))
                        x_values.append(resolution)

                # Set labels and legend differently for normalized mcf7 and mcf10 data
                if norm_status == "norm" and cell_line in ["mcf7", "mcf10"]:
                    normalized_color_index = color_index_map[(cell_line, norm_status)] / max_color_index
                    plt.plot(x_values, y_values, marker="o", linestyle=linestyle, color=color_map(normalized_color_index), label=f"{cell_line} (ICE)")  # ICE label
                elif norm_status == "raw":
                    normalized_color_index = color_index_map[(cell_line, norm_status)] / max_color_index
                    plt.plot(x_values, y_values, marker="o", linestyle=linestyle, color=color_map(normalized_color_index), label=f"{cell_line}")

        plt.xticks(resolutions, [f"{int(i/1000)}" if i < 1000000 else "1000" for i in resolutions])
        plt.xlabel("Resolution (kb)")
        plt.ylabel("Number of Nodes")
        plt.title("Number of nodes per network")
        plt.tick_params(axis='x', labelsize=6)
        plt.tick_params(axis='y', labelsize=8)
        plt.legend()
        if save_as:
            plt.savefig(save_as, dpi=600, format='png')
        plt.show()

    def plot_interactions(self, save_as=None):
        plt.clf()  # Clear the current figure
        size_interactions = self.get_interactions()

        unique_keys = set(size_interactions.keys())
        cell_lines = sorted(list(set(key[0] for key in unique_keys)))
        norm_statuses = sorted(list(set(key[1] for key in unique_keys)))
        resolutions = sorted(list(set(key[2] for key in unique_keys)))
        color_map = cm.get_cmap('jet', len(cell_lines) * len(norm_statuses))  # Multiply by the number of norm statuses to get unique colors

        color_index_map = {}  # Map for unique color index
        color_index = 0
        max_color_index = len(cell_lines) * len(norm_statuses) - 1
        for cell_line in cell_lines:
            for norm_status in norm_statuses:
                color_index_map[(cell_line, norm_status)] = color_index
                color_index += 1

        for cell_line in cell_lines:
            for norm_status in norm_statuses:
                y_values = []
                x_values = []
                linestyle = "--" if norm_status == "norm" else "-"
                for resolution in resolutions:
                    key = (cell_line, norm_status, resolution)
                    if key in size_interactions:
                        if key not in size_interactions:
                            print(key)
                        y_values.append(np.log10(size_interactions[key]))
                        x_values.append(resolution)

                # Set labels and legend differently for normalized mcf7 and mcf10 data
                if norm_status == "norm" and cell_line in ["mcf7", "mcf10"]:
                    normalized_color_index = color_index_map[(cell_line, norm_status)] / max_color_index
                    plt.plot(x_values, y_values, marker="o", linestyle=linestyle, color=color_map(normalized_color_index), label=f"{cell_line} (ICE)")  # ICE label
                elif norm_status == "raw":
                    normalized_color_index = color_index_map[(cell_line, norm_status)] / max_color_index
                    plt.plot(x_values, y_values, marker="o", linestyle=linestyle, color=color_map(normalized_color_index), label=f"{cell_line}")


        plt.xticks(resolutions, [f"{int(i / 1000)}" if i < 1000000 else "1000" for i in resolutions])
        plt.xlabel("Resolution (kb)")
        plt.ylabel("Number of Interactions (log10)")
        plt.title("Number of Interactions per network")
        plt.tick_params(axis='x', labelsize=6)
        plt.tick_params(axis='y', labelsize=8)
        plt.legend()
        if save_as:
            plt.savefig(save_as, dpi=300, format='png')
        plt.show()

    def plot_number_components(self, save_as=None):
        connected_components = self.get_components()
        unique_keys = set(connected_components.keys())
        cell_lines = sorted(list(set(key[0] for key in unique_keys)))
        norm_statuses = sorted(list(set(key[1] for key in unique_keys)))
        resolutions = sorted(list(set(key[2] for key in unique_keys)))
        color_map = cm.get_cmap("jet", len(cell_lines) * len(norm_statuses))

        color_index_map = {}  # Map for unique color index
        color_index = 0
        for cell_line in cell_lines:
            for norm_status in norm_statuses:
                color_index_map[(cell_line, norm_status)] = color_index
                color_index += 1

        max_color_index = len(cell_lines) * len(norm_statuses) - 1

        for cell_line in cell_lines:
            for norm_status in norm_statuses:
                y_values = []
                x_values = []
                linestyle = "--" if norm_status == "norm" else "-"
                for resolution in resolutions:
                    key = (cell_line, norm_status, resolution)
                    if key in connected_components:
                        y_values.append(np.log10(connected_components[key]))
                        x_values.append(resolution)

                # Set labels and legend differently for normalized mcf7 and mcf10 data
                normalized_color_index = color_index_map[(cell_line, norm_status)] / max_color_index
                if norm_status == "norm" and cell_line in ["mcf7", "mcf10"]:
                    plt.plot(x_values, y_values, marker="o", linestyle=linestyle, color=color_map(normalized_color_index), label=f"{cell_line} (ICE)")  # ICE label
                elif norm_status == "raw":
                    plt.plot(x_values, y_values, marker="o", linestyle=linestyle, color=color_map(normalized_color_index), label=f"{cell_line}")

        plt.xticks(resolutions, [f"{int(i / 1000)}" if i < 1000000 else "1000" for i in resolutions])
        plt.xlabel("Resolution (kb)")
        plt.ylabel("Number of Connected Components (log10)")
        plt.title("Number of connected components per network")
        plt.tick_params(axis='x', labelsize=6)
        plt.tick_params(axis='y', labelsize=8)
        plt.legend()
        if save_as:
            plt.savefig(save_as, dpi=300, format='png')
        plt.show()

    def plot_lcc_size(self, save_as=None):
        lcc_size = self.get_lcc()
        unique_keys = set(lcc_size.keys())
        cell_lines = sorted(list(set(key[0] for key in unique_keys)))
        norm_statuses = sorted(list(set(key[1] for key in unique_keys)))
        resolutions = sorted(list(set(key[2] for key in unique_keys)))
        color_map = cm.get_cmap("jet", len(cell_lines) * len(norm_statuses))

        color_index_map = {}  # Map for unique color index
        color_index = 0
        for cell_line in cell_lines:
            for norm_status in norm_statuses:
                color_index_map[(cell_line, norm_status)] = color_index
                color_index += 1

        max_color_index = len(cell_lines) * len(norm_statuses) - 1

        for cell_line in cell_lines:
            for norm_status in norm_statuses:  # Iterate over norm statuses
                y_values = []
                x_values = []
                linestyle = "--" if norm_status == "norm" else "-"  # Choose linestyle based on norm status
                for resolution in resolutions:
                    key = (cell_line, norm_status, resolution)  # Key now includes norm status
                    if key in lcc_size:
                        y_values.append(lcc_size[key])
                        x_values.append(resolution)

                # Set labels and legend differently for normalized mcf7 and mcf10 data
                normalized_color_index = color_index_map[(cell_line, norm_status)] / max_color_index
                if norm_status == "norm" and cell_line in ["mcf7", "mcf10"]:
                    plt.plot(x_values, y_values, marker="o", linestyle=linestyle, color=color_map(normalized_color_index), label=f"{cell_line} (ICE)")  # ICE label
                elif norm_status == "raw":
                    plt.plot(x_values, y_values, marker="o", linestyle=linestyle, color=color_map(normalized_color_index), label=f"{cell_line}")

        plt.xticks(resolutions, [f"{int(i / 1000)}" if i < 1000000 else "1000" for i in resolutions])
        plt.xlabel("Resolution (kb)")
        plt.ylabel("Number of Connected Components")
        plt.title("Number of connected components per network")
        plt.tick_params(axis='x', labelsize=6)
        plt.tick_params(axis='y', labelsize=8)
        plt.legend()
        if save_as:
            plt.savefig(save_as, dpi=300, format='png')
        plt.show()

    def plot_lcc_ratio(self, save_as=None):
        lcc_ratio = self.get_lcc_ratio()
        unique_keys = set(lcc_ratio.keys())
        cell_lines = sorted(list(set(key[0] for key in unique_keys)))
        norm_statuses = sorted(list(set(key[1] for key in unique_keys)))
        resolutions = sorted(list(set(key[2] for key in unique_keys)))
        color_map = cm.get_cmap("jet", len(cell_lines) * len(norm_statuses))

        color_index_map = {}

        color_index = 0
        for cell_line in cell_lines:
            for norm_status in norm_statuses:
                color_index_map[(cell_line, norm_status)] = color_index
                color_index += 1

        max_color_index = len(cell_lines) * len(norm_statuses) - 1

        for cell_line in cell_lines:

            for norm_status in norm_statuses:
                y_values = []
                x_values = []
                linestyle = "--" if norm_status == "norm" else "-"
                for resolution in resolutions:
                    key = (cell_line, norm_status, resolution)
                    if key in lcc_ratio:
                        y_values.append(lcc_ratio[key])
                        x_values.append(resolution)

                normalized_color_index = color_index_map[(cell_line, norm_status)] / max_color_index
                if norm_status == "norm" and cell_line in ["mcf7", "mcf10"]:
                    plt.plot(x_values, y_values, marker="o", linestyle=linestyle, color=color_map(normalized_color_index), label=f"{cell_line} (ICE)")
                elif norm_status == "raw":
                    plt.plot(x_values, y_values, marker="o", linestyle=linestyle, color=color_map(normalized_color_index), label=f"{cell_line}")

        plt.xticks(resolutions, [f"{int(i / 1000)}" if i < 1000000 else "1000" for i in resolutions])
        plt.xlabel("Resolution (kb)")
        plt.ylabel("LCC Ratio")
        plt.title("LCC Ratio per network")
        plt.tick_params(axis='x', labelsize=6)
        plt.tick_params(axis='y', labelsize=8)
        plt.legend()
        if save_as:
            plt.savefig(save_as, dpi=300, format='png')
        plt.show()

    def plot_component_size_distribution(self, save_as=None):
        component_size_distribution = self.get_component_size_distribution()
        unique_keys = set(component_size_distribution.keys())
        cell_lines = sorted(list(set(key[0] for key in unique_keys)))
        norm_statuses = sorted(list(set(key[1] for key in unique_keys)))
        resolutions = sorted(list(set(key[2] for key in unique_keys)))
        color_map = cm.get_cmap("jet", len(cell_lines) * len(norm_statuses))

        color_index_map = {}

        color_index = 0
        for cell_line in cell_lines:
            for norm_status in norm_statuses:
                color_index_map[(cell_line, norm_status)] = color_index
                color_index += 1

        max_color_index = len(cell_lines) * len(norm_statuses) - 1

        for cell_line in cell_lines:

            for norm_status in norm_statuses:
                y_values = []
                x_values = []
                linestyle = "--" if norm_status == "norm" else "-"
                for resolution in resolutions:
                    key = (cell_line, norm_status, resolution)
                    if key in component_size_distribution:
                        y_values.append(component_size_distribution[key])
                        x_values.append(resolution)

                normalized_color_index = color_index_map[(cell_line, norm_status)] / max_color_index
                if norm_status == "norm" and cell_line in ["mcf7", "mcf10"]:
                    plt.plot(x_values, y_values, marker="o", linestyle=linestyle, color=color_map(normalized_color_index), label=f"{cell_line} (ICE)")
                elif norm_status == "raw":
                    plt.plot(x_values, y_values, marker="o", linestyle=linestyle, color=color_map(normalized_color_index), label=f"{cell_line}")

        plt.xticks(resolutions, [f"{int(i / 1000)}" if i < 1000000 else "1000" for i in resolutions])
        plt.xlabel("Resolution (kb)")
        plt.ylabel("Number of Components")
        plt.title("Component Size Distribution")
        plt.legend()
        if save_as:
            plt.savefig(save_as, dpi=300, format='png')
        plt.show()


def get_size():
    graphs = gi.all_graphs()
    filter_instance = gm.FilterGraphs(graphs)
    filter_instance.filter_graphs(interaction_type="intra", condition="intra-split-raw")
    filter_instance2 = gm.FilterGraphs(graphs)
    filter_instance2.filter_graphs(cell_lines=["mcf10", "mcf7"], condition="intra-split-norm")
    filter_instance.graph_dict.update(filter_instance2.graph_dict)
    graph_dict = filter_instance.graph_dict
    size = GetMetrics(graph_dict)
    size.plot_size(save_as= Directories.networkstats_path / "intra_size_plot.png")
# get_size()

def get_interactions():
    graphs = gi.all_graphs()
    filter_instance = gm.FilterGraphs(graphs)
    filter_instance.filter_graphs(interaction_type="intra", condition="intra-split-raw")
    filter_instance2 = gm.FilterGraphs(graphs)
    filter_instance2.filter_graphs(cell_lines=["mcf10", "mcf7"], condition="intra-split-norm")
    filter_instance.graph_dict.update(filter_instance2.graph_dict)
    graph_dict = filter_instance.graph_dict
    size = GetMetrics(graph_dict)
    size.plot_interactions(save_as= Directories.networkstats_path / "intra_interactions_plot.png")
# get_interactions()

def get_components():
    graphs = gi.all_graphs()
    filter_instance = gm.FilterGraphs(graphs)
    filter_instance.filter_graphs(interaction_type="intra", condition="intra-split-raw")
    filter_instance2 = gm.FilterGraphs(graphs)
    filter_instance2.filter_graphs(cell_lines=["mcf10", "mcf7"], condition="intra-split-norm")
    filter_instance.graph_dict.update(filter_instance2.graph_dict)
    graph_dict = filter_instance.graph_dict
    size = GetMetrics(graph_dict)
    size.plot_number_components(save_as= Directories.networkstats_path / "intra_components_plot.png")
# get_components()

def get_lcc():
    graphs = gi.all_graphs()
    filter_instance = gm.FilterGraphs(graphs)
    filter_instance.filter_graphs(interaction_type="intra", condition="intra-split-raw")
    filter_instance2 = gm.FilterGraphs(graphs)
    filter_instance2.filter_graphs(cell_lines=["mcf10", "mcf7"], condition="intra-split-norm")
    filter_instance.graph_dict.update(filter_instance2.graph_dict)
    graph_dict = filter_instance.graph_dict
    size = GetMetrics(graph_dict)
    size.plot_lcc_size(save_as= Directories.networkstats_path / "intra_lcc_plot.png")
# get_lcc()

def get_lcc_ratio():
    graphs = gi.all_graphs()
    filter_instance = gm.FilterGraphs(graphs)
    filter_instance.filter_graphs(interaction_type="intra", condition="intra-split-raw")
    filter_instance2 = gm.FilterGraphs(graphs)
    filter_instance2.filter_graphs(cell_lines=["mcf10", "mcf7"], condition="intra-split-norm")
    filter_instance.graph_dict.update(filter_instance2.graph_dict)
    graph_dict = filter_instance.graph_dict
    size = GetMetrics(graph_dict)
    size.plot_lcc_ratio(save_as= Directories.networkstats_path / "intra_lcc_ratio_plot.png")

# get_lcc_ratio()

# TODO: Maybe
def component_size_distribution():
    graphs = gi.all_graphs()
    filter_instance = gm.FilterGraphs(graphs)
    filter_instance.filter_graphs(interaction_type="intra", condition="intra-split-raw")
    filter_instance2 = gm.FilterGraphs(graphs)
    filter_instance2.filter_graphs(cell_lines=["mcf10", "mcf7"], condition="intra-split-norm")
    filter_instance.graph_dict.update(filter_instance2.graph_dict)
    graph_dict = filter_instance.graph_dict
    size = GetMetrics(graph_dict)
    size.plot_component_size_distribution(save_as= None)  # Directories.networkstats_path / "intra_component_size_distribution_plot.png")

# component_size_distribution()


if __name__ == "__main__":
    print("RUN")












