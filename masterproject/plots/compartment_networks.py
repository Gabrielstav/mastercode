# Import modules
import pathlib as path
import igraph as ig
import matplotlib.pyplot as plt
from graph_processing import graph_metrics as gm
from graph_processing import graph_instances as gi
import numpy as np
import pandas as pd
import seaborn as sns


def compartment_paths(cell_line=None):

    comp_dir = path.Path("/Users/GBS/Master/HiC-Data/mcf_comps")
    comp_path = None

    if cell_line == "mcf10":
        comp_path = (comp_dir / "mcf10_compartments_noIC_no0.bed")

    elif cell_line == "mcf7":
        comp_path = (comp_dir / "mcf7_compartments_noIC_no0.bed")

    else:
        print("Please specify valid cell line (mcf10 or mcf7")

    return comp_path

class plot_dirs:
    base_path = path.Path("/Users/GBS/Master/Figures")
    comp_path = path.Path("/Users/GBS/Master/Figures/compartments")

    if not comp_path.exists():
        comp_path.mkdir(parents=True)


class CompartmentNetworkAnnotator:
    def __init__(self, bed_file, graph_dict):
        self.bed_file = bed_file
        self.graph_dict = graph_dict

    def annotate_networks(self):
        # Load BED file and create a dictionary to store compartments
        compartments = {}
        with open(self.bed_file, 'r') as f:
            for line in f:
                chrom, start, end, compartment, _ = line.strip().split('\t')
                compartments[(chrom, int(start), int(end))] = compartment

        # Iterate over the graph dictionary and annotate nodes
        for graph_name, graph in self.graph_dict.items():
            compartments_attr = [None] * graph.vcount()  # Initialize compartment attribute list

            for edge in graph.es:
                source = edge.source
                target = edge.target

                # Check if source node overlaps with any bin in the BED file
                for bin_range, compartment in compartments.items():
                    chrom, start, end = bin_range
                    if chrom == graph.vs[source]['chromosome']:
                        bin_start, bin_end = map(int, graph.vs[source]['location'].split('-'))
                        if start <= bin_start <= end or start <= bin_end <= end:
                            compartments_attr[source] = compartment
                            break

                # Check if target node overlaps with any bin in the BED file
                for bin_range, compartment in compartments.items():
                    chrom, start, end = bin_range
                    if chrom == graph.vs[target]['chromosome']:
                        bin_start, bin_end = map(int, graph.vs[target]['location'].split('-'))
                        if start <= bin_start <= end or start <= bin_end <= end:
                            compartments_attr[target] = compartment
                            break

            # Add compartment attributes to the network nodes
            graph.vs['compartment'] = compartments_attr

        return self.graph_dict

    def print_edgelist(self):
        for graph_name, graph in self.graph_dict.items():
            print(f"\nGraph: {graph_name}")
            for edge in graph.es:
                source = graph.vs[edge.source]
                target = graph.vs[edge.target]
                source_name = source['name']
                target_name = target['name']
                source_compartment = source['compartment']
                target_compartment = target['compartment']
                print(f"{source_name} ({source_compartment}) -- {target_name} ({target_compartment})")


class CompartmentNetworkPlotter:
    def __init__(self, graph_dict):
        self.graph_dict = graph_dict


    # @staticmethod
    # def plot_legend(graph, top_n):
    #     # Get top_n nodes with the highest degree
    #     top_nodes = sorted(graph.vs, key=lambda v: v["degree"], reverse=True)[:top_n]
    #
    #     # Generate colors for the legend
    #     colors = [plt.cm.viridis(node["degree"] / max(graph.degree())) for node in top_nodes]
    #
    #     # Change labels to Mega base pairs (Mb) format without decimal places and append degree
    #     labels = ["-".join(str(int(pos)//10**6) + 'Mb' for pos in node["location"].split('-')) + ': ' + str(node["degree"]) for node in top_nodes]
    #
    #     # Create legend handles
    #     legend_elements = [Line2D([0], [0], marker='o', color='w', label=label,
    #                               markerfacecolor=color, markersize=10) for color, label in zip(colors, labels)]
    #
    #     # Position the legend to avoid overlap with the network
    #     # Change the values in the tuple to adjust the position
    #     plt.legend(handles=legend_elements, loc=(1, 0))
    #
    #     # Return the legend handles
    #     return legend_elements

    def plot_networks(self, save_as=None):
        for graph_name, graph in self.graph_dict.items():

            # Create a dictionary to map compartments to colors
            compartment_colors = {"A": "red", "B": "blue"}

            # Create a list to store node colors
            node_colors = []

            # Iterate over nodes and assign colors based on compartment
            for node in graph.vs:
                compartment = node["compartment"]
                color = compartment_colors.get(compartment, "black")
                node_colors.append(color)

            # Set visual style
            fig, ax = plt.subplots()
            graph_name = graph_name.split("_")[1:2]
            node_size = 30 / graph.vcount()
            edge_width = 60 / graph.ecount()
            visual_style = {
                "layout": "kk",
                "vertex_size": node_size,
                "edge_width": edge_width,
                "bbox": (6000, 6000),
                "margin": 200,
                "vertex_color": node_colors,
                "vertex_label_size": 10,
                "vertex_label_dist": 10,
                "vertex_label_angle": 100,
                "vertex_label_color": "black",
                "vertex_frame_color": node_colors,  # Remove black outline, use in circle layout
                # "vertex_label": graph.vs["name"],  # Add location as label, can be used in small graphs
            }

            # Plot the graph
            ig.plot(graph, **visual_style, target=ax)
            plt.title(f"{graph_name}")

            if save_as:
                plt.savefig(save_as, dpi=300)

            plt.show()




def annotate_network():
    # Load and filter graphs
    graphs = gi.all_graphs()
    filter_instance = gm.FilterGraphs(graphs)
    filter_instance.filter_graphs(cell_lines=["mcf10"], interaction_type="intra", condition="intra-split-raw", resolutions=[1000000], chromosomes=["chr1"])
    graph_dict = filter_instance.graph_dict

    # Annotate networks with compartments
    annotater = CompartmentNetworkAnnotator(compartment_paths("mcf10"), graph_dict)
    annotated_graphs = annotater.annotate_networks()
    annotater.graph_dict = annotated_graphs  # Update the graph_dict attribute with the annotated graphs
    annotater.print_edgelist()

    # Plot annotated network
    # plot_instance = CompartmentNetworkPlotter(annotated_graphs)
    # plot_instance.plot_networks(save_as=plot_dirs.comp_path / "mcf10_chr16_compartments")

annotate_network()


# TODO: Make compartment switch plot for mcf10 and mcf7
#   Make network plot of same nodes that switched compartments, make the the nodes that went B to A red and A to B blue, no switch or overlap = gray? Or square? Kind of redundant..
#   And calcualte percentage switched per chromosome for overlapping nodes

class CompartmentComparison:
    def __init__(self, graph_dict_reference, graph_dict_perturbed):
        self.graph_ref = list(graph_dict_reference.values())[0]
        self.graph_change = list(graph_dict_perturbed.values())[0]

    def compare_compartments(self):
        no_overlap_count = 0
        no_change_count = 0
        active_count = 0
        inactive_count = 0

        for node_ref in self.graph_ref.vs:
            node_location_ref = node_ref['location']
            compartment_ref = node_ref['compartment']

            corresponding_nodes_change = [node for node in self.graph_change.vs if node['location'] == node_location_ref]

            if len(corresponding_nodes_change) == 0:
                no_overlap_count += 1
            else:
                corresponding_node_change = corresponding_nodes_change[0]
                compartment_change = corresponding_node_change['compartment']

                if compartment_change == compartment_ref:
                    no_change_count += 1
                elif compartment_ref == 'B' and compartment_change == 'A':
                    active_count += 1
                elif compartment_ref == 'A' and compartment_change == 'B':
                    inactive_count += 1

        print("No Overlap: ", no_overlap_count)
        print("No Change: ", no_change_count)
        print("Changed Active: ", active_count)
        print("Changed Inactive: ", inactive_count)

        return no_overlap_count, no_change_count, active_count, inactive_count

    def compare_compartments_chromosome(self):
        results = {}

        for node_ref in self.graph_ref.vs:
            node_location_ref = node_ref['location']
            compartment_ref = node_ref['compartment']
            chromosome = node_ref['chromosome']

            if chromosome not in results:
                results[chromosome] = {'No Overlap': 0, 'No Change': 0, 'Active': 0, 'Inactive': 0}

            corresponding_nodes_change = [node for node in self.graph_change.vs if node['location'] == node_location_ref]

            if len(corresponding_nodes_change) == 0:
                results[chromosome]['No Overlap'] += 1
            else:
                corresponding_node_change = corresponding_nodes_change[0]
                compartment_change = corresponding_node_change['compartment']

                if compartment_change == compartment_ref:
                    results[chromosome]['No Change'] += 1
                elif compartment_ref == 'B' and compartment_change == 'A':
                    results[chromosome]['Active'] += 1
                elif compartment_ref == 'A' and compartment_change == 'B':
                    results[chromosome]['Inactive'] += 1

        return results

    def plot_comparison(self, save_as=None):
        labels = ['No Overlap', 'No Change', 'Active', 'Inactive']
        counts = self.compare_compartments()

        x = np.arange(len(labels))
        plt.bar(x, counts)
        plt.xticks(x, labels)
        plt.xlabel('Compartments')
        plt.ylabel('Count')
        plt.title('Comparison of Compartment Changes')
        if save_as:
            plt.savefig(save_as, dpi=300)
        plt.show()

    def plot_comparison_chromosome(self, save_as=None):
        results = self.compare_compartments_chromosome()

        labels = ['No Change', 'Active', 'Inactive']
        colors = ['gray', 'red', 'blue']

        chromosomes = sorted(results.keys(), key=lambda x: int(x[3:]) if x[3:].isdigit() else float('inf'))
        chromosome_labels = [chrom.replace("chr", "") for chrom in chromosomes]  # Remove "chr" from chromosome name

        # Remove 'No Overlap' from the results
        for chromosome in results:
            results[chromosome].pop('No Overlap', None)

        x = np.arange(len(chromosomes))  # The label locations
        width = 0.25  # The width of the bars

        fig, ax = plt.subplots()

        # Create bars for each label
        for i, label in enumerate(labels):
            ax.bar(x - width / 2 + i * width, [results[chromosome][label] for chromosome in chromosomes], width, label=label, color=colors[i])

        # Add some text for labels, title, and custom x-axis tick labels, etc.
        ax.set_xlabel('Chromosome')
        ax.set_ylabel('Node Count')
        ax.set_title('Compartment Change from MCF-10A to MCF-7')
        ax.set_xticks(x)
        ax.set_xticklabels(chromosome_labels)
        ax.legend()

        fig.tight_layout()
        if save_as:
            plt.savefig(save_as, dpi=300)
        plt.show()

    def plot_comparison_chromosome_stacked(self, save_as=None):
        results = self.compare_compartments_chromosome()
        chromosomes = sorted(results.keys(), key=lambda x: int(x[3:]) if x[3:].isdigit() else float('inf'))
        labels = ['No Change', 'Active', 'Inactive']
        label_colors = {'No Change': 'gray', 'Active': 'red', 'Inactive': 'blue'}

        x = np.arange(len(chromosomes))
        width = 0.8

        bottom = np.zeros(len(chromosomes))
        for label in labels:
            counts = [results[chromosome][label] for chromosome in chromosomes]
            plt.bar(x, counts, width, bottom=bottom, label=label, color=label_colors[label])
            bottom += np.array(counts)

        plt.xticks(x, [chromosome[3:] for chromosome in chromosomes], rotation=45)
        plt.xlabel('Chromosome')
        plt.ylabel('Node Count')
        plt.title('Compartment Change from MCF-10A to MCF-7')
        plt.legend()

        if save_as:
            plt.savefig(save_as, dpi=300)
        plt.show()




def compare_compartments():

    # Load and filter reference and perturbed graphs
    graphs_ref = gi.all_graphs()
    filter_instance_ref = gm.FilterGraphs(graphs_ref)
    filter_instance_ref.filter_graphs(cell_lines=["mcf10"], interaction_type="intra", condition="intra-split-raw", resolutions=[1000000])  # , chromosomes=["chr1"])
    reference_graph = filter_instance_ref.graph_dict

    graphs_per = gi.all_graphs()
    filter_instance_perturbed = gm.FilterGraphs(graphs_per)
    filter_instance_perturbed.filter_graphs(cell_lines=["mcf7"], interaction_type="intra", condition="intra-split-raw", resolutions=[1000000])  # , chromosomes=["chr1"])
    perturbed_graph = filter_instance_perturbed.graph_dict

    # Call compartments
    annotater_ref = CompartmentNetworkAnnotator(compartment_paths("mcf10"), reference_graph)
    annotated_graphs_ref = annotater_ref.annotate_networks()
    annotater_ref.graph_dict = annotated_graphs_ref  # Update the graph_dict attribute with the annotated graphs

    annotater_per = CompartmentNetworkAnnotator(compartment_paths("mcf7"), perturbed_graph)
    annotated_graphs_per = annotater_per.annotate_networks()
    annotater_per.graph_dict = annotated_graphs_per

    # Compare compartments
    comparison = CompartmentComparison(reference_graph, perturbed_graph)
    comparison.plot_comparison_chromosome(save_as=plot_dirs.comp_path / "mcf10_mcf10_unstacked.png")

# compare_compartments()


def compare_compartments_overall():

    # Load and filter reference and perturbed graphs
    graphs_ref = gi.all_graphs()
    filter_instance_ref = gm.FilterGraphs(graphs_ref)
    filter_instance_ref.filter_graphs(cell_lines=["mcf10"], interaction_type="intra", condition="intra-split-raw", resolutions=[1000000])  # , chromosomes=["chr1"])
    reference_graph = filter_instance_ref.graph_dict

    graphs_per = gi.all_graphs()
    filter_instance_perturbed = gm.FilterGraphs(graphs_per)
    filter_instance_perturbed.filter_graphs(cell_lines=["mcf7"], interaction_type="intra", condition="intra-split-raw", resolutions=[1000000])  # , chromosomes=["chr1"])
    perturbed_graph = filter_instance_perturbed.graph_dict

    # Call compartments
    annotater_ref = CompartmentNetworkAnnotator(compartment_paths("mcf10"), reference_graph)
    annotated_graphs_ref = annotater_ref.annotate_networks()
    annotater_ref.graph_dict = annotated_graphs_ref  # Update the graph_dict attribute with the annotated graphs

    annotater_per = CompartmentNetworkAnnotator(compartment_paths("mcf7"), perturbed_graph)
    annotated_graphs_per = annotater_per.annotate_networks()
    annotater_per.graph_dict = annotated_graphs_per

    # Compare compartments
    comparison = CompartmentComparison(reference_graph, perturbed_graph)
    # comparison.plot_comparison(save_as=plot_dirs.comp_path / "mcf10_mcf10_overall_compartments_comparison")

# compare_compartments_overall()

class DegreeCompartmentComparison:

    def __init__(self, reference_graph, perturbed_graph):
        self.graph_reference = list(reference_graph.values())[0]
        self.graph_perturbed = list(perturbed_graph.values())[0]

    def plot(self, save_as=None):
        # Extract degrees and compartment status from each graph
        degrees_ref = [(node['degree'], node['compartment']) for node in self.graph_reference.vs]
        degrees_per = [(node['degree'], node['compartment']) for node in self.graph_perturbed.vs]

        # Convert to DataFrame
        df_ref = pd.DataFrame(degrees_ref, columns=['Degree', 'Compartment'])
        df_per = pd.DataFrame(degrees_per, columns=['Degree', 'Compartment'])

        # Add 'Graph' column
        df_ref['Graph'] = 'mcf10'
        df_per['Graph'] = 'mcf7'

        # Combine 'Graph' and 'Compartment' into one column
        df_ref['Legend'] = df_ref['Graph'] + ' ' + df_ref['Compartment']
        df_per['Legend'] = df_per['Graph'] + ' ' + df_per['Compartment']

        # Concatenate dataframes
        df = pd.concat([df_ref, df_per])

        # Define a color palette
        palette = {"mcf10 A": "red", "mcf10 B": "blue", "mcf7 A": "darkred", "mcf7 B": "darkblue"}

        # Plot
        plt.figure(figsize=(10, 8))

        for legend_group in df['Legend'].unique():
            data = df[df['Legend'] == legend_group]['Degree'].values
            hist, bin_edges = np.histogram(data, bins=50)
            bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
            plt.plot(bin_centers, hist, label=legend_group, color=palette[legend_group])

        plt.xscale('log')
        plt.yscale('log')
        plt.title('Degree Distribution by Compartment Status')
        plt.xlabel('Degree')
        plt.ylabel('Frequency')
        plt.legend(title='Compartment Status')

        if save_as:
            plt.savefig(save_as, dpi=300)
        plt.show()

class ClosenessCompartmentComparison:

    def __init__(self, reference_graph, perturbed_graph):
        self.graph_reference = list(reference_graph.values())[0]
        self.graph_perturbed = list(perturbed_graph.values())[0]

    def plot(self, save_as=None):
        # Extract degrees and compartment status from each graph
        degrees_ref = [(node['closeness'], node['compartment']) for node in self.graph_reference.vs]
        degrees_per = [(node['closenss'], node['compartment']) for node in self.graph_perturbed.vs]

        # Combine data and convert to pandas DataFrame for easier plotting
        degrees_all = degrees_ref + degrees_per
        df = pd.DataFrame(degrees_all, columns=['Closeness', 'Compartment'])

        # Plot degree distribution using seaborn for automatic histogram and kde
        plt.figure(figsize=(10, 8))
        sns.histplot(df, x='Degree', hue='Compartment', kde=True, stat="density", common_norm=False)

        plt.title('Degree Distribution by Compartment')
        plt.xlabel('Degree')
        plt.ylabel('Density')

        if save_as:
            plt.savefig(save_as, dpi=300)

        plt.show()


class BetweennessCompartmentComparison:

    def __init__(self, reference_graph, perturbed_graph):
        self.graph_reference = list(reference_graph.values())[0]
        self.graph_perturbed = list(perturbed_graph.values())[0]

    def plot(self, save_as=None):
        # Extract degrees and compartment status from each graph
        degrees_ref = [(node['betweenness'], node['compartment']) for node in self.graph_reference.vs]
        degrees_per = [(node['betweenness'], node['compartment']) for node in self.graph_perturbed.vs]

        # Combine data and convert to pandas DataFrame for easier plotting
        degrees_all = degrees_ref + degrees_per
        df = pd.DataFrame(degrees_all, columns=['Closeness', 'Compartment'])

        # Plot degree distribution using seaborn for automatic histogram and kde
        plt.figure(figsize=(10, 8))
        sns.histplot(df, x='Degree', hue='Compartment', kde=True, stat="density", common_norm=False)

        plt.title('Degree Distribution by Compartment')
        plt.xlabel('Degree')
        plt.ylabel('Density')

        if save_as:
            plt.savefig(save_as, dpi=300)

        plt.show()


def degree_distribution_of_compartment_nodes():

    # Load and filter reference and perturbed graphs
    graphs_ref = gi.all_graphs()
    filter_instance_ref = gm.FilterGraphs(graphs_ref)
    filter_instance_ref.filter_graphs(cell_lines=["mcf10"], interaction_type="intra", condition="intra-split-raw", resolutions=[1000000])  # , chromosomes=["chr1"])
    reference_graph = filter_instance_ref.graph_dict

    graphs_per = gi.all_graphs()
    filter_instance_perturbed = gm.FilterGraphs(graphs_per)
    filter_instance_perturbed.filter_graphs(cell_lines=["mcf7"], interaction_type="intra", condition="intra-split-raw", resolutions=[1000000])  # , chromosomes=["chr1"])
    perturbed_graph = filter_instance_perturbed.graph_dict

    # Call compartments
    annotater_ref = CompartmentNetworkAnnotator(compartment_paths("mcf10"), reference_graph)
    annotated_graphs_ref = annotater_ref.annotate_networks()
    annotater_ref.graph_dict = annotated_graphs_ref  # Update the graph_dict attribute with the annotated graphs

    annotater_per = CompartmentNetworkAnnotator(compartment_paths("mcf7"), perturbed_graph)
    annotated_graphs_per = annotater_per.annotate_networks()
    annotater_per.graph_dict = annotated_graphs_per

    for graph_name_ref, graph_ref in annotated_graphs_ref.items():
        graph_ref.vs['degree'] = graph_ref.degree()
        graph_ref.vs['closeness'] = graph_ref.closeness()
        graph_ref.vs['betweenness'] = graph_ref.betweenness()

    for graph_name_per, graph_per in annotated_graphs_per.items():
        graph_per.vs['degree'] = graph_per.degree()
        graph_per.vs['closeness'] = graph_per.closeness()
        graph_per.vs['betweenness'] = graph_per.betweenness()

    degree_distribution = DegreeCompartmentComparison(annotated_graphs_ref, annotated_graphs_per)
    degree_distribution.plot(save_as=plot_dirs.comp_path / "mcf10_mcf7_log_degree_compartments_comparison.png")

# degree_distribution_of_compartment_nodes()







