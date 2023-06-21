# Import modules
import pandas as pd
from Graph_Processing import graph_metrics as gm
from Graph_Processing import graph_instances as gi
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pathlib as path
import re
from Network_Analysis import centrality_metrics as cm
from Graph_Processing import graph_generator as gg

class Directories:
    base_path = path.Path("/Users/GBS/Master/Figures")
    chromosome_path = base_path / "chromosome_plots"

    if not chromosome_path.exists():
        chromosome_path.mkdir(parents=True)



def cytoband_path():
    cytoband = path.Path("/Users/GBS/Master/Reference/hg19/cytoBand_hg19.txt")
    return cytoband

def natural_sort_key(s):
    return [int(text) if text.isdigit() else text for text in re.split(r'(\d+)', s)]


class ChromosomePlotHelper:
    def __init__(self, graph_dict, chromosomes):
        self.graph_dict = graph_dict
        self.chromosomes = self.get_chromosomes(chromosomes)
        self.ideogram_data = self.get_ideogram_data()
        self.node_data = self.get_node_data()
        self.compartment_data = None

    @staticmethod
    def read_compartment_data(bed_file):
        df = pd.read_csv(bed_file, sep='\t', header=None, names=['chrom', 'start', 'end', 'compartment', 'eigenvalue'])
        # Convert compartment column to 'A' or 'B' based on the sign of eigenvalue
        df['compartment'] = df['eigenvalue'].apply(lambda x: 'A' if x >= 0 else 'B')
        # Convert the dataframe to a list of dictionaries and return
        return df.to_dict('records')

    def load_compartment_data(self, bed_file):
        self.compartment_data = self.read_compartment_data(bed_file)

    def get_chromosomes(self, chromosomes):
        if chromosomes == "all":
            all_chromosomes = set()
            for graph_name, graph in self.graph_dict.items():
                for node in graph.vs:
                    all_chromosomes.add(node["chromosome"])
            return sorted(list(all_chromosomes), key=natural_sort_key)
        else:
            return sorted(chromosomes, key=natural_sort_key)

    def get_ideogram_data(self):
        ideo_df = pd.read_table(
            cytoband_path(),
            skiprows=0,
            names=['chrom', 'start', 'end', 'name', 'gieStain']
        )
        ideo_df = ideo_df[ideo_df.chrom.apply(lambda x: x in self.chromosomes)]
        ideo_df['width'] = ideo_df.end - ideo_df.start

        color_lookup = {
            'gneg': (1., 1., 1.),
            'gpos25': (.6, .6, .6),
            'gpos50': (.4, .4, .4),
            'gpos75': (.2, .2, .2),
            'gpos100': (0., 0., 0.),
            'acen': (.8, .4, .4),
            'gvar': (.8, .8, .8),
            'stalk': (.9, .9, .9),
        }
        ideo_df['colors'] = ideo_df['gieStain'].apply(lambda x: color_lookup[x])

        return ideo_df.to_dict('list')

    def get_node_data(self):
        node_data = {
            'chrom': [],
            'start': [],
            'end': [],
            'colors': [],
        }

        for graph_name, graph in self.graph_dict.items():
            for node in graph.vs:
                chrom = node["chromosome"]
                start, end = map(int, node['location'].split('-'))
                node_data['chrom'].append(chrom)
                node_data['start'].append(start)
                node_data['end'].append(end)
                node_data['colors'].append("blue")  # "('#2243a8')  # Blue looks good (change to green?)

        return node_data

class LinearChromosomePlot:
    def __init__(self, plot_helper):
        self.plot_helper = plot_helper
        self.chromosomes = plot_helper.chromosomes

    @staticmethod
    def chromosome_collections(df, y_positions, height, alpha=1.0, linewidths=1.0, edgecolor='black'):
        if 'width' not in df.columns:
            df['width'] = df['end'] - df['start']
        for chrom, group in df.groupby('chrom'):
            if chrom not in y_positions:
                continue
            yrange = (y_positions[chrom], height)
            xranges = group[['start', 'width']].values
            yield plt.broken_barh(
                xranges, yrange, facecolors=group['colors'], alpha=alpha, linewidths=linewidths, edgecolor=edgecolor)

    @staticmethod
    def compartment_collections(df, y_positions, height, alpha=0.5, linewidths=0):
        """
        Add a collection for compartments. Compartment data (df) should include 'chrom', 'start', 'end', and 'compartment' columns,
        with 'compartment' being 'A' or 'B'. Compartment 'A' is colored red, and 'B' is colored blue.
        """
        df['width'] = df['end'] - df['start']
        df['colors'] = df['compartment'].map({'A': 'red', 'B': 'blue'})
        for chrom, group in df.groupby('chrom'):
            if chrom not in y_positions:
                continue
            yrange = (y_positions[chrom], height)
            xranges = group[['start', 'width']].values
            yield plt.broken_barh(
                xranges, yrange, facecolors=group['colors'], alpha=alpha, linewidths=linewidths)

    def plot(self, save_as=None):
        ideo_df = pd.DataFrame.from_dict(self.plot_helper.ideogram_data)
        nodes_df = pd.DataFrame.from_dict(self.plot_helper.node_data)
        compartments_df = pd.DataFrame.from_dict(self.plot_helper.compartment_data)

        chrom_height = 1
        chrom_spacing = 1
        node_height = 0.5
        node_padding = 0.1
        figsize = (6, 8)

        ybase = 0
        chrom_ybase = {}
        node_ybase = {}
        chrom_centers = {}

        for chrom in self.chromosomes[::-1]:
            chrom_ybase[chrom] = ybase
            chrom_centers[chrom] = ybase + chrom_height / 2.
            node_ybase[chrom] = ybase - node_height - node_padding
            ybase += chrom_height + chrom_spacing

        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111)

        for collection in self.chromosome_collections(ideo_df, chrom_ybase, chrom_height):
            ax.add_collection(collection)

        for collection in self.chromosome_collections(nodes_df, node_ybase, node_height, alpha=0.7, linewidths=0):
            ax.add_collection(collection)

        if self.plot_helper.compartment_data is not None:
            compartment_df = pd.DataFrame.from_dict(self.plot_helper.compartment_data)
            for collection in self.compartment_collections(compartment_df, chrom_ybase, chrom_height, alpha=0.4, linewidths=0):
                ax.add_collection(collection)

        # Convert x-axis labels to megabases
        ticks_loc = ax.get_xticks().tolist()
        ax.xaxis.set_major_locator(ticker.FixedLocator(ticks_loc))
        ax.set_xticklabels([str(int(each/1e6))+'M' for each in ticks_loc])

        ax.set_yticks([chrom_centers[i] for i in self.chromosomes])
        ax.set_yticklabels(self.chromosomes)
        ax.axis('tight')
        if save_as is not None:
            plt.savefig(save_as, dpi=300)
        plt.show()


def linear_plot():
    graphs = gi.all_graphs()
    filter_instance = gm.FilterGraphs(graphs)
    filter_instance.filter_graphs(cell_lines=["mcf7"], resolutions=[1000000], interaction_type="intra", condition="intra-split-norm")  # , chromosomes=["chr1", "chr2"])
    graph_dict = filter_instance.graph_dict
    plot_helper = ChromosomePlotHelper(graph_dict, "all")  # , ["chr1", "chr2"])
    lin_plot = LinearChromosomePlot(plot_helper)
    lin_plot.plot(save_as=Directories.chromosome_path / "mcf7_chrom.png")

# linear_plot()


def compartment_plot():
    graphs = gi.all_graphs()
    filter_instance = gm.FilterGraphs(graphs)
    filter_instance.filter_graphs(cell_lines=["imr90"], resolutions=[1000000], interaction_type="intra", condition="intra-split-raw")  # , chromosomes=["chr1", "chr2"])
    graph_dict = filter_instance.graph_dict
    plot_helper = ChromosomePlotHelper(graph_dict, "all")  # , ["chr1", "chr2"])

    # load compartment data
    compartment_file = "/Users/GBS/Master/HiC-Data/compartments/imr90_1Mb_multiples_compartments.bed"
    plot_helper.load_compartment_data(compartment_file)

    lin_plot = LinearChromosomePlot(plot_helper)
    lin_plot.plot()

# compartment_plot()

# circ plot with sorted chromosomes (just need to somehow add the ideogram data to the circ plot):
def plot_degree_network_circ():
    graphs = gi.all_graphs()
    filter_instance = gm.FilterGraphs(graphs)
    filter_instance.filter_graphs(cell_lines=["mcf7"], resolutions=[1000000], interaction_type="inter", condition="inter-nosplit-raw")  # , chromosomes=["chr1"])
    graph_dict = filter_instance.graph_dict
    sorter_instance = gg.Sorter(graph_dict)
    sorted_graph = sorter_instance.sort_graph()
    degree_instance = cm.PlotDegreeNetwork(sorted_graph)  # Use LCC for withtin-chromosome plots
    degree_instance.calculate_degree()
    degree_instance.normalize_degree()
    degree_instance.plot_graph(save_as=None, normalize=False, color_edges=True, layout="circle")

# plot_degree_network_circ()


if __name__ == "__main__":
    print("RUN")

