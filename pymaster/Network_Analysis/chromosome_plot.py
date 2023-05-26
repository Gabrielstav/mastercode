# Import modules
import pandas as pd
from Graph_Processing import graph_metrics as gm
from Graph_Processing import graph_instances as gi
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pathlib as path
import re


def cytoband_path():
    cytoband = path.Path("/Users/GBS/Master/Reference/hg19/cytoBand_hg19.txt")
    return cytoband

def natural_sort_key(s):
    return [int(text) if text.isdigit() else text for text in re.split(r'(\d+)', s)]


class ChromosomePlot:
    def __init__(self, graph_dict, chromosomes):
        self.graph_dict = graph_dict
        self.chromosomes = self.get_chromosomes(chromosomes)
        self.ideogram_data = self.get_ideogram_data()
        self.node_data = self.get_node_data()

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
                node_data['colors'].append('#2243a8')  # Blue looks good

        return node_data

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

    def plot(self):
        ideo_df = pd.DataFrame.from_dict(self.ideogram_data)
        nodes_df = pd.DataFrame.from_dict(self.node_data)

        chrom_height = 1
        chrom_spacing = 1
        node_height = 0.4
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

        for collection in self.chromosome_collections(nodes_df, node_ybase, node_height, alpha=0.5, linewidths=0):
            ax.add_collection(collection)

        # Convert x-axis labels to megabases
        ticks_loc = ax.get_xticks().tolist()
        ax.xaxis.set_major_locator(ticker.FixedLocator(ticks_loc))
        ax.set_xticklabels([str(int(each/1e6))+'M' for each in ticks_loc])

        ax.set_yticks([chrom_centers[i] for i in self.chromosomes])
        ax.set_yticklabels(self.chromosomes)
        ax.axis('tight')
        plt.show()


def ideogram():
    graphs = gi.all_graphs()
    filter_instance = gm.FilterGraphs(graphs)
    filter_instance.filter_graphs(cell_lines=["gsm2824367"], resolutions=[1000000], interaction_type="intra", condition="intra-split-raw")
    graph_dict = filter_instance.graph_dict
    ideogram_instance = ChromosomePlot(graph_dict, "all")
    ideogram_instance.plot()

ideogram()
