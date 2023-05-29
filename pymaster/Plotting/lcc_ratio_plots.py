# Import modules
import igraph as ig
from Graph_Processing import graph_metrics as gm
from Graph_Processing import graph_instances as gi
import matplotlib.pyplot as plt
import numpy as np


# Make LCC ratio plots

class LCC_ratio_plot:
    def __init__(self, graph_dict):
        self.graph_dict = graph_dict

    def get_lcc_ratio(self):
        inst = gm.LCC_Ratio(self.graph_dict)
        return inst.calculate_lcc_ratio_per_chromosome()

    def plot_stacked_bar(self):
        ratios = self.get_lcc_ratio()

        chromosomes = list(ratios.keys())
        chromosomes.sort(key=lambda x: int(x[3:]) if x.startswith('chr') and x[3:].isdigit() else float('inf'))

        lcc_ratios = []
        remaining_ratios = []

        for chromosome in chromosomes:
            ratio = ratios[chromosome]
            lcc_ratios.append(ratio[0])
            remaining_ratios.append(1 - ratio[0])

        x = np.arange(len(chromosomes))
        width = 0.35

        fig, ax = plt.subplots()
        lcc_bars = ax.bar(x, lcc_ratios, width, label='LCC Ratio')
        remaining_bars = ax.bar(x, remaining_ratios, width, bottom=lcc_ratios, label='Remaining Ratio')

        ax.set_xlabel('Chromosome')
        ax.set_ylabel('Proportion')
        ax.set_title('Proportion of Nodes in LCC per Chromosome')
        ax.set_xticks(x)
        ax.set_xticklabels(chromosomes)
        ax.legend()

        def autolabel(bars):
            for bar in bars:
                height = bar.get_height()
                ax.annotate(f'{height:.2f}', xy=(bar.get_x() + bar.get_width() / 2, height),
                            xytext=(0, 3), textcoords='offset points',
                            ha='center', va='bottom')

        autolabel(lcc_bars)
        autolabel(remaining_bars)

        plt.show()

def plot_lcc():
    graphs = gi.all_graphs()
    filter_instance = gm.FilterGraphs(graphs)
    filter_instance.filter_graphs(cell_lines=["gsm2824367"], resolutions=[1000000], interaction_type="intra", condition="intra-split-raw")
    graph_dict = filter_instance.graph_dict
    lcc_plot_instance = LCC_ratio_plot(graph_dict)
    lcc_plot_instance.plot_stacked_bar()

plot_lcc()


