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
        chromosomes.sort(
            key=lambda x: (int(x[3:]) if x.startswith('chr') and x[3:].isdigit() else float('inf'), x)
        )  # Sort chromosomes numerically with special cases
        lcc_ratios = np.array(list(ratios.values()))

        # Create stacked bar plot
        plt.figure()
        plt.bar(chromosomes, lcc_ratios[:, 0], label='Parent Graph Size')
        for i in range(1, lcc_ratios.shape[1]):
            plt.bar(chromosomes, lcc_ratios[:, i], bottom=np.sum(lcc_ratios[:, :i], axis=1),
                    label=f'LCC Ratio {i}')
        plt.xlabel('Chromosome')
        plt.ylabel('LCC Ratio')
        plt.title('LCC Ratios per Chromosome')
        plt.legend()
        plt.show()

def plot_lcc():
    graphs = gi.all_graphs()
    filter_instance = gm.FilterGraphs(graphs)
    filter_instance.filter_graphs(cell_lines=["gsm2824367"], resolutions=[1000000], interaction_type="intra", condition="intra-split-raw")
    graph_dict = filter_instance.graph_dict
    lcc_plot_instance = LCC_ratio_plot(graph_dict)
    lcc_plot_instance.plot_stacked_bar()

plot_lcc()


