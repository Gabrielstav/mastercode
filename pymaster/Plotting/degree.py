


# This works for combining graphs and filtering on them, but source identity is not preserved

@staticmethod
def combine_graphs(graphs):
    all_vertex_names = set()
    for graph in graphs:
        all_vertex_names.update(graph.vs["name"])

    combined_graph = ig.Graph()
    combined_graph.add_vertices(list(all_vertex_names))

    for graph in graphs:
        edge_list = [(graph.vs[e.source]['name'], graph.vs[e.target]['name']) for e in graph.es]
        combined_graph.add_edges(edge_list)

    return combined_graph


def combine_intra_inter_graphs(self):
    filtered_graph_dicts = []
    for graph_dict, filter_params in zip(self.graph_dicts, self.filter_parameters):
        chromosomes = filter_params.get("chromosomes", self.chromosomes)
        resolutions = filter_params.get("resolutions", self.resolutions)

        if chromosomes or resolutions:
            graph_filter = FilterGraphs(graph_dict)
            filtered_graph_dict = graph_filter.filter_graphs(chromosomes=chromosomes, resolutions=resolutions)
        else:
            filtered_graph_dict = graph_dict
        filtered_graph_dicts.append(filtered_graph_dict)

    print("Filtered graph dicts:", filtered_graph_dicts)  # Debugging

    graph_names = set.intersection(*(set(graph_dict.keys()) for graph_dict in filtered_graph_dicts))

    print("Graph names:", graph_names)  # Debugging

    for graph_name in graph_names:
        graphs_to_combine = [graph_dict[graph_name] for graph_dict in filtered_graph_dicts]
        combined_graph = self.combine_graphs(graphs_to_combine)
        self.combined_graphs[graph_name] = combined_graph

    return self.combined_graphs