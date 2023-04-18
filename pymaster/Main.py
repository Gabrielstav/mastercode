# Use fit hic instead of NCHG?
# Find papaer and read it.
# import os.path

# Old implementation of find files:

# @staticmethod
# def find_files(*root_directories):
#     raw_subdirectory_name = "raw"
#     iced_subdirectory_name = "iced"
#     bedfiles = []
#     matrixfiles = []
#     iced_matrixfiles = []
#     resolutions_provided = SetDirectories.get_resolutions()
#
#     # Find the raw data subdirectory in the root directory
#     raw_subdirectories = []
#     for root_directory in root_directories:
#         for root, _, _ in os.walk(root_directory):
#             if os.path.basename(root) == raw_subdirectory_name:
#                 raw_subdirectories.append(root)
#
#     # Find the ICE-normalized data subdirectory in the root directory
#     iced_subdirectories = []
#     for root_directory in root_directories:
#         for root, _, _ in os.walk(root_directory):
#             if os.path.basename(root) == iced_subdirectory_name:
#                 iced_subdirectories.append(root)
#
#     # Helper function to filter files based on resolution
#     def filter_files_on_resolution(files):
#         filtered_files = []
#         for file in files:
#             file_name = os.path.basename(file)
#             resolution_match = re.search(r'_(\d+)[_.]', file_name)
#             if resolution_match:
#                 resolution = int(resolution_match.group(1))
#                 if resolutions_provided is None or resolution in resolutions_provided:
#                     filtered_files.append(file)
#         return filtered_files
#
#     # Recursively search raw data subdirectory for bed and matrix files
#     for subdirectory_path in raw_subdirectories:
#         for root, _, files in os.walk(subdirectory_path):
#             bed_files = [os.path.join(root, file) for file in files if file.endswith(".bed")]
#             matrix_files = [os.path.join(root, file) for file in files if file.endswith(".matrix")]
#             bedfiles.extend(filter_files_on_resolution(bed_files, resolutions_provided))
#             # Debug line
#             print(f"Bedfiles found: {bedfiles}")
#             matrixfiles.extend(filter_files_on_resolution(matrix_files, resolutions_provided))
#             # Debug line
#             print(f"Matrixfiles found: {matrixfiles}")
#
#     # Recursively search ICE-normalized data subdirectory for matrix files
#     for subdirectory_path in iced_subdirectories:
#         for root, _, files in os.walk(subdirectory_path):
#             iced_matrix_files = [os.path.join(root, file) for file in files if file.endswith(".matrix")]
#             iced_matrixfiles.extend(filter_files_on_resolution(iced_matrix_files, resolutions_provided))
#             # Debug line
#             print(f"Iced Matrixfiles found: {iced_matrixfiles}")
#
#     return bedfiles, matrixfiles, iced_matrixfiles

# Old implementation of run nchg:

# old implementation
# nchg_run = sp.run(nchg_command, capture_output=True, check=True)
# print(f"Finished processing file: {bedpe_file}, PID: {os.getpid()}, TID: {threading.get_ident()}")
# return nchg_run.stdout.decode("utf-8").split("\t")

# Old implementation of find_files:

# @staticmethod
# def find_files(*root_directories):
#     """
#     Finds all files bed and matrix files in the raw data subdirectory of the root directory.
#     :param root_directories: one or more root directories to search in
#     :param resolutions: list of resolutions to search for
#     :return: a list of file paths for each BED and matrix file found
#     """
#
#     raw_subdirectory_name = "raw"
#     iced_subdirectory_name = "iced"
#     bedfiles = []
#     matrixfiles = []
#     iced_matrixfiles = []
#     resolutions_provided = SetDirectories.get_resolutions()
#
#     # Find the raw data subdirectory in the root directory
#     raw_subdirectories = []
#     for root_directory in root_directories:
#         for root, _, _ in os.walk(root_directory):
#             if os.path.basename(root) == raw_subdirectory_name:
#                 raw_subdirectories.append(root)
#
#     # Old implementation
#     # Recursively search raw data subdirectory for bed and matrix files
#     for subdirectory_path in raw_subdirectories:
#         for root, _, files in os.walk(subdirectory_path):
#             for file in files:
#                 if file.endswith(".bed"):
#                     bedfiles.append(os.path.join(root, file))
#                 if file.endswith(".matrix"):
#                     matrixfiles.append(os.path.join(root, file))
#
#     # Find the ICE-normalized data subdirectory in the root directory
#     iced_subdirectories = []
#     for root_directory in root_directories:
#         for root, _, _ in os.walk(root_directory):
#             if os.path.basename(root) == iced_subdirectory_name:
#                 iced_subdirectories.append(root)
#
#     # Old implementation
#     # # Recursively search ICE-normalized data subdirectory for matrix files
#     for subdirectory_path in iced_subdirectories:
#         for root, _, files in os.walk(subdirectory_path):
#             for file in files:
#                 if file.endswith(".matrix"):
#                     iced_matrixfiles.append(os.path.join(root, file))
#
#     return bedfiles, matrixfiles, iced_matrixfiles









# Set the path to the Fit-Hi-C executable
# fit_hic_executable = "path/to/fitHiC.py"
#
# # Define the input files and parameters
# hicpro_contact_matrix = "path/to/hicPro_contact_matrix.txt"
# fit_hic_output = "path/to/fit_hic_output.txt"
# resolution = 20000  # Set the desired resolution, e.g., 20 kb

# Run Fit-Hi-C with the specified parameters'
# subprocess.run([
#     "python",
#     fit_hic_executable,
#     "-f", hicpro_contact_matrix,
#     "-o", fit_hic_output,
#     "-r", str(resolution),
#     # Add any other required parameters
# ])

# Chromosome parallelization for NCHG:

# with open(file in input_dir) as input_files:

# chr_dict = {}
# or list with comprehension? tuple? set? queue? stack? tree? graph? array? Think dict.
# for file in input_files:
# find fast way to filter to chromosomes, regex? What is the best way to grep for changes in chr string, if col[0] == chr1?
# Can also filter on inter/intra here but whY?
# Then find way to split the file, use temp files? The input to nchg has to be files infortunately.
# The big thing again though is concatenating the output files into one file (per resolution), before passing the whole genome file to padj.

# Research first: Split files into temp files based on chromosome (fastest), base temp file name on chromosome + file anme (res and exp).
# The new files can be output to a new directory, dict or filter etc to make them per chromosome, and these files are then passed to the NCHG method by just changin the input dir.
# The NCHG method doesn't need to change, jus the input to nchg, since we need to take the input files, and then split the, by chrom, then write to new dir, and then save the output of nchg to new dir,
# and then concatenate the output files into one file (per resolution).

# Pass these files to the NCHG method
# Concatenate the output files into one file (per resolution), in the correct order.

# How does this work with the parallelization? Can we have thread and one process pool for the method? The NCHG method needs to be called with process pool.
# Try to implement process pool for "find siginificant interactions" method.


# nchg_input_file = os.path.abspath("/Users/GBS/Master/HiC-Data/testing_chrom_parallelization/output/temp_dir/input_to_nchg/chr18_inc_20000_no_blacklist_no_cytobands_split/chr18_inc_20000_no_blacklist_no_cytobands_chr18.bedpe")
# concatted_file = os.path.abspath("")
#
# def splitter(input_file):
#     with open(input_file) as input_file:
#         for line in input_file:
#             print(line.split())
#
# splitter(nchg_input_file)

split_chrom_nchg_file = os.path.abspath("/Users/GBS/Master/HiC-Data/testing_chrom_parallelization/eval/split_nchg.txt")

def find_sigs(nchg_out_file):
    sig_pvals = []
    num_interchrom = 0
    num_intrachrom = 0
    inter_list = []
    intra_list = []
    with open(nchg_out_file, "r") as nchg_file:
        for line in nchg_file:
            col = line.split()
            if float(col[6]) < 0.05:
                if col[0] != col[3]:
                    num_interchrom += 1
                    sig_pvals.append(line)
                    inter_list.append(line)
                else:
                    intra_list.append(line)
                    sig_pvals.append(line)
                    num_intrachrom += 1


    return print(f"Number of sig interactions total:, {len(sig_pvals)} in {sig_pvals}. \n"
                 f"Of these, {num_intrachrom} are intrachromosomal, and {num_interchrom} are interchromosmal. \n")

# find_sigs(split_chrom_nchg_file)

nchg_out_file = os.path.abspath("/Users/GBS/Master/HiC-Data/testing_chrom_parallelization/output/temp_dir/NCHG_output/chr18_inc_1000000_nchg_output.txt")


def find_cols(file):
    with open(file, "r") as nchg_file:
        for line in nchg_file:
            col = line.split()
            print(len(col))

find_cols(nchg_out_file)
# find_cols(split_chrom_nchg_file)


# Stable-ish version of input to nchg

# @staticmethod
# def input_to_nchg():
#     """
#     Calls the NCHG script on all files in the no_cytobands directory to find significant interactions
#     """
#
#     no_cytobands_dir_path = SetDirectories.get_temp_dir() + "/no_cytobands"
#     no_cytobands_dir = os.listdir(no_cytobands_dir_path)
#
#     # Create the output directory if it doesn't exist
#     output_dir = os.path.join(SetDirectories.get_temp_dir(), "NCHG_output")
#     if not os.path.exists(output_dir):
#         os.mkdir(output_dir)
#     else:
#         shutil.rmtree(output_dir)
#         os.mkdir(output_dir)
#
#     # Create a directory to store split chromosome files
#     chr_split_base_dir = os.path.join(SetDirectories.get_temp_dir(), "input_to_nchg")
#     if not os.path.exists(chr_split_base_dir):
#         os.mkdir(chr_split_base_dir)
#     else:
#         shutil.rmtree(chr_split_base_dir)
#         os.mkdir(chr_split_base_dir)
#
#     # Split input files by chromosome
#     all_chr_files = []
#     # keep original input file name for each chromosome file
#     input_file_map = {}
#
#     for file in no_cytobands_dir:
#         full_path = os.path.join(no_cytobands_dir_path, file)
#         if not os.path.isfile(full_path):
#             print(f"File {file} does not exist.")
#
#         # Create a subdirectory for each input file inside the chr_split_base_dir
#         file_split_dir = os.path.join(chr_split_base_dir, f"{file[:-len('.bedpe')]}_split")
#         if not os.path.exists(file_split_dir):
#             os.mkdir(file_split_dir)
#         else:
#             shutil.rmtree(file_split_dir)
#             os.mkdir(file_split_dir)
#
#         chr_files = Pipeline.split_bedpe_by_chromosome(full_path, file_split_dir)
#
#         for chr_file in chr_files:
#             input_file_map[chr_file] = file
#         all_chr_files.extend(chr_files)
#
#     # Run find_significant_interactions on chromosome-specific files in parallel
#     output_file_data = defaultdict(list)
#     with concurrent.futures.ProcessPoolExecutor(max_workers=SetDirectories.get_threads()) as executor:
#         futures = list(executor.map(Pipeline.find_significant_interactions, all_chr_files))
#         for bedpe_file, future in zip(all_chr_files, futures):
#             try:
#                 nchg_output = future
#                 input_file = input_file_map[bedpe_file]
#                 output_file_data[input_file].append(nchg_output)
#             except Exception as e:
#                 tid = threading.get_ident()
#                 print(f"Error processing {bedpe_file}: {e}, PID: {os.getpid()}, TID: {tid}")
#
#     # Merge output files back together
#     for input_file, nchg_outputs in output_file_data.items():
#         output_filename = f"{os.path.basename(input_file)[:-len('_no_blacklist_no_cytobands.bedpe')]}_nchg_output.txt"
#         output_filepath = os.path.join(output_dir, output_filename)
#         with open(output_filepath, "w") as f:
#             for nchg_output in nchg_outputs:
#                 f.writelines(line + "\n" for line in nchg_output)


# class NetworkMetrics:
#     """
#     Class to calculate network metrics for a given graph object.
#     """
#
#     def __init__(self, graph_dict_or_function, metrics=None):
#         self.graph_dict = graph_dict
#         # if isinstance(graph_dict_or_function, types.FunctionType):
#         #     self.graph_dict = graph_dict_or_function()
#         # else:
#         #     self.graph_dict = graph_dict_or_function
#         # self.metrics_dict = metrics if metrics is not None else {}
#
#     available_metrics = {
#         "size": lambda g: g.vcount(),
#         "edges": lambda g: g.ecount(),
#         "density": lambda g: g.density(),
#         "fg_communities": lambda g: g.community_fastgreedy(),
#         "community": lambda g: g.community_multilevel(),
#         "assortativity": lambda g: g.assortativity_degree(),
#         "betweenness": lambda g: g.betweenness(),
#         "degree": lambda g: g.degree(),
#         "closeness": lambda g: g.closeness(),
#         "eigen_centrality": lambda g: g.eigenvector_centrality(),
#         "shortest_path_length": lambda g: g.distances(),
#         "radius": lambda g: g.radius(),
#         "diameter": lambda g: g.diameter(),
#         "average_path_length": lambda g: g.average_path_length(),
#         "clustering_coefficient": lambda g: g.transitivity_undirected(),
#         "jaccard_coefficient": lambda g: g.similarity_jaccard(),
#         "pagerank": lambda g: g.pagerank(),
#         "fg_modularity": lambda g: g.community_fastgreedy().as_clustering(),
#         "average_clustering": lambda g: g.transitivity_avglocal_undirected(),
#         "average_local_clustering": lambda g: g.transitivity_local_undirected(),
#         "transitivity_undirected": lambda g: g.transitivity_undirected(),
#     }
#
#     # Calculate metrics for all graphs passed to class
#     def calculate_metrics(self, selected_metrics=None):
#         if selected_metrics is None:
#             selected_metrics = self.available_metrics.keys()
#
#         calculated_metrics = {}
#
#         for graph_name, graph in self.graph_dict.items():
#             metrics_for_graph = {}
#             for metric in selected_metrics:
#                 metric_function = self.available_metrics[metric]
#                 metrics_for_graph[metric] = metric_function(graph)
#             calculated_metrics[graph_name] = metrics_for_graph
#
#         return calculated_metrics
#
#     def get_metrics(self, graph_dict_or_function=None, cell_lines=None, chromosomes=None, resolutions=None, metrics: Union[None, dict, list] = None, root_dir=None):
#         """
#         Filter metrics from dict, function returning dict or from root directory containing edgelists
#         :param cell_lines: cell lines to filter on
#         :param graph_dict_or_function: input as dictionary or function returning dictionary
#         :param chromosomes: chromosome to filter on
#         :param resolutions: specific resolution to filter on
#         :param metrics: network metrics to calculate
#         :param root_dir: root dir containing edge lists (if calculating metrics from edge lists)
#         :return: dict containing graph objects and metrics
#         """
#
#         graph_dict = {}
#
#         # If root_dir is provided, create graph_dict from the directory
#         if root_dir is not None:
#             graph_creator = CreateGraphsFromDirectory(root_dir)
#             graph_creator.from_edgelists()
#             graph_dict = graph_creator.graph_dict
#         else:
#             # If graph_dict_or_function is a function, call it to get the graph_dict
#             if callable(graph_dict_or_function):
#                 graph_dict = graph_dict_or_function()
#             # If graph_dict_or_function is a dictionary, use it directly
#             elif isinstance(graph_dict_or_function, dict):
#                 graph_dict = graph_dict_or_function
#
#         # Filter the graph_dict based on the filtering criteria
#         if cell_lines or chromosomes or resolutions:
#             graph_dict = CreateGraphsFromDirectory("").filter_graphs(graph_dict=graph_dict, chromosomes=chromosomes, resolutions=resolutions)
#
#         # Filter the edges within the graph objects based on the chromosome information
#         if chromosomes:
#             for graph_name, graph in graph_dict.items():
#                 df = pd.DataFrame([(e.source_vertex['name'], e.target_vertex['name']) for e in graph.es], columns=['source', 'target'])
#                 df = df[df['source'].str.startswith(tuple(chromosomes)) & df['target'].str.startswith(tuple(chromosomes))]
#                 graph = ig.Graph.TupleList(df.itertuples(index=False), directed=False)
#                 graph_dict[graph_name] = graph
#
#         # If metrics is None, use all available metrics
#         if metrics is None:
#             metrics = cls.available_metrics
#         # If metrics is a list, filter the default_metrics dictionary
#         elif isinstance(metrics, list):
#             metrics = {metric.lower(): cls.available_metrics[metric.lower()] for metric in metrics}
#
#         metrics_data = {}
#         for graph_name, graph in graph_dict.items():
#             metrics_data[graph_name] = {}
#             for metric_name, metric_function in metrics.items():
#                 metric_value = metric_function(graph)
#                 metrics_data[graph_name][metric_name] = metric_value
#
#         return metrics_data
#
#     @classmethod
#     def print_metrics(cls, metrics_dict, metric_names=None):
#         if metric_names is None:
#             metric_names = []
#
#         for graph_name, metrics in metrics_dict.items():
#             print(f"Metrics for {graph_name}:")
#             for metric_name, metric_value in metrics.items():
#                 if metric_name in metric_names and callable(metric_value):
#                     try:
#                         metric_value = metric_value()
#                     except TypeError:
#                         if metric_name == "fg_modularity":
#                             membership = metric_value.community_multilevel().membership
#                             metric_value = metric_value.modularity(membership)
#                         else:
#                             raise ValueError(f"Unsupported metric '{metric_name}'")
#                 print(f"  {metric_name}: {metric_value}")
#             print()



# class InteractivePlotNetworkxGraphs:
#
#     def __init__(self, networkx_graph_dict):
#         self.networkx_graph_dict = networkx_graph_dict
#
#     def plot(self):
#         for graph_key, nx_graph in self.networkx_graph_dict.items():
#             # Use a spring layout to position the nodes
#             pos = nx.spring_layout(nx_graph, seed=42)
#
#             # Create node trace
#             node_trace = go.Scatter(mode='markers+text', textposition='bottom center',
#                                     hoverinfo='text', marker=dict(showscale=False, size=20, colorscale='Viridis', reversescale=True, colorbar=dict(thickness=15, title='Node Connections', xanchor='left', titleside='right')))
#
#             # Create edge trace
#             edge_trace = go.Scatter(line=dict(width=0.5, color='#888'), hoverinfo='none', mode='lines')
#
#             node_x, node_y, node_text, node_colors = [], [], [], []
#             edge_x, edge_y = [], []
#
#             # Add nodes and edges to the traces
#             for node, adjacencies in enumerate(nx_graph.adjacency()):
#                 x, y = pos[node]
#                 node_x.append(x)
#                 node_y.append(y)
#                 node_colors.append(len(adjacencies[1]))
#                 node_info = f"{node} - # of connections: {len(adjacencies[1])}"
#                 node_text.append(node_info)
#
#                 for neighbor in adjacencies[1].keys():
#                     x0, y0 = pos[node]
#                     x1, y1 = pos[neighbor]
#                     edge_x.extend([x0, x1, None])
#                     edge_y.extend([y0, y1, None])
#
#             node_trace.update(x=node_x, y=node_y, text=node_text, marker=dict(color=node_colors))
#             edge_trace.update(x=edge_x, y=edge_y)
#
#             # Customize the plot appearance
#             fig = go.Figure(data=[edge_trace, node_trace],
#                             layout=go.Layout(title=graph_key, showlegend=False, hovermode='closest',
#                                              margin=dict(b=20, l=5, r=5, t=40),
#                                              xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
#                                              yaxis=dict(showgrid=False, zeroline=False, showticklabels=False)))
#             fig.show()


# interactive_plotter = InteractivePlotNetworkxGraphs(ConvertIgraphToNetworkx(chr18_inc_norm_graphs_50kb()).convert())
# interactive_plotter.plot()