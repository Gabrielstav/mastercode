from dataclasses import dataclass
from Node_Standardization import SetDirectories
import os as os

# TODO: Refactor to read in edge lists from Pipeline
# TODO: Finish chromosome and resolution filter methods
# TODO: Make "as list" and "as dict" methods
# TODO: Calculate overlap between nodes across cell lines, might not have any overlap in nodes

# Set directories containing edge lists (theycan be names anything, but if they come from the Pipeline class they follow the standard naming convention.
# So check if they follow the standard naming convention, if they do,
# Standard naming: experiment_resolution_edgelist.txt, eg: IMR90_10000_edgelist.txt or chr18_lowres_500000_edgelist.txt
# The root folder containing each set of raw data is the experiment name, eg in out case chr18_lowres (can be for celline as well).

# Where experiment is the Job name in HiC-Pro, and resolution is the resolution of the Hi-C data in bp.
# But IDK what I'll name things yet when I run things on HPC yet, so just leave as is now I think, and then refactor later.

# After setting input directories, we can read in the files, and process them to nodes.
# The Node class has the following attributes: ID, edges, chromosome, start, end. Think start and end is nice for calculating overla with other nodes from other experiments later.
# Or maybe redundant since Bedtools might be possible to use? IDK yet.

# The next level of abstraction should be experiment dataclass (Cell line, chromosome_lowres etc)
# containing list of all nodes for that experiment, as well as metadata like resolution, experiment name (cell line or chromosome_hires etc)
# Methods can be used to filter on resolution, chromosome, cell line etc.

# The next level of abstraction should be a class containing all experiments, and can be used to calculate overlap between nodes across experiments.
# We also need to call filtering methods on the experiments, to filter on resolution, chromosome etc.
# It makes no sense to compare data across resolution since it compares different genomic regions defined by different contraints, so we need to filter on resolution.
# Method can be used to calculate overlap between nodes across experiments, and the filtering things defined in the experiment/Node class.

# TODO: The goal of this (and the only thing I NEED to do) is to be able to calculate overlap between nodes across experiments. So start with that.

class SetInputDir:
    """
    Sets input directory for processing, by default same as output directory used in Pipeline.
    Overwrite default input directory by defining user defined input directory.
    """

    # Why does this run Pipeline? Shouldn't it just set the input directory?

    user_defined_input_dirs = os.path.abspath("")

    @classmethod
    def default_input_dir(cls):
        default_input_dir = os.path.abspath(SetDirectories.get_output_dir())
        return default_input_dir

    @classmethod
    def defined_input_dir(cls, user_defined_input_dir):
        defined_input_dir = os.path.abspath(user_defined_input_dir)
        return defined_input_dir

    @classmethod
    def get_input_dir(cls):
        if SetInputDir.default_input_dir() is None and SetInputDir.defined_input_dir() is None:
            return print("No input directory defined. Please define an input directory in SetInputDir class.")
        if SetInputDir.defined_input_dir() is not None and SetInputDir.default_input_dir() is not None:
            return print("Both default and user defined input directories are set. Please remove one of them.")
        if SetInputDir.defined_input_dir() is None:
            return [os.path.abspath(dirs) for dirs in SetInputDir.user_defined_input_dirs.split(";")]
        else:
            return [os.path.abspath(SetInputDir.default_input_dir())]


@dataclass(frozen=True, eq=True)
class Node_test:
    id: str
    edges = list[str]
    chr: str
    start: int
    end: int

    @classmethod
    def process_edgelist_to_node(cls, *input_directories):
        """
        This method processes the files to nodes, and stores them in a list.
        Input files are edge lists processed by the Pipeline class.
        :param args: Directory containing the edge lists
        :return: Instances of Node class
        """

        nodes = []
        for input_dir in input_directories:
            for file in os.listdir(input_dir):
                if file.endswith(".txt"):
                    with open(os.path.join(input_dir, file)) as f:
                        for line in f:
                            columns = line.strip().split("\t")
                            # node = Node_test(columns[])
                            node = Node_test(columns[3], columns[6].split(";"), columns[0], int(columns[4]), int(columns[5]))
                            nodes.append(node)
        return nodes

    """
    chr18:200000-250000  chr18-2650000:2700000
    chr18:200000-250000  chr18-2700000:2750000
    chr18:250000-300000  chr18-1000000:1050000
    """

    """ 
    edgelists/chr18_lowres_50000_edgelist.txt
    """


@dataclass(frozen=True, eq=True)
class Weighted_Node:
    id: str
    edges: list[str]
    chr: str
    start: int
    end: int
    weight: float

    @classmethod
    def process_weighted_edgelist_to_node(cls, *args):
        """
        This method processes the files to nodes, and stores them in a list.
        Input files are weighted edge lists processed by the Pipeline class.
        :param args:
        :return:
        """

        # Same as Node_test class but with weight attribute added

@dataclass(frozen=True, eq=True)
class Experiment:
    identifier: str  # cell line, experiment name, etc
    resolution: int






@dataclass(frozen=True, eq=True)
class Node:
    id: str
    edges: list[str]
    chr: str

    def as_list(self):
        for edge in self.edges:
            return edge.as_list()

    def all_nodes(self):
        return self.id + self.edges + self.chr

    def is_isolated(self):
        return self.edges == "."

    def is_connected(self):
        return self.edges != "."

    def get_chromosome(self):
        return self.chr

    @classmethod
    def process_file_to_node(cls, *args):
        nodelist = []
        for arg in args:
            with open(arg, encoding="latin-1") as file:
                file_content: list[str] = file.readlines()
        for index, line in enumerate(file_content):
            if line.startswith("chr"):
                columns = line.strip("\n").split("\t")
                if len(columns) != 7:
                    print(f"Bad line format: {columns} with index: {index}")
                else:
                    if columns[6] == ".":
                        nodelist.append(Node(columns[3], columns[6], columns[0]))
                    else:
                        nodelist.append(Node(columns[3], columns[6].split(";"), columns[0]))
        return nodelist

    @classmethod
    def print_all_nodes(cls, nodelist):
        newline = "\n"
        for nodes in nodelist:
            print(f"Node: {nodes.id} with edges: {nodes.edges} on chromosome: {nodes.chr} {newline}")

    # @classmethod
    # def nodelist_oit(cls, nodelist): skriv nodelist_out funksjon for empty nodes

    @classmethod
    def nodelist_out(cls, nodelist):
        for nodes in nodelist:
            if nodes.edges != ".":
                print(f"Node: {nodes.id} with edges: {nodes.edges} on chromosome: {nodes.chr} is connected")  # need to write this as list

    # create function that outputs list of all nodes on a chromosome
    # @classmethod
    # def chromosome_nodes(cls, *args):
    #     for arg in args:
    #         get_chroms = list(filter(lambda x: arg == str(x.chr), Node.as_list()))
    #         return get_chroms



# @dataclass(frozen=True, eq=True)
# class Celline:
#     strain: str
#     nodes: list[Node]  # Not just a node, slå sammen nodelist og celline?
#
#     def nodes(self):
#         return self.nodes
#
#     def number_of_nodes_string(self):
#         return f"{self.nodes.__sizeof__()} nodes"
#
#     def celline_string(self):
#         return f"celline {self.strain} with {self.nodes}"
#
#     def print_strain(self):
#         return f"Cellines: {self.strain}"
#
#     def as_list(self):
#         return self.nodes
#
#     def only_iso_nodes(self):
#         iso_nodes = list(filter(lambda n: n.is_isolated(), self.nodes))
#         return Celline(self.strain, iso_nodes)
#
#     def only_con_nodes(self):
#         con_nodes = list(filter(lambda n: n.is_connected(), self.nodes))
#         return Celline(self.strain, con_nodes)
#
#     def group_by_chromosome(self, *args):
#         for arg in args:
#             selected_chromosomes = list(filter(lambda n: arg == n.get_chromosome(), self.nodes))
#             return Celline(self.strain, selected_chromosomes)
#
#
# @dataclass(frozen=True, eq=True)
# class Cellines:
#     cellines: list[Celline]
#
#     # defualt dir to make calling methods faster
#     @classmethod
#     def from_default(cls):
#         return cls.from_dir("/Users/GBS/Master/HiC-Data/Hi-C_data_fra_Jonas/4linescopy")
#
#     # instantiating Celline class
#     # and pre-processing gtrack files in directory
#     @classmethod
#     def from_dir(cls, directory):
#         paths = []
#         for file in os.listdir(directory):
#             paths.append(os.path.join(directory, file))
#         return Cellines.from_files(paths)
#
#     @classmethod
#     def from_files(cls, files):
#         celline_list = []
#         for arg in files:
#             with open(arg):
#                 strain = arg.split("/")[7].split("_")[0]
#                 if not strain.startswith("."):
#                     celline_list.append(Celline(str(strain), Node.process_file_to_node(arg)))
#         return Cellines(celline_list)
#
#     @classmethod
#     def by_strain_dictionary(cls, celline_list):
#         cellines_by_strain = {}
#         for cline in celline_list:
#             strain_dict = cellines_by_strain[cline.strain, cline]
#         return strain_dict
#
#     def print(self):
#         newline = "\n"
#         for cline in self.cellines:
#             print(f"Strain: {cline.strain} with Nodelist: {cline.nodes} {newline}")
#
#     def with_strain(self, *args):
#         cellines_with_strain = []
#         for strain in args:
#             for celline in self.cellines:
#                 if strain == celline.strain:
#                     cellines_with_strain.append(celline)
#         return Cellines(cellines_with_strain)
#
#     def only_iso(self):
#         iso_list = list(map(lambda c: c.only_iso_nodes(), self.cellines))
#         return Cellines(iso_list)
#
#     def only_con(self):
#         con_list = list(map(lambda c: c.only_con_nodes(), self.cellines))
#         return Cellines(con_list)
#
#     def with_chromosome(self, *args):
#         for arg in args:
#             chromosomes_gotten = list(map(lambda c: c.group_by_chromosome(arg), self.cellines))
#         return Cellines(chromosomes_gotten)
#
#
# # list of all cellines
# all_cellines = Cellines.from_default()
#
# # examples of use
# K562_iso = Cellines.from_default().with_strain("K562").only_iso()
# HMEC_con = Cellines.from_default().with_strain("HMEC").only_con()
# HMEC_iso = Cellines.from_default().with_strain("HMEC").only_iso()  # .chromosome("1-10").as_dict() OR as.list()
#
# print(HMEC_con)
#

# Need to add functionality for filtering on chromosomes
# and maybe write function that creates objects of the cell lines automatically
# and maybe write function that can "export" the cell line objects in a correct format, check the iGraph API for instance

# Write method for chromosome filtering and transposing input
# HUVEC_chr1_to_10 = Cellines.from_default().with_strain("HUVEC").only_con().with_chromosome("1-10")

