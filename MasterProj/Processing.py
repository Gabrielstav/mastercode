from dataclasses import dataclass, astuple
import os

""" Creating classes """


@dataclass
class Node:
    id: str
    edges: list[str]
    chr: str

    def __iter__(self):
        return iter(astuple(self))

    def is_isolated(self):
        return self.edges == "."

    def is_connecte(self):
        return self.edges != "."


# @dataclass
# class NodeList:
#     nodes: list
#
#     def nodes(self):
#         return self.nodes
#
#     def as_list(self):
#         return self.nodes
#
#     def only_iso_nodes(self):
#         return NodeList(list(filter(lambda node: node.is_empty(), self.nodes)))


# (frozen=True, eq=True)
@dataclass(frozen=True, eq=True)
class Celline:
    strain: str
    nodes: list

    def nodes(self):
        return self.nodes

    def as_list(self):
        return self.nodes

    def only_iso_nodes(self):
        return Celline(list(filter(lambda node: node.is_isolated(), self.nodes)))

    def only_con_nodes(self):
        return Celline(list(filter(lambda node: node.is_connected(), self.nodes)))


    def to_string(self):
        return f"{self.nodes.__sizeof__()} nodes"

    @classmethod
    def print_cellines(cls, celline_list):
        newline = "\n"
        for celline in celline_list.values():
            print(f"Strain: {celline.strain} with nodelist: {celline.nodes} {newline}")

    @classmethod
    def print_cellines(cls, celline_list):
        for celline in celline_list.values():
            print(celline.nodes.only_iso_nodes().to_string())

    def to_string(self):
        return f"celline '" + self.strain + " with " + nodes.to_string()




def process_file_to_node(*args):
    nodes = []
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
                    nodes.append(Node(columns[3], columns[6], columns[0]))
                else:
                    nodes.append(Node(columns[3], columns[6].split(";"), columns[0]))
    return Celline(nodes)


nodes = process_file_to_node("/Users/GBS/Master/HiC-Data/Hi-C_data_fra_Jonas/4linescopy/IMR90_50kb.domain.RAW.no_cen.NCHG_fdr.o_by_e5_to_plot.gtrack")


# Seelcting directory containing files to process and instantiating Celline class
def process_directory_to_celline(directory):
    paths = []
    for file in os.listdir(directory):
        paths.append(os.path.join(directory, file))
    return process_files_to_celline(paths)


def process_files_to_celline(files):
    cellines = {}
    for arg in files:
        with open(arg):
            strain = arg.split("/")[7].split("_")[0]
            if not strain.startswith("."):
                celline = Celline(str(strain), process_file_to_node(arg))
                cellines[celline.strain] = celline
    return cellines


cellines = process_directory_to_celline("/Users/GBS/Master/HiC-Data/Hi-C_data_fra_Jonas/folder_test")

Celline.print_cellines(cellines)








"""    def group_chromosomes(self):
        pass

    def group_connected_nodes(self):
        pass

    def group_isolated_nodes(self):
        pass

    def isolated_degree(self):
        pass

    def node_overlap(self):
        pass

    def standardize_nodes(self):
        pass
"""





""" Writing pre-processed output to new files """

# Use to OS module to check what directory we are in, read the files put in a specific directory (to be pre-processed),
# checks their format (only .gtrack compatibility as of now) and then run the pre-processing parts of the code,
# and outputs pre-processed files in the edge list format to a specific directory with specific file names,
# maybe just the original file names with "processed" or "edge-list" behind the name.


# Pseudocode-function to automatically process files and write to a new dir if we provide the dir:
# count = 0
# directory = "/Users/GBS/Master/HiC-Data/Processed_Data/HiC_from_Jonas/FullGenome"
# for filename in os.listdir(directory):
#     if filename.endswith(".txt"):
#         count += 1
#         with open(directory + filename, "r") as read_file:
#            return_of_your_function = do_something_with_data() # Putt pre-processinga over inn en funksjon som calles her
#         with open(directory + count + filename, "w") as write_file:
#             write_file.write(return_of_your_function)  !!! igjen bare caller funksjonen her !!!

# Selecting full genome and writing edge list to new file
# with open("/Users/GBS/Master/HiC-Data/Processed_Data/HiC_from_Jonas/FullGenome/TESTING", "w") as fp:
#     for node in nodes:
#         for edge in node.edges:
#             fp.write((node.id + " " + edge + "\n"))


# Selecting specific chromsome and writing edge list to new file
# with open("/Users/GBS/Master/HiC-Data/Processed_Data/Processed_J_data/HiC_from_Jonas/Chr18/IMR90_processed_chr18.txt", "w") as fp:
#     for l in nodes:
#         for edge in l.edges:
#              if l.chr == "chr18":  # Make a function that we can call to select id + edges from chromosome(s)?
#                  fp.write((l.id + " " + edge + "\n"))

# This outputs a file with the correct format but with header: id, edges (which should be removed).
# Use with open(file) that points to file opened earlier, should be able to find the path to the file opened when
# pre-processing the data, so that we only need to point to the full path once and the rest of the code knows what
# file(s) we are working with.
