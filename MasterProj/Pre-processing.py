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

@dataclass
class NodeList:
    nodes: list

    def __iter__(self):
        return iter(astuple(self))

    def print_class(self):
        print(f"Nodes: {self.nodes}")


    def as_list(self):
        return self.nodes

    def group_cellines(self):
        pass

    def group_chromosomes(self):
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


@dataclass
class Celline:
    strain: str
    nodes: NodeList

    def __iter__(self):
        return self

    def print_self(self):
        print(f"Strain: {self.strain} with nodelist: {self.nodes}")

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
    return NodeList(nodes)

# nodes = process_file_to_node("/Users/GBS/Master/HiC-Data/Hi-C_data_fra_Jonas/4linescopy/IMR90_50kb.domain.RAW.no_cen.NCHG_fdr.o_by_e5_to_plot.gtrack")


# Secelting directory containing files to process and instantiating Celline class
def process_directory_to_celline(directory):
    paths = []
    for file in os.listdir(directory):
        paths.append(os.path.join(directory, file))
    return process_files_to_celline(paths)

def process_files_to_celline(files):
    cellines = []
    for arg in files:
        with open(arg):
            strain = arg.split("/")[7].split("_")[0]
            if not strain.startswith("."):
                cellines.append(Celline(strain, process_file_to_node(arg)))
    return cellines

cellines = process_directory_to_celline("/Users/GBS/Master/HiC-Data/Hi-C_data_fra_Jonas/4cell_lines_Hi-C")

# # Calling repr method
# for celline in cellines:
#     celline.print_self()
