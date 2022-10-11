from dataclasses import dataclass, astuple, fields

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


for celline in cellines:
    celline.print_self()
