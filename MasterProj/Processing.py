from dataclasses import dataclass, astuple
import os
import matplotlib as mpl
import matplotlib.pyplot as plt


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

    # Instantiating Node class
    # and pre-processing the gtrack data

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


@dataclass(frozen=True, eq=True)
class Celline:
    strain: str
    nodes: list[Node]

    def nodes(self):
        return self.nodes

    def number_of_nodes_string(self):
        return f"{self.nodes.__sizeof__()} nodes"

    def celline_string(self):
        return f"celline {self.strain} with {self.nodes}"

    def print_strain(self):
        return f"Cellines: {self.strain}"

    def as_list(self):
        return self.nodes

    def only_iso_nodes(self):
        iso_nodes = list(filter(lambda n: n.is_isolated(), self.nodes))
        return Celline(self.strain, iso_nodes)

    def only_con_nodes(self):
        con_nodes = list(filter(lambda n: n.is_connected(), self.nodes))
        return Celline(self.strain, con_nodes)

    def group_by_chromosome(self, *args):
        for arg in args:
            selected_chromosomes = list(filter(lambda n: arg == n.get_chromosome(), self.nodes))
            return Celline(self.strain, selected_chromosomes)




@dataclass(frozen=True, eq=True)
class Cellines:
    cellines: list[Celline]

    # defualt dir to make calling methods faster
    @classmethod
    def from_default(cls):
        return cls.from_dir("/Users/GBS/Master/HiC-Data/Hi-C_data_fra_Jonas/4linescopy")

    # instantiating Celline class
    # and pre-processing gtrack files in directory
    @classmethod
    def from_dir(cls, directory):
        paths = []
        for file in os.listdir(directory):
            paths.append(os.path.join(directory, file))
        return Cellines.from_files(paths)

    @classmethod
    def from_files(cls, files):
        celline_list = []
        for arg in files:
            with open(arg):
                strain = arg.split("/")[7].split("_")[0]
                if not strain.startswith("."):
                    celline_list.append(Celline(str(strain), Node.process_file_to_node(arg)))
        return Cellines(celline_list)

    @classmethod
    def by_strain_dictionary(cls, celline_list):
        cellines_by_strain = {}
        for cline in celline_list:
            strain_dict = cellines_by_strain[cline.strain, cline]
        return strain_dict

    # @classmethod
    # def print_celline(cls, celline_list):
    #     newline = "\n"
    #     for cline in celline_list:
    #         print(f"Strain: {cline.strain} with Nodelist: {cline.nodes} {newline}")

    def print(self):
        newline = "\n"
        for cline in self.cellines:
            print(f"Strain: {cline.strain} with Nodelist: {cline.nodes} {newline}")

    # @classmethod
    # def with_strain(cls, *args, celline_list):
    #     cellines_with_strain = []
    #     for arg in args:
    #         for cline in celline_list:
    #             if arg == cline.strain:
    #                 cellines_with_strain.append(cline)
    #     return cellines_with_strain

    def with_strain(self, *args):
        cellines_with_strain = []
        for strain in args:
            for celline in self.cellines:
                if strain == celline.strain:
                    cellines_with_strain.append(celline)
        return Cellines(cellines_with_strain)

    def only_iso(self):
        iso_list = list(map(lambda c: c.only_iso_nodes(), self.cellines))
        return Cellines(iso_list)

    def only_con(self):
        con_list = list(map(lambda c: c.only_con_nodes(), self.cellines))
        return Cellines(con_list)

    def with_chromosome(self, *args):
        for arg in args:
            chromosomes_gotten = list(map(lambda c: c.group_by_chromosome(arg), self.cellines))
        return Cellines(chromosomes_gotten)

    # @classmethod
    # def with_chromosome(cls, *args, celline_list):
    #     for arg in args:
    #         chromosomes_gotten = list(map(lambda c: c.group_by_chromosome(arg), celline_list))
    #         return chromosomes_gotten


all_cellines = Cellines.from_default()


K562_iso = Cellines.from_dir("/Users/GBS/Master/HiC-Data/Hi-C_data_fra_Jonas/4linescopy").with_strain("K562").only_iso()
HMEC_con = Cellines.from_dir("/Users/GBS/Master/HiC-Data/Hi-C_data_fra_Jonas/4linescopy").with_strain("HMEC").only_con()

# print(K562_iso)
#print(HMEC_con)



# Returner både connected OG isolated noder for K562 cellelinjen

