from dataclasses import dataclass
import os


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
                print(f"Node: {nodes.id} with edges: {nodes.edges} on chromosome: {nodes.chr} is connected") # need to write this as list


    # create function that outputs list of all nodes on a chromosome
    @classmethod
    def chromosome_nodes(cls, *args):
        for arg in args:
            get_chroms = list(filter(lambda x: arg == str(x.chr), nodes.as_list()))
            return get_chroms





@dataclass(frozen=True, eq=True)
class Celline:
    strain: str
    nodes: Node.nodelist # Not just a node, slå sammen nodelist og celline?

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

    #


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

    def print(self):
        newline = "\n"
        for cline in self.cellines:
            print(f"Strain: {cline.strain} with Nodelist: {cline.nodes} {newline}")

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


# list of all cellines
all_cellines = Cellines.from_default()

# examples of use
K562_iso = Cellines.from_default().with_strain("K562").only_iso()
HMEC_con = Cellines.from_default().with_strain("HMEC").only_con()
HMEC_iso = Cellines.from_default().with_strain("HMEC").only_iso()

# Need to add functionality for filtering on chromosomes
# and maybe write function that creates objects of the cell lines automatically
# and maybe write function that can "export" the cell line objects in a correct format, check the iGraph API for instance

# Write method for chromosome filtering and transposing input
# HUVEC_chr1_to_10 = Cellines.from_default().with_strain("HUVEC").only_con().with_chromosome("1-10")

# what do i google to find solution to sdk not configured error?
# how do i get the sdk to be configured?
# how do i get the sdk to be configured on my mac?
# and the solution is?
# how do i get the sdk to be configured on my mac?
# and the solution is?
