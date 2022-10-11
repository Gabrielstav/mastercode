### TIPS

# zip function loops over two lists at once, so maybe save nodes in one and edges in another
# and save edges so the index corresponds to the index of nodes but each nodelist is a sublist

# Unpacking with a, b, *c, d = (1, 2, 3, 4, 5, 6) = (a=1, b=2, c= 3,4,5, d=6)

# setaddr and getaddr sets attributes, useful for looping over things and adding them to attributes


### Old non-args-version of function:

# chromlist = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
#              "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX"]

# def get_nodes_from_chromosome(arg):
#     args = list(arg)
#
#     """ creates lists containing empty and connected nodes from selected chromosomes: chromlist[0:22]
#         and calculates some statistics on empty vs connected nodes for the selected chromosome(s)"""
#
#     # creating empty and connected node lists
#     for x in nodes_and_edges:
#         if x.edges == "." and x.chr in arg:
#             chromlist_empty_nodes.append(Nodes(x.id, x.edges, x.chr))
#             # print(chromlist_empty_nodes)
#         elif x.edges != "." and x.chr in arg:
#             chromlist_connected_nodes.append(Nodes(x.id, x.edges, x.chr))
#             # print(chromlist_connected_nodes)
#
#     # total nodes
#     total_nodes = len(chromlist_empty_nodes) + len(chromlist_connected_nodes)
#     if len(arg) <= 1:
#         print(f"Total nodes in selected chromosome {arg} is {total_nodes}")
#     if len(arg) >= 2:
#         print(f"Total nodes in selected chromosomes {arg} is {total_nodes}")
#
#     # sub-chromosomal nodes
#
#     for x in chromlist_empty_nodes and chromlist_connected_nodes:
#         if arg[x] == x.chr in chromlist_empty_nodes and chromlist_connected_nodes:
#             print(f"Chromosome {arg[x]} has {len(chromlist_empty_nodes)+len(chromlist_connected_nodes)}")

### Function to read in mulitple files at once and automate the pre-processing process:


# def empty_node_stat(chromosomes):
#
#     """
#     :param chromosomes: what chromosomes we want to group by
#     :param dataset: our dataset that is pre-processed and contans empty nodes only
#     :return: this should return empty nodes with edges grouped by chromosome
#     """
#
#     chrlist = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11" "chr12",
#                "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX"]
#     dataset_empty = empty_nodes
#     dataset_connected = nodes_and_edges
#     chromosome_prefix = "chr"

    # Dette "henter" ut empty edges for spesifikke kromosom, men koden må iterere over listen og gruppere
    # alle kromosomene. Trenger ikke å printe det, men må gruppere på kromosom for å gjøre stats.
    # for x in dataset_empty:
    #     if x.chr == chromosome_prefix + str(chromosome_numbers):
    #         print(x.id + " " + x.edges)


# Filtering and string comprehension + functions

# Kanskje adde en method til EmptyNodes som outputter det vi skal ha? overrider __eq__?
# Den sier at "chr18" ikke er i listen, men det er det jo.
# class EmptyNodes:
#     if __name__ == "__main__":
#         idx = empty_nodes.index('chr18')
#         print(f"Empty nodes in chr18 {empty_nodes[idx].chr_empty}")
#
# print(empty_nodes)

# Kanskje lag en funksjon som gjør det jeg vil (ikke fått til enda) og calle den i filter:
# list(x for x in empty_nodes if chr_empty "our function"(x))

# Hva med fields da? type object. 'EmptyNodes' has no attribute 'id_empty'. Joda, det har den.
# for field in fields(EmptyNodes):
#     print(field.name, getattr(EmptyNodes, field.name))


# How about this? Nope, EmptyNodes is not iterable.
# next(chr_empty for node in EmptyNodes if node.chr_empty == "chr18")


# How to use filter function:
# scores = [1, 2, 3, 4, 5, 6]
# filtered = filter(lambda s: s >= 4, scores)
# print(list(filtered))

# Testing filter on empty_nodes but with length instead:
# But this doesn't work because EmptyNodes (class) has no length, is not iterable and not subscriptable.
# Any iterable can be the second arg of filter, and the function checks for boolean statements
# but we want multiple "if" checks depending on how many chromosomes we are filtering for.

# for x in file_content:
#     for y in empty_nodes:
#         filtered_list2 = [x for x in empty_nodes if x == "chr18"]
#         print(filtered_list2)




# for x in file_content:
#     for y in empty_nodes:
#         filteredlist = list(filter(lambda x: x == "chr18:", empty_nodes))
#         print(filteredlist)

# But this doesn't work. Empty_Nodes class object isn't iterable and when I use the list I get an empty list returned:
# filtered_chr18 = filter(lambda s: s == "chr18", empty_nodes)
# print(list(filtered_chr18))

# Maybe we need to write EmptyNodes() or something, so that the filter function doesn't return
# any string being equal to "chrx" but the whole object where the "chrx" string is found?


# This doesn't work either because EmptyNodes (class) is not iterable.
# This would give the wrong output even if it worked, because we just look for any string containing "chr18" and "chr2".
# process_these_chromosomes = ["chr18", "chr2"]
# def Filter(process_these_chromosomes, empty_nodes):
#     return [n for n in empty_nodes if
#             any(m in n for m in process_these_chromosomes)]
# print("filter data:", Filter(process_these_chromosomes, empty_nodes))

# Ok trying this method then, but I don't understand how it'll help, it just returns "chr18":
# sublist = ["chr18", "chr2"]
# def Filter2(list):
#     return [val for val in empty_nodes if re.search(r"^chr18", val)]
# print(Filter2(sublist))

# So I can't use this method, and instead need to write a function with different amounts of arguments?
# Either I can write "chrx" and get stats for that chromosome, or write "chrx-y" and get info for those chromosomes.

# chrom: list[str] = f"[chr2, chr18]"

# print(chrom[0])
# counter = 0
# for index, counter in enumerate(chrom):
#     print(chrom[counter], index)
#     counter += 1

# def get_empty_nodes(args):
#     chrlist = []
#     while counter >= len(chrom):
#         for z in empty_nodes:
#             for edge inx z.edges_empty:
#                 if z.chr_empty == chromosomes:
#                     chrlist.append(z.id_empty + " " + edge)
#         counter += 1
#         if counter > len(list(chrom)):
#             break
#     return print(chrlist)
#
# get_empty_nodes("chr2", "chr18")


# We also want the function to output metrics like this:

#     print("Total nodes:", len(empty_nodes) + len(nodes_and_edges))
#     print("Empty nodes:", len(empty_nodes))
#     print("Connected nodes:", len(nodes_and_edges))
#     print("Empty percentage:", len(empty_nodes)/(len(nodes_and_edges)+len(empty_nodes)), "%")
#     print("Connected percentage:", len(nodes_and_edges) / (len(nodes_and_edges) + len(empty_nodes)), "%")

# Total nodes in "cell line/chr": 5821
# Empty nodes: 3332
# Connected nodes: 2489
# Proportion of empty nodes: 0.5724102387905858
# Proportion of connected nodes: 0.42758976120941417 (round and give as %?)
# Nodes with x edges (and which nodes with names).

# So maybe we can use kwargs to create a dict where each chromosome is the key
# adn the values are connected/empty ratio or something? Each list can be a value for a key.
# So maybe we can do this for the function above, splitting list for each chromosome passed as an argument
# and creating a list for each chromosome. Then the chromosome list can be used as arg for the new function which
# creates a dict where each key is a chromosome and each vlaue is a list of the connected and empty nodes in tht
# chromosome (maybe one dict for empty and one for connected). Then after this we can use this dict (or just a function
# that uses the lists) to compute node statistics.
# Maybe just create dict with keys for chrom and empty values as lists, and append to them with iter in func.
# Extend dict by .extend or += / append / or extend: where dict["keyname"].extend(["values"]).
# add method is the setdefaultmethod which uses extend method is another way.
# Append creates sublists withtin the list of keys tho, so kinda bad maybe, but just idex twice so it works.


# New function that outputs some stats on empty nodes.
# It can loop over the input/output from the get empty nodes function,
# only selecting wanted chromosomes and storing the infomration as
# data that can be plotted. Maybe a dict or a list, we want to make a graph of
# empty nodes ratio for each chromosoem and for each cell line.


# Our functions should output something like this:
# Total nodes in "cell line/chr": 5821
# Empty nodes: 3332
# Connected nodes: 2489
# Proportion of empty nodes: 0.5724102387905858
# Proportion of connected nodes: 0.42758976120941417 (round and give as %?)
# Nodes with x edges (and which nodes with names).
# A list of the chromosomes we selected, so that we can write this to a file or something.
# Because we want to examine the statistics of the connectedness intra- and interchromosomally.

#     print("Total nodes:", len(empty_nodes) + len(nodes))
#     print("Empty nodes:", len(empty_nodes))
#     print("Connected nodes:", len(nodes))
#     print("Empty percentage:", len(empty_nodes)/(len(nodes)+len(empty_nodes)), "%")
#     print("Connected percentage:", len(nodes) / (len(nodes) + len(empty_nodes)), "%")

# We want to save each chr as a key and each ratio count of empty/connected nodes as values.
# So we get two dictionaries, one for the empty and one for the connected nodes.
# Each chromosome key in both dictionaries values (node ratio) should sum to 1 (they're complimentary).
# So maybe we only need one dict? Idk.


# Open all files in a directory and and process them? This needs to open files one by one in a directory,
# and process them by calling the functions above.


# TESTING FUNCTION
# chrlist = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11" "chr12",
#            "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX"]
# dataset_empty = empty_nodes
# dataset_connected = nodes_and_edges
# chromosome_prefix = "chr"

# Jeg vil ha all data i en dataklasse: Nodes eller noe.
# Så m jeg har ulike constructors for denne dataklassen, så jeg kan lese inn data
# og lage objekter med ulike argumenter. Jeg vil lage en attribute som er empty nodes,
# og objektene fra den skal bare ta argumentene id, edges (tomme) og chr, samme som
# objektene som er connected nodes.

# Printing all nodes with edges
# for x in nodes_and_edges:
#     for edge in x.edges:
#         print(x.id + " " + edge)

# Selecting nodes and edges for specific chromosomes (should be a better way to do this?)
# for l in nodes_and_edges:
#     for edge in l.edges:
#          if l.chr == "chr18":  # Make a function that we can call to select id + edges from chromosome(s)?
#              print(l.id + " " + edge)

# counter = 0
# for x in dataset_empty:
#     if chrlist[counter] == x.chr:
#         print(x.id + " " + x.edges)
#     counter += 1
#     if counter > len(chrlist):
#         break

# counter = 0
# for x in dataset_empty:
#     while counter >= len(chrlist):
#         if x.chr == chrlist[counter]:
#             print(x.id + " " + x.edges)
#         counter = counter+1


# counter = 0
# while counter <= len(chrlist):
#     for x in empty_nodes:
#         for edge in x.edges:
#             if x.chr in empty_nodes == chrlist[counter]:
#                 print(x.id + " " + x.edges)
#     counter += 0
#     if counter >= len(chrlist):
#         break

# for x in empty_nodes:
#     for edge in x.edges:
#         #if x.chr == "chr18":
#         print(x.id + " " + edge)


# path = ""
# process_input = os.listdir(path)
# def process(process_input):
#     for files in process_input:
#         if files in process_input == "gtrack":
            # do the pre-processing
            # for each file separately
            # and write separate outputs for each file

### Function to write output file

# Writing output to file (is there a better way to do this, where you can call a function instead?)
# Here we need to first create a file in the dir we want, and then copy the full path and put it here
# I want to be able to call a function (write_output) that just automatically creates a textfile
# with the correct edge list format, where we just specify the directory in which it should be placed.

# Might need to use the OS module to avoid cross-platform issues? If the dir is not there it returns an exception.
# The current working dir resolves the relative paths, meaning that we need to find the absolute path to the
# folder we want:
# import os
# dirname = os.path.dirname(__file__)
# filename = os.path.join(dirname, "relative path to our output folder")

#def write_to_file(path, input):


    # :param path: to directory we want to write to
    # :param input: pre-processed data from Node dataclass
    # :return: output is writing edge list to new file in directory
    #



### On classes:


"""
class Node:

    def __init__(self, id: str, edges: str):  # Constructor
        self.id = str(id)  # Attribute
        self.edges = list(edges)  # Attributes are passed to the initialize function

        # Methods go here (like functions), which decide what operations we can do with objects of class 'Node'

        def __repr__(self):
            return print({self.id, self.edges}) #Manually implement representation of attributes (use dataclass)

        def number_of_nodes(cls):
            if range(len(nodes)) > 0:
                print(len(nodes))
            return cls.number_of_nodes

        # def add_nodes_to_list:
        # return nodes.append()
"""









"""
file_content = file.readlines()
for lines in file_content:
    columns.append(str((lines.split("\t")))
for index, lines in enumerate(file_content):
    if len(columns) != 7:
        print("Bad line format")
    else:
        #print("length of columns:", len(columns)) #test
        #print("Column index 3:", columns[3]) #test
        edges = columns[6].split(";")
        nodes_and_edges = Node(columns[3], edges)


for line in file_content:
    columns = line.split("\t")
    if len(columns) != 7:
        print("Bad line format")
    edges = columns[6].split(";")
    nodes.append(Node(columns[3], edges))"""

# So in the class Node we have ID (node) and its edges as attributes of the class. This means for every node object
# we need to pass an ID (str) and a node (list of strings) So when we create new instances of the Node class,
# all we need to do is read each line in the file and create ID (column 4) and the corresponding nodes (column 7).
# Once we have looped through each line of the file, we have stored all IDs and Nodes. Then we can pick out only
# nodes containing edges, and then we need to create the correct format that we want, which is the edge list format.

# This can be done with a for-loop (?) that splits the edges lists and joins each element to the corresponding ID
# attribute. Leaving us with an edge list format where each line is one ID attribute with one corresponding edge
# attribute (?). So we need to use class methods to add instances of the id and edges to the Node class..?

# OR I'm just unfamiliar with objects and classes, and there is no need for a for loop, we can just
# use a method in the Node class and pick out the data we want, and use another method to process the data into an
# edge list format.

"""file_content = file.readlines()
for index, line in enumerate(file_content):
    content_stripped = line.strip("\n")
    columns = content_stripped.split("\t")
    #print(columns)
    if len(columns) != 7:
        print("Bad line format", len(columns), index)
    else:
        edges = columns[6].split(";")
        nodes_and_edges.append(Node(columns[3], edges))"""

"""file_content = file.readlines()
for lines in file_content:
        content_stripped = lines.strip("\n")
        columns.append(content_stripped.split("\t"))
        #print(columns)
for index, line in enumerate(file_content):
    if len(columns) != 7:
        print("Bad line format", len(columns), index)
    else:
        edges = columns[6].split(";")
        nodes_and_edges.append(Node(columns[3], edges))"""


# y = 0
# for x in range(len(columns)):
#     if columns[x] in ("chr:", "-"):
#         test_list_edges.append(columns[y:x])
#         x += 1

#print(edges)
#print(nodes_and_edges)

"""file_content = file.readlines()
for index, lines in enumerate(file_content):
            content_stripped = lines.strip("\n")
            columns.append(content_stripped.split("\t"))
            #real_columns += columns
            #columns.append(content_stripped.split("\t"))
            #print(columns)
            if len(columns) != 7:
                print("Bad line format", len(columns), index)
            else:
                edges = columns[6].split(";")
                nodes_and_edges.append(Node(columns[3], edges))

print(columns[10])"""


#If we use append(content_stripped) to columns , we cant split list. Np.array.split  can split lists but not on delimitors.
#So thats a big problem. Maybe find fix for that.
#Or we need to keep the list as is (why did this work earlier, now I can split the indexed list???) and find a work around
# to keep the globa llst (columns = []) contain all columns crated when to loop is run, not just keeping the last
# line that is assigned to it.
# OG om dette hade funket, så kunne vi ikke splitte lister basert på index, siden hver linje korresponderer til en index, så
# vi må splitte strenger inni hver liste?
# Dette er rotete som faen.

"""if str in real_columns == "chr, -, ;, .":
    newlist.append(what we need)"""




#print("Length of columns outside of loop:", len(columns))




"""for x in range(len(content)):
    subunits = content[x].strip().split()
    new_node = Node(subunits[4], subunits[7])
    nodes.append(Node(new_node))
    if x > len(content):
        break

subunits[0]
"""
# We need to create a new object Node and Edge, which are contained in the Node (?) class.
# Then we need to instanciate objects with a for loop that reads every line of the file, and:
#   - Only picks out column 4 (ID = Node) and column 7 (edges)
#   - Adds ID to ID and all its edges to a list containing ID + one edge
#   - Removes unwanted characters, like newline, ; and headers.

# node1 without edgelist = Node("chr1:900000-1300000", "chr1:154350000-155100000; chr1:228100000-228550000")
# Example of an instance of class 'Node'




# while (reading_from_file_however_that_is_done_in_python) {
#  nodes.append(Node(column4, split_edges_into_a_list(column7)))

# node1 = Node("Chrx:y Chrz:c") Må loope over filen og adde hver node som objekt til Node classen



### Random:

# i = 1
# print(f)
#for line in range(len(f)):
    #if ";" == line[i]:
    #x = re.split(";", f)
    #grid.append(x)
    #else:
    #grid.append(split(";")[i])
    #i += 1

#print(grid)
#ROWS, COLS = len(f), len(f[0])
#for r in range(ROWS):
#    for line in range(COLS):
#        if (str(line)[r]) == ";":
#            grid

#les over hver linje
#for hver linje, velg ut første del før ; og appende biten etter ; til en matrise?
#Så velg ut første del før ; og appende biten mellom ;2 og ;3
"""
Lines er en liste hvor hver linje i fila er en subliste, den kan splittes ved å bruke numpy


#with open("/Users/GBS/Master/HiC-Data/Hi-C_data_fra_Jonas/4linescopy/test_chr18_K562.txt") as f:
#    test18 = f.readlines()
#print(test18)


empty = pd.array((0, 0))

i = 0

#with open("/Users/GBS/Master/HiC-Data/Hi-C_data_fra_Jonas/4linescopy/test_chr18_K562.txt") as f:
   #for line in f:
     #  print(line)
      # if ";" in line:
          # split_list = f.split(";")
#print(split_list)

#with open("/Users/GBS/Master/HiC-Data/Hi-C_data_fra_Jonas/4linescopy/test_chr18_K562.txt") as txt:
#    for lines in txt:
#        if not lines in txt:
#            break
#        txt_split = [i.split(";") for i in txt]

#open(input('\\Users\\GBS\\Master\\HiC-Data\\Hi-C_data_fra_Jonas\\4linescopy\\test_chr18_K562.txt'), "r")


f = open("/Users/GBS/Master/HiC-Data/Hi-C_data_fra_Jonas/4linescopy/test_chr18_K562.tsv")
f = f.read()
grid = []
x = re.split(";", "\n", f)
print(x)

#test rammeverk
#enhetstesting
#python enhetstestingrammeverk
#-- assert at antall reads er 7:
#unit test -- class (liste av stuf med dataen din så du kan eksportere) -- gg
"""

### Lag klasser (@dataclass?), datastruktur og importer filen.


# Eller CSV som dict? Men da må vi konverte hver fil til CSV først da

"""with open("fil", "r") as csv_fil:
    csv_reader = csv.DictReader(csv_fil)

    for line in csv_reader:
        IDcsv = line["ID"] #for å feks få ut alle IDs, siden hver field name blir en key med en value
        print(line["edges"]) #må legge inn: != "."

    # Også outputte det til ny fil? ---- Nei bare lagre ID og edges i et objekt
        #with open("only_nodes_edges", "w") as new_file:
            #csv_writer = csv.DictWriter(new_file, fieldnames = fieldnmaes, delimiter = "", )
"""

# Eller dataclasses?

"""@dataclass
class Node:
    id: str
    node: str
"""
