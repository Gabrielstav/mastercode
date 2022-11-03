# Find out if we can do the Hi-C processing here instead of using the command line

# NOTES ATM

# I don't think we can go straight from Hi-C data to an edge list we cancrate networks from.
# There is soemthing deep I do not understand about this process, why cant we just take the #
# contact map and matrix from the HiC output and generate networks with this?
# WE NEED to find out why, beacuse this is at the core of a piece of understanding we need to have.
# What in the pipeline ADDS information to this process, we know Armatus does, but what exactly?

# The first three files in the pipeline are parts of the HI-C processing python pipeline scripts,
# the node connectedness needs to be integrated into the processin part, in the classes, so we can output
# connectedness. Another core problem is outputting data from the dataclasses that we can use further on,
# for plotting different things and for outputting edge lists. I think iGraph needs to import edge lists
# from files, so we need to integrate a "write to dir" method in the classes that does this automatically.
# Does this mean that we need make methods that output data for every single type, say if we need a plot
# we need to write methods for this? Perhaps a better way to do this would be to write methods that outputs
# data in general formats: On Nodes, NodeList, Celline(s) levels = list of nodes, dictionary of nodes,
# edge list of nodes, and we can call what class methods we want to output in these formats.
# Do we also need to write custom iterators for out classes?




