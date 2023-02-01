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


# Ok, so this is the first pre-processing step that needs to be done after running HiC-Pro.
# We need to take the output from HiC-Pro and aggregate the contacts between the HiC interactions
#


# The old way of using the GTrack files with differing bin sizes:
# K562_c = pbt.BedTool(Cellines.from_default().with_strain("K562").only_con())
# This doesn't work, need to write custom methods for representing our instances
# such that iGraph understands this, if it's possible. Low pri RN.
# HMEC = pbt.BedTool("/Users/GBS/Master/HiC-Data/Hi-C_data_fra_Jonas/4linescopy/HMEC_50kb.domain.RAW.no_cen.NCHG_fdr.o_by_e5_to_plot.gtrack")
# K562 = pbt.BedTool("/Users/GBS/Master/HiC-Data/Hi-C_data_fra_Jonas/4linescopy/K562_50kb.domain.RAW.no_cen.NCHG_fdr.o_by_e5_to_plot.gtrack")
# IMR90 = pbt.BedTool("/Users/GBS/Master/HiC-Data/Hi-C_data_fra_Jonas/4linescopy/IMR90_50kb.domain.RAW.no_cen.NCHG_fdr.o_by_e5_to_plot.gtrack")
# HUVEC = pbt.BedTool("/Users/GBS/Master/HiC-Data/Hi-C_data_fra_Jonas/4linescopy/K562_50kb.domain.RAW.no_cen.NCHG_fdr.o_by_e5_to_plot.gtrack")

# We want to use Pybedtools to create a 500 kb binned genome, and then map the nodes to this binned genome
# Use tutorial and files from INC to reverse engineer the correct format.
# We take HiC-Pro, aggregate interactions into beads using some string comprehension, which in the INC tut is done using awk
# but I want to do this in python instead, use BED and matrix output to create BEDPE format.
# One of these steps in the tutorial creates BEDPE files that have varying bin sizes, it is unclear if this is due to
# Armatus (where the bin size specifies is the same as the bin size in HiC-Pro) or if this is due to Bedtools complementing the gaps,
# giving varying regions of interaction/no interaction, and then later these are used as input to the statistical scripts, cuasing the
# non-standardized node lengths.



