import pybedtools as pbt
from Processing import Node
from Processing import Celline
from Processing import Cellines



# 4 cell lines

# K562_c = pbt.BedTool(Cellines.from_default().with_strain("K562").only_con())
# This doesn't work, need to write custom methods for representing our instances
# such that iGraph understands this, if it's possible. Low pri RN.

HMEC = pbt.BedTool("/Users/GBS/Master/HiC-Data/Hi-C_data_fra_Jonas/4linescopy/HMEC_50kb.domain.RAW.no_cen.NCHG_fdr.o_by_e5_to_plot.gtrack")
K562 = pbt.BedTool("/Users/GBS/Master/HiC-Data/Hi-C_data_fra_Jonas/4linescopy/K562_50kb.domain.RAW.no_cen.NCHG_fdr.o_by_e5_to_plot.gtrack")
IMR90 = pbt.BedTool("/Users/GBS/Master/HiC-Data/Hi-C_data_fra_Jonas/4linescopy/IMR90_50kb.domain.RAW.no_cen.NCHG_fdr.o_by_e5_to_plot.gtrack")
HUVEC = pbt.BedTool("/Users/GBS/Master/HiC-Data/Hi-C_data_fra_Jonas/4linescopy/K562_50kb.domain.RAW.no_cen.NCHG_fdr.o_by_e5_to_plot.gtrack")

# Reference genome (change dir location to ref or something later)
reference_bed = pbt.BedTool("/Users/GBS/Master/reference/GCA_000001405.15_GRCh38_GRC_exclusions.bed")

# Segmented 500 bins (file createed by bedtools)
bins = pbt.BedTool("/Users/GBS/Master/testing_node_standardization/chrom_hg19_test.sizes")

# But I need to crete the 500 kb bins here

# So make ethe 500 kb binned genome in pbt
# Make objects of all cell lines (do this in a function that calls the Processing module, automatically making the pbt bedtool objects)
# Next we need to write a mapping function, that loops over/uses lambda functions to
#   - grab each node (position, length, chromosome) and hold this info
#   - compares the position of the node with the binned genome corresponding to the position of the node
#   - then comes the logic of the function;
#       -

# BEDmap to map between the BED objects? Map has multiple statistical operations, check BEDtools documentation.
# 


# # Find overlap between cellines (obsolete?)
# tester = HUVEC.intersect(HMEC+IMR90+K562, names=HMEC+IMR90+K562, f=0.80, r=True, sorted=True, wa=True, wb=True)
# tester.head(20)





""" 
Node read length statistics 
"""

# We need to create positional information from the reads in the nodes, so we can sort them e.g by length and
# plot them from start to finish or something? Later we want to integrate RNAseq data so we can overlay this on the nodes?

# Use pybedtools or do it organically in python, think we should avoid hte command line if we can.

# Can we just write a function in python that takes each connected node as input from the data class. The nodes are already ordered from start to end of chromosome/genome.
# Then we calculate the length of the first node for the 4 different cell lines.
# Then we quantify their positions, by using the start/end point.
# Maybe the largest node is the anchor node, and the smaller nodes are matched against it?
# After we have length and positions, we can calculate the percentage of overlap.
# If nodes overlap less than say 50%, the node is discarded (not deleted from the data class, but just not included in the function output).

# We want to examine node overlap both on the "node level" (whatever that means) and on the bp level.
# Find overlap using pybadtools intersect? or just write the function to calculate the percentage overlap/ bp overlap.
# Do this for every node, iterating through each instance of the data class for the 4 cell lines.
# Before we do this we need to write functionality for reading in multiple cell lines and store their data as instances of the Node class.
# Then we can iterate over the instances, for each cell line, all the way through the data class.
# Inside this iterator, we need to do the percentage/bp/node overlap calculations
# OR
# We can output the mapping of the nodes (start and end to find length) for each isntance to a dict or something,
# and then use this for further coverage calculations.