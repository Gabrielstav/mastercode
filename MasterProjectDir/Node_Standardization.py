import pybedtools as pbt
from Processing import Node
# 562_all
from Processing import Celline
from Processing import Cellines

HMEC = pbt.BedTool("/Users/GBS/Master/HiC-Data/Hi-C_data_fra_Jonas/4linescopy/HMEC_50kb.domain.RAW.no_cen.NCHG_fdr.o_by_e5_to_plot.gtrack")
K562 = pbt.BedTool("/Users/GBS/Master/HiC-Data/Hi-C_data_fra_Jonas/4linescopy/K562_50kb.domain.RAW.no_cen.NCHG_fdr.o_by_e5_to_plot.gtrack")
IMR90 = pbt.BedTool("/Users/GBS/Master/HiC-Data/Hi-C_data_fra_Jonas/4linescopy/IMR90_50kb.domain.RAW.no_cen.NCHG_fdr.o_by_e5_to_plot.gtrack")
HUVEC = pbt.BedTool("/Users/GBS/Master/HiC-Data/Hi-C_data_fra_Jonas/4linescopy/K562_50kb.domain.RAW.no_cen.NCHG_fdr.o_by_e5_to_plot.gtrack")

K562_c = pbt.BedTool(Cellines.from_default().with_strain("K562").only_con())
HUVEC_c = pbt.BedTool(Cellines.from_default().with_strain("HUVEC").only_con())
reference_bed = pbt.BedTool("/Users/GBS/Master/reference/GCA_000001405.15_GRCh38_GRC_exclusions.bed")
# reference_fastq = pbt.bedtool("/Users/GBS/Master/reference/GRCh38_latest_genomic.fna.gz") need to unzip
# HMEC_500kb = HMEC.intersect(reference_fastq, u=True, f=0.5, header=True)
# print(HMEC)

# maps = K562.map(HMEC, concat=True, c=4, f=0.8, r=True)
# print(maps)

# Trying node overlap for all 4 cellines:
# her er koden i Bedtools:
# bedtools intersect -a HMEC_50kb.domain.RAW.no_cen.NCHG_fdr.o_by_e5_to_plot.gtrack -b HUVEC_50kb.domain.RAW.no_cen.NCHG_fdr.o_by_e5_to_plot.gtrack
# IMR90_50kb.domain.RAW.no_cen.NCHG_fdr.o_by_e5_to_plot.gtrack K562_50kb.domain.RAW.no_cen.NCHG_fdr.o_by_e5_to_plot.gtrack
# -f 0.80 -r -sorted -wa -wb -names HUVEC IMR90 K562 | head -n 1000

# OK at first glance this seems to almost work, but can look into this later?
# loj=True is something?
# where are the names? We still need to know which genomes the nodes that overlap orginate from
# since this is the information we want to end up with in the end somehow (strain + nodes).

# Find overlap between cellines
tester = HUVEC.intersect(HMEC+IMR90+K562, names=HMEC+IMR90+K562, f=0.80, r=True, sorted=True, wa=True, wb=True)
tester.head(20)





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