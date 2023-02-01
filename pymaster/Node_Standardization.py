from dataclasses import dataclass
import pandas as pd
import numpy as np
import pybedtools as pbt
import os as os
from Processing import Node
from Processing import Celline
from Processing import Cellines



def default_path_to_raw_data(*args):
    default_path_to_raw_data_directory = "/Users/GBS/Master/Pipeline/INC-tutorial/hicpro_results/hic_results/matrix/chr18/raw/50000"
    return os.path.join(default_path_to_raw_data_directory, *args)


def default_output_dir(*args):
    default_path_to_output_directory = "/Users/GBS/Master/Pipeline/python_pipe_test/bedpe_testing"
    return os.path.join(default_path_to_output_directory, *args)

# Clean opp senere og lagre matrix som armatus input, kanskje skriv amratus output end-inlcusive function så TAD annotation kan gjøres

def read_hicpro_output():
    """
    Reads raw HiC-Pro output and output dataframes
    Write dataframes to specified directory
    """

    matrix_file = os.path.join(default_path_to_raw_data("chr18_50000.matrix"))
    bed_file = os.path.join(default_path_to_raw_data("chr18_50000_abs.bed"))

    matrix_df = pd.read_csv(matrix_file, sep="\t", header=None)
    bed_df = pd.read_csv(bed_file, sep="\t", header=None)

    def write_to_dir(output_directory):
        output_path_matrix = os.path.join(output_directory, "chr18_50000.matrix")
        output_path_bed = os.path.join(output_directory, "chr18_50000_abs.bed")
        matrix_path = matrix_df.to_csv(output_path_matrix, sep="\t", header=None, index=False)
        bed_path = bed_df.to_csv(output_path_bed, sep="\t", header=None, index=False)

        return matrix_path, bed_path

    # write_to_dir(default_output_dir())

    return matrix_df, bed_df

read_hicpro_output()

# THIS WORKS:

def create_bedpe(bed_file, matrix_file, output_file):
    """
    Function to create a BEDPE file from the HiC-Pro output (BED and matrix)
    and write it to a specified directory
    """

    with open(bed_file, "r") as bed, open(matrix_file, "r") as matrix, open(output_file, "w") as bedpe:
        bed_lines = bed.readlines()
        matrix_lines = matrix.readlines()
        bed_dict = {}
        for line in bed_lines:
            line = line.strip().split("\t")
            bed_dict[line[3]] = line
        for line in matrix_lines:
            line = line.strip().split("\t")
            bedpe.write(f"{bed_dict[line[0]][0]}"
                        f"\t{bed_dict[line[0]][1]}"
                        f"\t{bed_dict[line[0]][2]}"
                        f"\t{bed_dict[line[1]][0]}"
                        f"\t{bed_dict[line[1]][1]}"
                        f"\t{bed_dict[line[1]][2]}"
                        f"\t{line[2]}\n")

    return bedpe

create_bedpe(default_path_to_raw_data("chr18_50000_abs.bed"), default_path_to_raw_data("chr18_50000.matrix"), default_output_dir("standardized_nodes_bedpe"))




# Use pybedtools to make the beads?, and then to blacklist, then statistical tests and then write to gtrack











# function that combines all regions of hg38 that should be excluded from the analysis
# GRC declared regions of hg38 which contain false duplication or contamination
grc_hg38_exclusions = pbt.BedTool("/Users/GBS/Master/reference/GCA_000001405.15_GRCh38_GRC_exclusions.bed")
encode_blacklist = pbt.BedTool("/Users/GBS/Master/reference/ENCFF356LFX.bed")

# find overlap betwen grc_hg38_exclusions and encode_blacklist
overlap_test = pbt.BedTool.intersect(grc_hg38_exclusions, encode_blacklist, wa=True, wb=True)
print(overlap_test)


# Make function that combines all overlapping "blacklisted" regions, which outputs this as a BED file?
# Or abstract this into an object? In either case we want to subtract this from significant interactions later on


def hg38_exclusions():
    """Combine all regions of hg38 that should be excluded from the analysis"""
    # Maybe someting like this but also including blacklisted regions, 5cap, telomeres etc:
    # hg38_exclusions = grc_hg38_exclusions.cat(hg38_centromeres, postmerge=False)
    # return hg38_exclusions


def bin_genome(binsize):
    """Create a 500 kb binned genome, using the reference genome as a template"""

    bin_size = 500000
    bins = reference_bed.window(b=reference_bed, w=bin_size, u=True)
    return bins


# Now we need to create a function that maps the beads to the binned genome

def map_beads_to_bins(beads, bins):
    """Map the beads to the binned genome"""
    mapped_beads = beads.map(b=bins, c=4, o="count")
    return mapped_beads
