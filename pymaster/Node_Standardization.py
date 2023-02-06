import pathlib
from dataclasses import dataclass
import pandas as pd
import numpy as np
import pybedtools as pbt
import os as os



# I could see all of this being in a class called "Pre-processing"


# OK, make everything into a static class called Pipeline

class Pipeline:

    @staticmethod
    def default_path_to_raw_data(*args):
        default_path_to_raw_data_directory = "/Users/GBS/Master/Pipeline/INC-tutorial/hicpro_results/hic_results/matrix/chr18/raw/50000"
        return os.path.join(default_path_to_raw_data_directory, *args)

    @staticmethod
    def default_path_to_output(*args):
        default_path_to_output_directory = "/Users/GBS/Master/Pipeline/python_pipe_test/bedpe_testing"
        return os.path.join(default_path_to_output_directory, *args)

    @staticmethod
    def read_hicpro_output():
        """
        Reads raw HiC-Pro output and output dataframes
        Write dataframes to specified directory
        """

        matrix_file = os.path.join(Pipeline.default_path_to_raw_data("chr18_50000.matrix"))
        bed_file = os.path.join(Pipeline.default_path_to_raw_data("chr18_50000_abs.bed"))

        matrix_df = pd.read_csv(matrix_file, sep="\t", header=None)
        bed_df = pd.read_csv(bed_file, sep="\t", header=None)

        return matrix_df, bed_df

    @staticmethod
    def print_hicpro_output():
        """
        Prints raw HiC-Pro output
        """
        matrix_df, bed_df = Pipeline.read_hicpro_output()
        print(matrix_df)
        print(bed_df)


Pipeline.print_hicpro_output()


# @dataclass(frozen=True, eq=True) #doesnt make sense to use dataclasses here
# class Pipeline:
#
#
#     @classmethod
#     def from_default_raw(cls):
#         return cls.from_dir("/Users/GBS/Master/Pipeline/INC-tutorial/hicpro_results/hic_results/matrix/chr18/raw/50000")
#
#     @classmethod
#     def from_default_output(cls):
#         return cls.from_dir("/Users/GBS/Master/Pipeline/python_pipe_test/bedpe_testing")
#
#     @classmethod
#     def from_dir(cls, param):
#         pass
#
#     @classmethod
#     def instantiate(cls):
#         for file in os.listdir(self.rawdata):
#             if file.endswith(".matrix"):
#                 matrix_file = os.path.join(self.rawdata, file)
#             elif file.endswith(".abs.bed"):
#                 bed_file = os.path.join(self.rawdata, file)
#         return Pipeline(matrix_file, bed_file)
#
#     def read_hic_output(self):
#         matrix_file = os.path.join(self.rawdata, "chr18_50000.matrix")
#         bed_file = os.path.join(self.rawdata, "chr18_50000_abs.bed")
#
#         matrix_df = pd.read_csv(matrix_file, sep="\t", header=None)
#         bed_df = pd.read_csv(bed_file, sep="\t", header=None)
#
#         return matrix_df, bed_df
#
#     def write_hic_output(self):
#         matrix_df, bed_df = self.read_hic_output()
#         matrix_path = matrix_df.to_csv(os.path.join(self.output, "chr18_50000.matrix"), sep="\t", header=None, index=False)
#         bed_path = bed_df.to_csv(os.path.join(self.output, "chr18_50000_abs.bed"), sep="\t", header=None, index=False)
#
#         return matrix_path, bed_path
#
#     def make_bedpe(self):
#         matrix_df, bed_df = self.read_hic_output()
#         bedpe_df = pd.DataFrame()
#         bedpe_df["chr1"] = bed_df[0]
#         bedpe_df["start1"] = bed_df[1]
#         bedpe_df["end1"] = bed_df[2]
#         bedpe_df["chr2"] = bed_df[0]
#         bedpe_df["start2"] = bed_df[1]
#         bedpe_df["end2"] = bed_df[2]
#         bedpe_df["score"] = matrix_df[3]
#
#         return bedpe_df
#
#     def write_bedpe(self):
#         bedpe_df = self.make_bedpe()
#         bedpe_path = bedpe_df.to_csv(os.path.join(self.output, "chr18_50000.bedpe"), sep="\t", header=None, index=False)
#
#         return bedpe_path
#



def default_path_to_raw_data(*args):
    default_path_to_raw_data_directory = "/Users/GBS/Master/Pipeline/INC-tutorial/hicpro_results/hic_results/matrix/chr18/raw/50000"
    return os.path.join(default_path_to_raw_data_directory, *args)


def default_output_dir(*args):
    default_path_to_output_directory = "/Users/GBS/Master/Pipeline/python_pipe_test/bedpe_testing"
    return os.path.join(default_path_to_output_directory, *args)

def hg38_blacklist(*args):
    grc_hg38_exclusions = pbt.BedTool("/Users/GBS/Master/reference/GCA_000001405.15_GRCh38_GRC_exclusions.bed")
    encode_blacklist = pbt.BedTool("/Users/GBS/Master/reference/ENCFF356LFX.bed")
    combined_blacklist = pbt.Bedtools(grc_hg38_exclusions, encode_blacklist, )
    return grc_hg38_exclusions, encode_blacklist

def hg19_blacklist(*args):
    """ ENCODE hg19/Grc37 blacklisted regions """




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
# But needs to make it able to store/process multiple cell lines and feed them into the "Processing" classes
# Clean opp senere og lagre matrix som armatus input, kanskje skriv amratus output end-inlcusive function så TAD annotation kan gjøres


# Old function:

# def create_bedpe(bed_file, matrix_file, output_file):
#     """
#     Function to create a BEDPE file from the HiC-Pro output (BED and matrix)
#     and write it to a specified directory
#     """
#
#     # Split the function: One part for processing, another for writing to file
#
#     with open(bed_file, "r") as bed, open(matrix_file, "r") as matrix, open(output_file, "w") as bedpe:
#         bed_lines = bed.readlines()
#         matrix_lines = matrix.readlines()
#         bed_dict = {}
#         for line in bed_lines:
#             line = line.strip().split("\t")
#             bed_dict[line[3]] = line
#         for line in matrix_lines:
#             line = line.strip().split("\t")
#             bedpe.write(f"{bed_dict[line[0]][0]}"
#                         f"\t{bed_dict[line[0]][1]}"
#                         f"\t{bed_dict[line[0]][2]}"
#                         f"\t{bed_dict[line[1]][0]}"
#                         f"\t{bed_dict[line[1]][1]}"
#                         f"\t{bed_dict[line[1]][2]}"
#                         f"\t{line[2]}\n")
#
#     return bedpe


def make_BEDPE():

    """This function takes the BED and matrix files from HiC-Pro and makes a BEDPE object"""

    bedpe = []

    bed_file = os.path.join(default_path_to_raw_data("chr18_50000_abs.bed"))
    matrix_file = os.path.join(default_path_to_raw_data("chr18_50000.matrix"))

    bed_lines = open(bed_file, "r").readlines()
    matrix_lines = open(matrix_file, "r").readlines()

    bed_dict = {}
    for line in bed_lines:
        line = line.strip().split("\t")
        bed_dict[line[3]] = line
    for line in matrix_lines:
        line = line.strip().split("\t")
        bedpe.append(f"{bed_dict[line[0]][0]}"
                     f"\t{bed_dict[line[0]][1]}"
                     f"\t{bed_dict[line[0]][2]}"
                     f"\t{bed_dict[line[1]][0]}"
                     f"\t{bed_dict[line[1]][1]}"
                     f"\t{bed_dict[line[1]][2]}"
                     f"\t{line[2]}\n")

    return bedpe

def write_bedpe(bedpe, output_file):

    """This function writes the BEDPE list to a specified directory"""

    with open(output_file, "w") as bedpe_file:
        for line in bedpe:
            bedpe_file.write(line)

    return bedpe_file

write_bedpe(make_BEDPE(), default_output_dir("standardized_nodes_bedpe_writingtest1"))

# This works, but I want to make it into a class, so I can make its better, but later when I have done comp-net-analysis




# This all means that we have the correct thought process, and we finally understand, it is revealed to my by myself.
# So to summarize:
# 1. Create BEDPE file (done)
# 2. Remove regions overlapping blacklisted regions
# 3. Run NCHG script to find significant interactions
# 4. Run FDR script to correct for multiple testing
# 5. Make GTrack format (unclear if I can use one file (beads))

# So, first we need to make a function that removes regions overlapping blacklisted regions


# 1. Download and store hg38 and hg19 blacklisted regions (find out how reliable they are), make a function
# 2. Integrate the blacklisted script (look at it, is it open source? Also find out how reliable it is)
# 3. Downlaod and store reference genome (hg38 and hg19) and chromosome sizes (hg38 and hg19)
# 4. Take the BEDPE format and find overlap with the correct blacklisted regions (for the same ref genome)
# 5. Make a function that removes the overlap and outputs a BEDPE file




def remove_blacklisted_regions(bedpe_file, output_file):
    """This function takes the combined regions of the relevant reference genome
    containing unmappable regions, caps, centromeres etc. and removes the parts of the
    BEDPE file that overlap with these regions. The output is a BEDPE file"""




# function that combines all regions of hg38 that should be excluded from the analysis
# GRC declared regions of hg38 which contain false duplication or contamination
# We also need Hg19 (if that is the genome we are working with, which in the Jonas data it is at least I tihnk)
# grc_hg38_exclusions = pbt.BedTool("/Users/GBS/Master/reference/GCA_000001405.15_GRCh38_GRC_exclusions.bed")
# encode_blacklist = pbt.BedTool("/Users/GBS/Master/reference/ENCFF356LFX.bed")

# find overlap betwen grc_hg38_exclusions and encode_blacklist
# overlap_test = pbt.BedTool.intersect(grc_hg38_exclusions, encode_blacklist, wa=True, wb=True)
# print(overlap_test)


# Make function that combines all overlapping "blacklisted" regions, which outputs this as a BED file?
# Or abstract this into an object? In either case we want to subtract this from significant interactions later on


def hg38_exclusions():
    """Combine all regions of hg38 that should be excluded from the analysis"""
    # Maybe someting like this but also including blacklisted regions, 5cap, telomeres etc:
    # hg38_exclusions = grc_hg38_exclusions.cat(hg38_centromeres, postmerge=False)
    # return hg38_exclusions


def bin_genome(binsize):
    """Create a 500 kb binned genome, using the reference genome as a template"""

    # bin_size = 500000
    # bins = reference_bed.window(b=reference_bed, w=bin_size, u=True)
    # return bins


# Now we need to create a function that maps the beads to the binned genome

def map_beads_to_bins(beads, bins):
    """Map the beads to the binned genome"""
    mapped_beads = beads.map(b=bins, c=4, o="count")
    return mapped_beads
