import pandas as pd
import pybedtools as pbt
import os as os
import bioframe as bf


# OK, make everything into a static class called Pipeline
# Question now is how read in multiple files (from either multiple chromosomes or multiple cellines)

class Pipeline:

    @staticmethod
    def default_rawdata_path(*args):
        default_path_to_raw_data_directory = "/Users/GBS/Master/Pipeline/INC-tutorial/hicpro_results/hic_results/matrix/chr18/raw/50000"
        return os.path.join(default_path_to_raw_data_directory, *args)

    @staticmethod
    def default_output_path(*args):
        default_path_to_output_directory = "/Users/GBS/Master/Pipeline/python_pipe_test/bedpe_testing"
        return os.path.join(default_path_to_output_directory, *args)

    @staticmethod
    def default_reference_path(*args):
        default_path_to_reference_dicrectory = "/Users/GBS/Master/reference"
        return os.path.join(default_path_to_reference_dicrectory, *args)

    @staticmethod
    def read_hicpro_output():
        """
        Reads raw HiC-Pro output and output dataframes
        Write dataframes to specified directory
        """

        matrix_file = os.path.join(Pipeline.default_rawdata_path("chr18_50000.matrix"))
        bed_file = os.path.join(Pipeline.default_rawdata_path("chr18_50000_abs.bed"))

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

    @staticmethod
    def make_bedpe():
        """
        Makes bedpe file from HiC-Pro output
        """
        bedpe = []

        bed_file = os.path.join(Pipeline.default_rawdata_path("chr18_50000_abs.bed"))
        matrix_file = os.path.join(Pipeline.default_rawdata_path("chr18_50000.matrix"))

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

    # @staticmethod
    # def make_bedpe_df():
    #     bedpe_df = Pipeline.make_bedpe()
    #     return bedpe_df

    @staticmethod
    def print_bedpe():
        print(Pipeline.make_bedpe())

    @staticmethod
    def write_bedpe(*args):
        """
        Writes bedpe file to specified directory
        """
        bedpe_file = os.path.join(Pipeline.default_output_path(*args))
        with open(bedpe_file, "w") as f:
            f.writelines(Pipeline.make_bedpe())

    @staticmethod
    def blacklisted_regions():
        """
        Returns blacklisted regions
        """
        # Path to blacklist file
        # Use pandas to read in file
        # Use pybedtools or bioframe to calculate overlap
        # Add check to see if there is any overlap, if there is, remove whole bin (meaning whole entry?)
        # Return new bedpe file

        blacklisted = Pipeline.default_reference_path("hg19/unmappable_blacklist.bed")
        blacklisted_regions = pbt.BedTool(blacklisted)
        bedpe_pbt = pbt.BedTool(Pipeline.make_bedpe())
        overlap = blacklisted_regions.intersect(bedpe_pbt)
        overlap_df = overlap.to_dataframe()
        return overlap

    @staticmethod
    def write_overlap(*args):
        # overlap = os.path.join(Pipeline.default_output_path(*args))
        # with open(overlap, "w") as f:
        #     f.writelines(Pipeline.blacklisted_regions())

        out = Pipeline.blacklisted_regions().saveas("/Users/GBS/Master/Pipeline/python_pipe_test/bedpe_testing/overlap.bed")
        return out


    @staticmethod
    def print_bedpe_pbt():
        bedpe_pbt = pbt.BedTool(Pipeline.make_bedpe())
        print(bedpe_pbt)


    @staticmethod
    def cytobands():
        """
        returns cytoband regions
        """
        pass

    @staticmethod
    def repeats_or_something():
        """
        Returns unmappable regions
        """
        pass

    @staticmethod
    def remove_unmappable_regions():
        """
        Removes unmappable regions from bedpe file combining all blacklisted
        regions into one file and finding overlap with BEDPE file.
        Any overlap removes the whole interaction (bin), retaining equal bin size acros the genome.
        """

        pass

    @staticmethod
    def find_siginificant_interactions():
        """
        Find significant interactions using a p-value cutoff using the non-central hypergeometric ditribution
        """
        pass

    @staticmethod
    def adjust_pvalues():
        """
        Adjust p-values using the Benjamini-Hochberg method
        """
        pass

    @staticmethod
    def make_gtrack():
        """
        Make gtrack file from bedpe file with significant interactions
        """
        pass


# Pipeline.write_overlap("testing_overlap")



# 1. Create BEDPE file (done)
# 2. Remove regions overlapping blacklisted regions
# 3. Run NCHG script to find significant interactions
# 4. Run FDR script to correct for multiple testing
# 5. Make GTrack format (unclear if I can use one file (beads))
# 6. Make the Pipeline able to read in multiple files (cell lines) and process them


# Blacklisted regions:

# 1. Make reference dir containing:
    # A. Download and store hg38 and hg19 blacklisted regions (find out how reliable they are), make a function
    # B. Downlaod and store reference genome (hg38 and hg19)
    # C. Download and store chromosome sizes (hg38 and hg19)
# 2. Combina all the different unmappable regions into one BED file for hg19 and one for hg38?
# 3. Take the BEDPE format and find overlap with the correct blacklisted regions (for the same ref genome) using pybedtools
# 4. Write the output to a BEDPE file/BED file? Check.





# Now write logic to see if the overlap overlaps with the bedpe file, if it does, remove the whole bin (interaction)


def remove_overlap():

    blacklisted = open("/Users/GBS/Master/reference/hg19/unmappable_blacklist.bed")
    blacklised_pbt = pbt.BedTool(blacklisted)
    bedpe_pbt = pbt.BedTool(Pipeline.make_bedpe())
    # overlap = bedpe_pbt.intersect(blacklised_pbt, r=True, v=False, w=50000)

    overlap = bedpe_pbt.window(blacklised_pbt, w=50000, r=False, v=True)
    pbt.BedTool(overlap).saveas("/Users/GBS/Master/Pipeline/python_pipe_test/bedpe_testing/overlap.bed")
    return print(overlap), print(type(overlap))

remove_overlap()

# I think this actually work, just by looking at the output file.
# To verify, do the reverse, output only the overlap and use diff to see if the output is the same as the overlap file.
# But this looks promising.



