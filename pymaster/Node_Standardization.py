import pandas as pd
import pybedtools as pbt
import os as os
import subprocess as sp
import tempfile as tf
import numpy as np
import scipy as scp
import math
import sys
import statsmodels.api as sm
from statsmodels.sandbox.stats import multicomp

# Pre-processing pipline for Hi-C data:

# Output by default should be files containing edge lists, named by chromosome/celline and resolution
# Optionally can also output GTrack files with the same naming convention with or without p-values (for weights)

# TODO:
# 1. Make automatic file naming, reading and writing in separate module, that feeds the pipeline module with the correct files
# 2. Make Pipeline handle Hg38 as well, and make it so that the user can specify which reference genome to use and the methods adapt accordingly
# 3. Make GTRack methods with- and without p-values and map between GTrack and edge list files for visualization purposes



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
    def path_to_NCHG():
        return "/Users/GBS/Master/Pipeline/INC-tutorial/processing_scripts/NCHG_hic/NCHG"

    @staticmethod
    def window_size():
        return 50000

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

    @staticmethod
    def print_bedpe():
        print(Pipeline.make_bedpe())

    @staticmethod
    def print_bedpe_pbt():
        bedpe_pbt = pbt.BedTool(Pipeline.make_bedpe())
        print(bedpe_pbt)

    @staticmethod
    def write_bedpe(*args):
        """
        Writes bedpe file to specified directory
        """
        bedpe_file = os.path.join(Pipeline.default_output_path(*args))
        with open(bedpe_file, "w") as f:
            f.writelines(Pipeline.make_bedpe())

    @staticmethod
    def remove_blacklist_hg19():

        """
        Removes blacklisted regions from bedpe file
        This blacklist is from the ENCODE blacklist of problematic regions in hg19
        https://github.com/Boyle-Lab/Blacklist

        """

        blacklisted = os.path.join(Pipeline.default_reference_path("hg19/hg19-blacklist.v2.bed"))
        blacklised_pbt = pbt.BedTool(blacklisted)
        bedpe_pbt = pbt.BedTool(Pipeline.make_bedpe())

        no_overlap_bedpe = bedpe_pbt.window(blacklised_pbt, w=Pipeline.window_size(), r=False, v=True)

        return no_overlap_bedpe

    @staticmethod
    def remove_blacklist_hg38():

        """
        Removes blacklisted regions from bedpe file
        This blacklist is from the ENCODE blacklist of problematic regions in hg38
        """

        blacklisted = os.path.join(Pipeline.default_reference_path("hg38 blacklist path"))
        blacklised_pbt = pbt.BedTool(blacklisted)
        bedpe_pbt = pbt.BedTool(Pipeline.make_bedpe())

        no_overlap_bedpe = bedpe_pbt.window(blacklised_pbt, w=Pipeline.window_size(), r=False, v=True)

        return no_overlap_bedpe

        pass

    @staticmethod
    def write_blacklist_hg19(*args):
        Pipeline.remove_blacklist_hg19().saveas(Pipeline.default_output_path(*args))

    @staticmethod
    def write_blacklist_hg38(*args):
        Pipeline.remove_blacklist_hg38().saveas(Pipeline.default_output_path(*args))

    # Instead of writing a new function for each blacklist, I could write a check for the genome version

    @staticmethod
    def remove_cytobands():
        """
        Cytoband locations are determined in this case by Giemsa staining (I think)
        and are located and removed from the BEDPE file
        """

        cytobands = os.path.join(Pipeline.default_reference_path("hg19/cytoBand_hg19.txt"))

        centromeric_regions = []
        with open(cytobands, "r") as f:
            cytobands = f.readlines()
            for line in cytobands:
                line = line.strip().split("\t")
                if line[4] == "acen":
                    centromeric_regions.append(line[0:5])

        centromeric_regions_pbt = pbt.BedTool(centromeric_regions)
        no_cytobands_bedpe = Pipeline.remove_blacklist_hg19().window(centromeric_regions_pbt, w=Pipeline.window_size(), r=False, v=True)

        return no_cytobands_bedpe

    @staticmethod
    def print_remove_cytobands():
        return print(Pipeline.remove_cytobands())

    @staticmethod
    def write_remove_cytobands(*args):

        with open(Pipeline.default_output_path(*args), "w") as f:
            f.writelines(Pipeline.remove_cytobands().to_dataframe().to_csv(sep="\t", header=False, index=False))

    @staticmethod
    def remove_cap():

        # Refactor to work with pbt in pipeline

        """
        Sets a cap on the length of the chromosomes, any interactions longer than this are removed
        this is only applicable to interactions between chromosomes, this ensures that the interactions
        do not go beyond the end of the chromosome? I think. But we have a constant bin size, we do not mix
        bin sizes, so this is not necessary maybe? IDK skip for now. Refactor this to work with pbt in pipeline?:

        chrsize = {}
        with open(chrsize_file) as f:
            for length in f:
                chrsize[length.split()[0]] = int(float(length.split()[1]))

        for n in sys.stdin:
            n = n.split()
            if int(n[2]) > chrsize[n[0]]:
                n[2] = chrsize[n[0]]
            if int(n[5]) > chrsize[n[3]]:
                n[5] = chrsize[n[3]]
            n = map(str,n)
            print('\t'.join(n))
        """
        pass

    @staticmethod
    def input_to_nchg():
        """
        This is the input to the NCHG script
        """
        os.chdir(Pipeline.default_output_path())
        file = Pipeline.default_output_path("nchg_input.txt")
        with open(file, "w") as f:
            f.writelines(Pipeline.remove_cytobands().to_dataframe().astype(str).to_csv(sep="\t", header=False, index=False))
        return file

    @staticmethod
    def find_siginificant_interactions():
        """
        NCHG script to calculate the significance of interactions:
        m = minimum interaction length in bp, should be same as window size used to make the bedpe file
        p = input file, which is the output of the remove_cytobands function but reformatted to be compatible with NCHG
        """

        nchg_run = sp.run([Pipeline.path_to_NCHG(), "-m", str(Pipeline.window_size()), "-p", Pipeline.input_to_nchg()])

        nchg_out = nchg_run.stdout
        return nchg_out

    @staticmethod
    def write_sig_interactions(*args):
        file = os.path.join(Pipeline.default_output_path(*args))
        with open(file, "w") as f:
            f.writelines(Pipeline.find_siginificant_interactions())

    @staticmethod
    def adjust_pvalues():
        """
        Adjust p-values using the Benjamini-Hochberg method
        """

        # Take output from NCHG
        # set method for mutliple testing correction (BH), log ratio and FDR threshold
        # Read in NCHG output and strip on new line, store in list
        # Compute the log ratio for each line in list, add it to line as new field the log ratio is the line[9] / line[10] in base 2
        # Use statmodels multicomp to adjust p-values
        # Filter p values based on FDR threshold and log ratio (less than FDR, same or greater than log ratio)
        # Return the filtered lines with only the edge and adjusted p-value

        # NCHG = os.path.join(Pipeline.default_output_path("NCHG_chr18_test1"))

        # Dict where keys are chromosome names and values chromosome sizes
        # Then iterate over NCHG file

        # # Get file name
        # file_name = input("Enter file name: ")
        #
        # with open(file_name) as infileH:
        #     infile = infileH.readlines()
        #
        # # Get method
        # method = input("Enter method: BH")
        # # Get log ratio
        # log_ratio = float(input("Enter log ratio: "))
        # # Get FDR threshold
        # fdr_thres = float(input("Enter FDR threshold: "))
        #
        # infile = [x.strip() for x in infile]
        #
        # pval = []
        # padj = []
        # processline = []
        # for line in infile:
        #     line = line.split()
        #     pval.append(float(line[6]))
        #     if float(line[9]) == 0 or float(line[10]) == 0:
        #         logratio = 0.0
        #     else:
        #         logratio = math.log(float(line[9]), 2) - math.log(float(line[10]), 2)
        #
        #     line = " ".join(line) + " " + str(logratio)
        #     processline.append(line)
        #
        # padj = list(multicomp.multipletests(pval, method=method))
        #
        # for i in range(len(processline)):
        #     line = processline[i] + " " + str(padj[1][i])
        #     line = line.split()
        #     if float(line[11]) >= log_ratio and float(line[12]) <= fdr_thres:
        #         print(line[0] + "\t" + line[1] + "\t" + line[2] + "\t" + line[3] + "\t" + line[4] + "\t" + line[5] + "\t" + line[12])

    @staticmethod
    def make_gtrack():
        """
        Make gtrack file from bedpe file with significant interactions
        Make this able to output a gtrack file with weights (p-values) and without weights (p-values)
        Make the pipeline output edge list by default, but optionally GTrack with and without weights
        Also refactor processing to be Gtrack to edge list converter class
        """

        pass

    @staticmethod
    def make_edge_list():
        """
        Make edge list from bedpe file with significant interactions and padj
        """

        padj = os.path.join(Pipeline.default_output_path("FDR_chr18_test1"))  # This is the output from the FDR script
        # NCHG = os.path.join(Pipeline.default_output_path("NCHG_chr18_test1"))  # This is the output from the NCHG script

        edge_list_nop = []
        edge_list_withp = []

        with open(padj) as file:
            for line in file:
                if len(line) == 6:
                    line = line.split()
                    edge_list_nop.append(line[0] + ":" + line[1] + "-" + line[2] + " " + line[3] + "-" + line[4] + ":" + line[5])
                else:
                    line = line.split()
                    edge_list_withp.append(line[0] + ":" + line[1] + "-" + line[2] + " " + line[3] + "-" + line[4] + ":" + line[5] + " " + line[6])

        return edge_list_nop, edge_list_withp

    # this is nooby noob, make better!!!:
    # Split withp and nop into separate methods
    @staticmethod
    def make_edgelist_withp():
        with open(os.path.join(Pipeline.default_output_path("FDR_chr18_test1"))) as file:
            for line in file:
                print(line)

    @staticmethod
    def write_edge_list():

        edge_list_nop, edge_list_withp = Pipeline.make_edge_list()

        with open("edge_list_nop.txt", "w") as file:
            for line in edge_list_nop:
                file.write(line + "\n")

        with open("edge_list_withp.txt", "w") as file:
            for line in edge_list_withp:
                file.write(line + "\n")

    @staticmethod
    def print_edge_list_withp():

        edge_list_withp = Pipeline.make_edge_list()
        for line in edge_list_withp:
            print(line)

    @staticmethod
    def print_edge_list_nop():

        edge_list_nop = Pipeline.make_edge_list()
        for line in edge_list_nop:
            print(line)




# print(type(Pipeline.remove_cytobands()))
# print(type(Pipeline.NCHG_input()))
# print(Pipeline.NCHG_input())


Pipeline.find_siginificant_interactions()


# Pipeline.remove_cytobands()
# Pipeline.write_remove_cytobands("testing_cytob")

# exec(open("/Users/GBS/Master/Pipeline/INC-tutorial/processing_scripts/NCHG_hic/NCHG").read())


# 1. Create BEDPE file (done)
# 2. Remove regions overlapping blacklisted regions (done)
# 3. Run NCHG script to find significant interactions (done)
# 4. Run FDR script to correct for multiple testing (done)
# 5. Make GTrack format (unclear if I can use one file (beads) nope, need two files)
# 6. Make the Pipeline able to read in multiple files (cell lines) and process them

# finish code, clean it up bby, make it nice
# wrap the executable, so it runs here
# make the class able to read in data from multiple cell lines
# make network from test data
# make code run as standalone script, either with sys.argv, argparse or input
# make ylm file for the pipeline

# automatisk les inn filer: run N-times where N is the number of files in the folder (or a specified number) and run the pipeline on each file
# ha check for at den ene filen er ferdig processert, før den neste starter
