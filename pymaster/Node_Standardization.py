import pandas as pd
import pybedtools as pbt
import os as os
import numpy as np
import scipy as sp
import math
import sys
from statsmodels.sandbox.stats import multicomp


# OK, make everything into a static class called Pipeline
# Question now is how read in multiple files (from either multiple chromosomes or multiple cellines)
# Abstract away the file reading and writing to separate methods later

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
    def write_blacklist_hg19(*args):
        Pipeline.remove_blacklist_hg19().saveas(Pipeline.default_output_path(*args))

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
        """
        This is the input to the NCHG script
        """
        Pipeline.remove_cytobands().saveas(Pipeline.default_output_path(*args))


    @staticmethod
    def remove_cap():
        """
        Sets a cap on the length of the chromosomes, any interactions longer than this are removed
        this is only applicable to interactions between chromosomes, this ensures that the interactions
        do not go beyond the end of the chromosome? I think. But we have a constant bin size, we do not mix
        bin sizes, so this is not necessary maybe? IDK skip for now. Maybe this?:

        chrsize = {}
        with open(chrsize_file) as f:
            for l in f:
                chrsize[l.split()[0]] = int(float(l.split()[1]))

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

# Then I think we just need to use the NCHG script from the terminal for now
# and write to a dir, then use python again to adjust pvalues and make the GTrack format

    @staticmethod
    def find_siginificant_interactions():
        """
        Find significant interactions using a p-value cutoff using the non-central hypergeometric ditribution
        Seems hard to port this to python, source is in C++ and I don't want to spend time on that right now
        """
        pass

    @staticmethod
    def adjust_pvalues():
        """
        Adjust p-values using the Benjamini-Hochberg method
        """

        NCHG = os.path.join(Pipeline.default_output_path("NCHG_chr18_test1"))

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
        Make this later?
        """
        pass


# Pipeline.remove_cytobands()
# Pipeline.write_remove_cytobands("testing_cytob")

# exec(open("/Users/GBS/Master/Pipeline/INC-tutorial/processing_scripts/NCHG_hic/NCHG").read())


# 1. Create BEDPE file (done)
# 2. Remove regions overlapping blacklisted regions (done
# 3. Run NCHG script to find significant interactions
# 4. Run FDR script to correct for multiple testing
# 5. Make GTrack format (unclear if I can use one file (beads))
# 6. Make the Pipeline able to read in multiple files (cell lines) and process them


Pipeline.adjust_pvalues()