import pathlib

import pandas as pd
import pybedtools as pbt
import os as os
import subprocess as sp
import math
from statsmodels.sandbox.stats import multicomp
from File_Handler import FileHandler



# Pre-processing pipline for Hi-C data:

# Output by default should be files containing edge lists, named by chromosome/celline and resolution
# Optionally can also output GTrack files with the same naming convention with or without p-values (for weights)

# TODO:
# 1. Make automatic file naming, reading and writing in separate module, that feeds the pipeline module with the correct files
# 2. Make Pipeline handle Hg38 as well, and make it so that the user can specify which reference genome to use and the methods adapt accordingly
# 3. Make GTRack methods with- and without p-values and map between GTrack and edge list files for visualization purposes
# 4. Gtrack processing methods allow for selection of specific chromosomes, and also allow for selection of specific cell lines, we want this functionality in the edge list methods as well
# but without converting to GTrack first, so either make the methods that does this in the Pipeline class with a flag to select between GTrack and edge list, or make a separate class for edge list processing
# so we do not need to convert to GTrack first. Also include/refactor GTrack methods into Pipeline class.
# 5. Export p values and padj for plotting

# Maybe best to make the Processing classes handle both GTrack and edge list files, and then make a separate class for file handling that handles the reading and writing of files.





# TODO: Call this class from main later on, and make it so that the user can specify input and output dirs

class SetDirectories:

    """
    SET INPUT-, OUTPUT- AND REFERENCE DIRS AND FULLPATH TO NCHG HERE
    """

    input_dir = os.path.abspath("/Users/GBS/Master/Pipeline/diff_res")
    output_dir = os.path.abspath("/Users/GBS/Master/Pipeline/diff_res/output")
    reference_dir = os.path.abspath("/Users/GBS/Master/reference")
    nchg_path = os.path.abspath("/Users/GBS/Master/Pipeline/INC-tutorial/processing_scripts/NCHG_hic/NCHG")

    @classmethod
    def set_input_dir(cls, input_dir):
        cls.input_dir = os.path.abspath(input_dir)

    @classmethod
    def get_input_dir(cls):
        return cls.input_dir

    @classmethod
    def set_output_dir(cls, output_dir):
        cls.output_dir = os.path.abspath(output_dir)

    @classmethod
    def get_output_dir(cls):
        return cls.output_dir

    @classmethod
    def set_reference_dir(cls, reference_dir):
        cls.reference_dir = os.path.abspath(reference_dir)

    @classmethod
    def get_reference_dir(cls):
        return cls.reference_dir

    @classmethod
    def set_NCHG_path(cls, nchg_path):
        cls.NCHG_path = os.path.abspath(nchg_path)

    @classmethod
    def get_NCHG_path(cls):
        return cls.nchg_path

    @classmethod
    def set_temp_dir(cls, temp_dir):
        cls.temp_dir = os.path.abspath(temp_dir)

    @staticmethod
    def get_temp_dir():
        temp_dir = SetDirectories.get_output_dir() + "/temp_dir"
        if not os.path.exists(temp_dir):
            os.makedirs(temp_dir)
        return temp_dir

print(SetDirectories.get_input_dir())
print(type(SetDirectories.get_input_dir()))

class Pipeline_Input:

    @staticmethod
    def find_files(*root_directories):
        """
        Finds all files with the specified suffixes in the specified subdirectory of each root directory.
        Currently, this does not match several files with the same resolition, so only works for one file per resolution.
        :param root_directories: one or more root directories to search in
        :return: a list of file paths for each BED and matrix file found
        """

        subdirectory_name = "raw"
        bedfiles = []
        matrixfiles = []

        # Find the raw data subdirectory in the root directory
        raw_subdirectories = []
        for root_directory in root_directories:
            for root, _, _ in os.walk(root_directory):
                if os.path.basename(root) == subdirectory_name:
                    raw_subdirectories.append(root)

        # Recursively search raw data subdirectory for bed and matrix files
        for subdirectory_path in raw_subdirectories:
            for root, _, files in os.walk(subdirectory_path):
                for file in files:
                    if file.endswith(".bed"):
                        bedfiles.append(os.path.join(root, file))
                    if file.endswith(".matrix"):
                        matrixfiles.append(os.path.join(root, file))


        return bedfiles, matrixfiles

    @staticmethod
    def group_files(*args):
        """
        Groups files by resolution and experiment
        :param args: one or more root directories containing raw data from HiC-Pro
        :return: dict of file paths for each BED and matrix file found, grouped by resolution and experiment
        """

        bedfiles = Pipeline_Input.find_files(*args)[0]
        matrixfiles = Pipeline_Input.find_files(*args)[1]
        grouped_files = {}

        # Extract resolution and experiment name from file path
        for matrixfile in matrixfiles:
            resolution = int(matrixfile.split("/")[-2])
            experiment = matrixfile.split("/")[-4]
            key = f"{experiment, resolution}"

            # Group bed file to matrix file
            for bedfile in bedfiles:
                if bedfile.startswith(matrixfile[:-len(".matrix")]):
                    if key not in grouped_files:
                        grouped_files[key] = (bedfile, matrixfile)
                    else:
                        grouped_files[key] += (bedfile, matrixfile)

        return grouped_files


print(Pipeline_Input.group_files(SetDirectories.get_input_dir()))

class test:
    @staticmethod
    def make_bedpe():
        """
        Makes bedpe file from HiC-Pro output
        """
        bedpe = []

        bed = Pipeline_Input.group_files(SetDirectories.get_input_dir())[0]
        matrix = Pipeline_Input.group_files(SetDirectories.get_input_dir())[1]

        for range in range(len(bed)):
            bed_lines = open(bed, "r").readlines()
            matrix_lines = open(matrix, "r").readlines()

        bed_dict = {}

        for line in bed:
            line = line.strip().split("\t")
            bed_dict[line[3]] = line

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
    def window_size():  # TODO: Make this automatically adjust to resulution of current input file
        # For resoultion in input files, adjust window size to equal resolution for each file
        return 50000

    @staticmethod
    def read_hicpro_output():
        """
        Reads raw HiC-Pro output and output dataframes
        Write dataframes to specified directory
        """

        matrix_file = os.path.join(Pipeline.default_rawdata_path("chr18_50000.matrix"))  # TODO: Make this automatic
        bed_file = os.path.join(Pipeline.default_rawdata_path("chr18_50000_abs.bed"))  # TODO: Make this automatic

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

        bed_file = os.path.join(Pipeline.default_rawdata_path("chr18_50000_abs.bed"))  # TODO: Make this automatic
        matrix_file = os.path.join(Pipeline.default_rawdata_path("chr18_50000.matrix"))  # TODO: Make this automatic

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
    def write_bedpe(*args):  # TODO: Make this automatic
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
    def write_blacklist_hg19(*args):  # TODO: Make this automatic
        Pipeline.remove_blacklist_hg19().saveas(Pipeline.default_output_path(*args))


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
    def write_remove_cytobands(*args):  # TODO: Make this automatic

        with open(Pipeline.default_output_path(*args), "w") as f:
            f.writelines(Pipeline.remove_cytobands().to_dataframe().to_csv(sep="\t", header=False, index=False))

    @staticmethod
    def cap_chromosomes():

        """
        Ensures that chromosome interactions do not go beyond the end of the chromosome
        """

        # read in chromosome sizes
        chromosome_size_hg19 = Pipeline.default_reference_path("hg19/chrom_hg19_test.sizes")
        chromosome_size_dict = {}
        with open(chromosome_size_hg19) as f:
            for line in f:
                chromosome_size_dict[line.split()[0]] = int(line.split()[1])
        interactions = pbt.BedTool.to_dataframe(Pipeline.remove_cytobands()).to_csv(sep="\t", header=False, index=False)
        interactions_in = interactions.split("\n")

        # find interactions exceeding chromosome size
        empty_missing_lines = []
        for line in interactions_in:

            fields = line.strip().split()
            if not line.strip():
                empty_missing_lines.append(line)
                continue
            if len(fields) < 7:
                empty_missing_lines.append(line)
                continue

            chrom1, start1, end1 = fields[0], int(fields[1]), int(fields[2])
            chrom2, start2, end2 = fields[3], int(fields[4]), int(fields[5])

            if end1 > chromosome_size_dict.get(chrom1, end1):
                excess_length = end1 - chromosome_size_dict.get(chrom1, end1)
                end1 = chromosome_size_dict.get(chrom1, end1)
                start1 = max(start1 - excess_length, 0)

            if end2 > chromosome_size_dict.get(chrom2, end2):
                excess_length = end2 - chromosome_size_dict.get(chrom2, end2)
                end2 = chromosome_size_dict.get(chrom2, end2)
                start2 = max(start2 - excess_length, 0)

        # cap interactions exceeding chromosome size
        print(f"{chrom1}:{start1}-{end1} -- {chrom2}:{start2}-{end2}")
        # 78 077 248 chr18
        # 78 000 000-78 050 000 chr18

        # Test this with full genome and write to file/output to pbd for input_to_nchg method (or bypass this method by writing as df.str.csv directly?)

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

        nchg_run = sp.run([Pipeline.path_to_NCHG(), "-m", str(Pipeline.window_size()), "-p", Pipeline.input_to_nchg()], stdout=sp.PIPE, stderr=sp.PIPE, text=True)
        nchg_out = nchg_run.stdout.splitlines()
        return nchg_out

    @staticmethod
    def write_sig_interactions(*args):
        file = os.path.join(Pipeline.default_output_path(*args))
        with open(file, "w") as f:
            f.writelines(Pipeline.find_siginificant_interactions())

    @staticmethod
    def adjust_pvalues(log_ratio_threshold=2, fdr_threshold=0.05, method="fdr_bh"):
        """
        Adjust p-values using the Benjamini-Hochberg method from
        The log-ratio threshold is the observed/expected ratio of interactions
        The fdr threshold is the false discovery rate threshold for the adjusted p-values
        The method is the method used to adjust the p-values (here BH, other methods are available)
        """

        data = Pipeline.find_siginificant_interactions()

        pval = []
        processed = []
        for line in data:
            line = line.split()
            pval.append(float(line[6]))
            if float(line[9]) == 0 or float(line[10]) == 0:
                logratio = 0.0
            else:
                logratio = math.log(float(line[9]), 2) - math.log(float(line[10]), 2)

            line = ' '.join(line) + ' ' + str(logratio)
            processed.append(line)

        padj = list(multicomp.multipletests(pval, method=method))
        padj_out = []

        for i in range(len(processed)):
            line = processed[i] + " " + str(padj[1][i])
            line = line.split()
            if float(line[11]) >= log_ratio_threshold and float(line[12]) <= fdr_threshold:
                padj_out.append(line[0] + "\t" + line[1] + "\t" + line[2] + "\t" + line[3] + "\t" + line[4] + "\t" + line[5] + "\t" + line[12])

        return padj_out


    @staticmethod
    def make_edgelist():

        padj = Pipeline.adjust_pvalues()
        edge_list = []

        for line in padj:
            line = line.split()
            edge_list.append(line[0] + ":" + line[1] + "-" + line[2] + " " + line[3] + "-" + line[4] + ":" + line[5])

        return edge_list

    @staticmethod
    def make_weighted_edgelist():

        padj = Pipeline.adjust_pvalues()
        edge_list = []

        for line in padj:
            line = line.split()
            edge_list.append(line[0] + ":" + line[1] + "-" + line[2] + " " + line[3] + "-" + line[4] + ":" + line[5] + " " + line[6])

        return edge_list

    @staticmethod
    def write_edgelist():

        os.chdir(Pipeline.default_output_path())
        file = Pipeline.default_output_path("edgelist.txt")  # TODO: this file needs to be handled by filehandler for automatic naming (args from filehanlder), name from input file
        with open(file, "w") as f:
            f.writelines(Pipeline.make_edgelist())
        return file

    @staticmethod
    def write_weighted_edgelist():
        os.chdir(Pipeline.default_output_path())
        file_path = Pipeline.default_output_path("weighted_edgelist.txt")
        lines = Pipeline.make_weighted_edgelist()
        FileHandler.write_lines_to_file(file_path, lines)
        return file_path


# Calling Pipeline on data?


# Just call Pipeline above with root dir(s) containing data and output dir
# End pipeline with outputting the edgelist for each experiment/resolution
# and then write all edgelists renamed to shortened name to output directory

# Make main module for running pipeline when its done
# if __name__ == "__main__":
#     main()



