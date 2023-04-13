# Import modules
import concurrent.futures
import shutil
import pybedtools as pbt
import os as os
import subprocess as sp
import math
from statsmodels.sandbox.stats import multicomp
import time as time
from tqdm import tqdm
import re as re
import argparse as argparse
import pickle as pickle
import threading as threading
from collections import defaultdict
import concurrent.futures

#########################################################################################
# Pre-processing pipeline for Hi-C data from HiC-Pro (command line version, parallelized)
#########################################################################################

help_message = "Pipeline for processing Hi-C data from HiC-Pro to statistically significant edge-lists for HG19. Data is run on one output folder from Hi-C pro at a time. \n\n" \
               "INPUT DIR: -i\n" \
               "Directory containing HiC-Pro output folders (bed and matrix files) is set as input directory. Any folder can be the input, as long as it contains the HiC-Pro output folders (raw, matrix) for one HiC-Pro run. \n\n" \
               "OUTPUT DIR: -o \n" \
               "Any directory to output processed data is set as output directory. The output directory will contain extremenly large temp files if running whole genome analysis (At least 60 GB needed). \n\n" \
               "REFERENCE DIR: -r \n" \
               "Directory containing reference genome files is set as reference directory. Put the reference files in a folder called hg19. \n" \
               "This directory should contain the following files: \n" \
               "    cytoBand_hg19.txt (USCS cytoband reference file: https://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/ \n" \
               "    hg19-blacklist.v2.bed (Encode hg19 blacklisted regions: https://github.com/Boyle-Lab/Blacklist/tree/master/lists \n\n" \
               "NCHG PATH: -n \n" \
               "Path to NCHG executable is set as NCHG path. NCHG is a C++ tool for calculating p-values for Hi-C interactions using the NCHG distribution. \n" \
               "It can be found here: https://github.com/Chrom3D/preprocess_scripts/blob/master/NCHG_hic.zip \n\n" \
               "NORM OPTION: -m \n" \
               "Specifies if normalized data or raw data is processed to edge lists. Options: raw, iced, norm, normalized. \n" \
               "If raw is selected, the script will look for raw data in the HiC-Pro output folder. If iced is selected, the script will look for ICE normalized data in the HiC-Pro output folder. \n" \
               "\n\n" \
               "If no arguments are given, the script will run with the hardcoded paths set in the SetDirectories class. Meaning, it's possible to run without providing arguments. \n" \
               "For instance, set the NCHG and reference paths hardcoded, and provide input and output directories for each run. \n\n" \
               "INTER-CHROMOSOMAL INTERACTIONS: -inter. \n" \
               "If specified, run NCHG with the -i flag, meaning interchromosmal interactions will be included in calculating significance. This is set to true by default." \
               "INTRA-CHROMOSOMAL INTERACTIONS: -intra. \n" \
               "If specified, run NCHG without the -i flag, meaning only intrachromosmal interactions will be used to calcualte significance. " \
               "MIXED INTERACTIONS: -mixed \n" \
               "If this flag is set, the NCHG script will consider if inter- or intrachromosomal interactions are appropriate per input file. \n" \
               "Allows for running mixed input files, some of which contain inter-chromosomal interactions and some that do not. This increases time complexity. \n\n" \
               "THEADS: -t \n" \
               "Int: Number of threads to use for processing. Default is cores available on machine. Always specify on HPC cluster. \n\n" \
               "FDR THRESHOLD: -f \n" \
               "Float: FDR threshold for significance. Default is 0.05. \n\n" \


parser = argparse.ArgumentParser(description=help_message, formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument("-i", "--input_dir", help="Directory containing HiC-Pro output folders (bed and matrix files)", required=False)
parser.add_argument("-o", "--output_dir", help="Any directory to output processed data", required=False)
parser.add_argument("-r", "--reference_dir", help="Directory containing reference genome files.", required=False)
parser.add_argument("-np", "--nchg_path", help="Path to NCHG executable", required=False)
parser.add_argument("-n", "--norm_option", help="Normalization option", choices=["raw", "iced", "norm", "normalized"], required=False)
parser.add_argument("-inter", "--interchromosomal_interactions", help="Consider inter-chromosomal interactions to find significant interactions. This is the default setting.", action="store_true", required=False)
parser.add_argument("-intra", "--intrachromosomal_interactions", help="Consider only intra-chromosomal interactions to find significant interactions.", action="store_true", required=False)
parser.add_argument("-mixed", "--mixed_interactions", help="Consider inter-chromosomal and intra-chromosmal interactions on a file-by-file basis.", action="store_true", required=False)
parser.add_argument("-t", "--threads", help="Int: Number of threads to use for processing. Default is cores available on machine. Always specify on HPC cluster.", required=False)
parser.add_argument("-f", "--fdr_threshold", help="Float: FDR threshold for significance. Default is 0.05.", required=False, type=float)
args = parser.parse_args()

# Sets args to None if not provided in command line
input_directory = args.input_dir
output_directory = args.output_dir
reference_directory = args.reference_dir
nchg_executable_path = args.nchg_path
norm_option = args.norm_option
interchromosomal_interactions = args.interchromosomal_interactions
intrachromosomal_interactions = args.intrachromosomal_interactions
mixed_interactions = args.mixed_interactions
fdr_threshold = args.fdr_threshold
threads = args.threads


def check_file_exists(file_path):
    if os.path.exists(file_path):
        print(f"File {file_path} exists")
    else:
        print(f"File {file_path} does not exist")


class SetDirectories:
    """
    SET INPUT-, OUTPUT- AND REFERENCE DIRS AND FULLPATH TO NCHG HERE IF USING HARD CODED PATHS

    For each run, change the input and output directories to the appropriate directories
    Set input dir to root dir containing HiC-Pro output folders (raw, matrix).
    Set normalized data = True to process ICE matrices, False to process raw data.
    """

    input_dir = os.path.abspath("/Users/GBS/Master/HiC-Data/HiC-Pro_out/chr18_inc/chr18_inc")
    output_dir = os.path.abspath("/Users/GBS/Master/HiC-Data/testing_chrom_parallelization/output_terminal")
    reference_dir = os.path.abspath("/Users/GBS/Master/Reference")
    nchg_path = os.path.abspath("/Users/GBS/Master/Scripts/NCHG_hic/NCHG")
    normalized_data = True  # Checks for ICE normalized data in matrix folder
    inter_interactions = False  # If true, considers both intra- and inter-chromosomal interactions for statistical testing
    intra_interactions = False  # If true, considers only intra-chromosomal interactions for statistical testing
    mixed_interactions = True  # If true, considers both inter- and intra-chromosomal interactions for statistical testing on file-by-file basis, using in the NCHG script.
    threads = os.cpu_count()  # Sets threads to number of cores on machine, can be overwritten by user input in command line
    fdr_threshold = 0.05  # Sets FDR threshold for significance, can be overwritten by user input in command line

    @classmethod
    def set_normalized_data(cls, normalized_data):
        cls.normalized_data = normalized_data

    @classmethod
    def get_normalized_data(cls):
        return cls.normalized_data

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
    def set_nchg_path(cls, nchg_path):
        cls.nchg_path = os.path.abspath(nchg_path)

    @classmethod
    def get_nchg_path(cls):
        return cls.nchg_path

    @classmethod
    def set_mixed_interactions(cls, mixed_inter):
        cls.mixed_interactions = mixed_inter

    @classmethod
    def get_mixed_interactions(cls):
        return cls.mixed_interactions

    @classmethod
    def set_interchromosomal_interactions(cls, inter_interactions):
        cls.interchromosomal_interactions = inter_interactions

    @classmethod
    def get_interchromosomal_interactions(cls):
        return cls.inter_interactions

    @classmethod
    def set_intrachromosomal_interactions(cls, intra_interactions):
        cls.intrachromosomal_interactions = intra_interactions

    @classmethod
    def get_intrachromosomal_interactions(cls):
        return cls.intra_interactions

    @classmethod
    def set_fdr_threshold(cls, fdr_thresh):
        cls.fdr_threshold = fdr_thresh

    @classmethod
    def get_fdr_threshold(cls):
        return cls.fdr_threshold

    @classmethod
    def set_temp_dir(cls, temp_dir):
        cls.temp_dir = os.path.abspath(temp_dir)

    @staticmethod
    def get_temp_dir():
        temp_dir = SetDirectories.get_output_dir() + "/temp_dir"
        if not os.path.exists(temp_dir):
            os.makedirs(temp_dir)
        return temp_dir

    @staticmethod
    def set_threads(threads_used):
        SetDirectories.threads = int(threads_used)

    @staticmethod
    def get_threads():
        return SetDirectories.threads

    @staticmethod
    def set_pbt_temp_dir(cls, pbt_temp_dir):
        cls.pbt_temp_dir = os.path.abspath(pbt_temp_dir)

    @staticmethod
    def get_pbt_temp_dir():
        pbt_temp_dir = SetDirectories.get_temp_dir() + "/pbt_temp_dir"
        if not os.path.exists(pbt_temp_dir):
            os.makedirs(pbt_temp_dir)
        return pbt_temp_dir


# Sets input, output and reference directories and flags if to values held in SetDirectories class if not provided by user on command line
if input_directory is not None:
    SetDirectories.set_input_dir(input_directory)

if output_directory is not None:
    SetDirectories.set_output_dir(output_directory)

if reference_directory is not None:
    SetDirectories.set_reference_dir(reference_directory)

if nchg_executable_path is not None:
    SetDirectories.set_nchg_path(nchg_executable_path)

if norm_option is not None:
    if norm_option in ["raw"]:
        SetDirectories.set_normalized_data(False)
    elif norm_option in ["iced", "norm", "normalized"]:
        SetDirectories.set_normalized_data(True)

if mixed_interactions:
    SetDirectories.set_mixed_interactions(mixed_interactions)

if fdr_threshold is not None:
    SetDirectories.set_fdr_threshold(fdr_threshold)

if interchromosomal_interactions is not None:
    SetDirectories.set_interchromosomal_interactions(interchromosomal_interactions)

if intrachromosomal_interactions is not None:
    SetDirectories.set_intrachromosomal_interactions(intrachromosomal_interactions)

if threads is not None:
    SetDirectories.set_threads(threads)

# Sets temporary directory for pbt without cleanup
pbt.set_tempdir(SetDirectories.get_pbt_temp_dir())
pbt.cleanup(False)


class Pipeline_Input:

    @staticmethod
    def find_files(*root_directories):
        """
        Finds all files bed and matrix files in the raw data subdirectory of the root directory.
        :param root_directories: one or more root directories to search in
        :return: a list of file paths for each BED and matrix file found
        """

        raw_subdirectory_name = "raw"
        iced_subdirectory_name = "iced"
        bedfiles = []
        matrixfiles = []
        iced_matrixfiles = []

        # Find the raw data subdirectory in the root directory
        raw_subdirectories = []
        for root_directory in root_directories:
            for root, _, _ in os.walk(root_directory):
                if os.path.basename(root) == raw_subdirectory_name:
                    raw_subdirectories.append(root)

        # Recursively search raw data subdirectory for bed and matrix files
        for subdirectory_path in raw_subdirectories:
            for root, _, files in os.walk(subdirectory_path):
                for file in files:
                    if file.endswith(".bed"):
                        bedfiles.append(os.path.join(root, file))
                    if file.endswith(".matrix"):
                        matrixfiles.append(os.path.join(root, file))

        # Find the ICE-normalized data subdirectory in the root directory
        iced_subdirectories = []
        for root_directory in root_directories:
            for root, _, _ in os.walk(root_directory):
                if os.path.basename(root) == iced_subdirectory_name:
                    iced_subdirectories.append(root)

        # Recursively search ICE-normalized data subdirectory for matrix files
        for subdirectory_path in iced_subdirectories:
            for root, _, files in os.walk(subdirectory_path):
                for file in files:
                    if file.endswith(".matrix"):
                        iced_matrixfiles.append(os.path.join(root, file))

        return bedfiles, matrixfiles, iced_matrixfiles

    @staticmethod
    def group_files(*dirs):
        """
        Groups bed and matrix files by resolution and experiment.
        :param dirs: one or more root directories containing raw data from HiC-Pro
        :return: dict of file paths for each BED and matrix file found, grouped by resolution and experiment
        """

        bedfiles = Pipeline_Input.find_files(*dirs)[0]
        matrixfiles = Pipeline_Input.find_files(*dirs)[1]
        iced_matrixfiles = Pipeline_Input.find_files(*dirs)[2]
        inted_iced_matrixfiles = []

        # Round floats in ICE-normalized matrix files to integers if using ICE normalization
        if SetDirectories.get_normalized_data():
            # Create output directory
            output_dir = os.path.join(SetDirectories.get_temp_dir(), "inted_matrixfiles")
            if not os.path.exists(output_dir):
                os.mkdir(output_dir)
            else:
                shutil.rmtree(output_dir)
                os.mkdir(output_dir)

            # Round floats to integers
            for iced_matrixfile in iced_matrixfiles:
                with open(iced_matrixfile) as f:
                    lines = f.readlines()

                    new_lines = []
                    for line in lines:
                        cols = line.strip().split()
                        cols[2] = str(round(float(cols[2])))
                        new_line = cols[0] + "\t" + cols[1] + "\t" + cols[2] + "\n"
                        new_lines.append(new_line)

                    file_name = os.path.basename(iced_matrixfile)

                    # Create the output file path with the same name as the original file
                    output_file_path = os.path.join(output_dir, file_name)

                    # Save the modified matrix file to the new directory
                    with open(output_file_path, "w") as f_out:
                        f_out.writelines(new_lines)

                    # Add the output file path to the rounded_iced_matrixfiles list
                    inted_iced_matrixfiles.append(output_file_path)

        grouped_raw_files = {}
        grouped_iced_files = {}

        # Extract resolution and experiment name from raw file path
        for matrixfile in matrixfiles:
            resolution = int(matrixfile.split("/")[-2])
            experiment = matrixfile.split("/")[-4]
            key = f"{experiment, resolution}"

            # Group raw bed file to raw matrix file
            for bedfile in bedfiles:
                if bedfile.startswith(matrixfile[:-len(".matrix")]):
                    if key not in grouped_raw_files:
                        grouped_raw_files[key] = (bedfile, matrixfile)
                    else:
                        grouped_raw_files[key] += (bedfile, matrixfile)

        # Extract resolution and experiment name from ICE-normalized file path
        for inted_matrixfile in inted_iced_matrixfiles:
            file_name = os.path.basename(inted_matrixfile)
            experiment, resolution, _ = file_name.rsplit("_", 2)
            resolution = int(resolution)
            key = f"{experiment, resolution}"

            # Group ICE-normalized matrix file
            for bedfile in bedfiles:
                bedfile_name = os.path.basename(bedfile)
                bedfile_experiment, bedfile_resolution, _ = bedfile_name.rsplit("_", 2)
                bedfile_resolution = int(bedfile_resolution)

                if bedfile_experiment == experiment and bedfile_resolution == resolution:
                    if key not in grouped_iced_files:
                        grouped_iced_files[key] = (bedfile, inted_matrixfile)
                    else:
                        grouped_iced_files[key] += (bedfile, inted_matrixfile)

        # Checks if Pipeline should be run on raw or ICE-normalized data
        grouped_files_checked = None
        if SetDirectories.get_normalized_data():
            grouped_files_checked = grouped_iced_files
        if not SetDirectories.get_normalized_data():
            grouped_files_checked = grouped_raw_files

        return grouped_files_checked

first_print = True
def custom_print(*argss, **kwargs):
    global first_print
    if first_print:
        print("\n", end="")
        first_print = False
    print(*argss, **kwargs)

class Pipeline:

    @staticmethod
    def make_bedpe(bed_file, matrix_file):
        """
        Makes bedpe file from HiC-Pro output
        """

        try:
            bedpe = []

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

        except Exception as e:
            tid = threading.get_ident()
            print(f"Error processing file: {bed_file, matrix_file}, {e}, PID: {os.getpid()}, TID: {tid}")
            raise

    @staticmethod
    def input_to_make_bedpe(grouped_files):
        """
        Makes bedpe files from HiC-Pro output and saves them to the temp directory
        :param grouped_files: Output of Pipeline_Input.group_files()
        :return: BEDPE files saved to temp directory
        """

        global first_print
        first_print = True

        os.chdir(SetDirectories.get_temp_dir())
        if not os.path.exists("bedpe"):
            os.mkdir("bedpe")
        else:
            shutil.rmtree("bedpe")
            os.mkdir("bedpe")
        os.chdir("bedpe")

        for key, val in grouped_files.items():

            split_key = key.split(",")
            experiment = split_key[0].replace("'", "").replace("(", "").replace(")", "").replace(" ", "_")
            resolution = split_key[1].replace("'", "").replace("(", "").replace(")", "").replace(" ", "_")
            # Remove trailing underscore or space
            if experiment[-1] in ("_", " "):
                experiment = experiment[:-1]
            if resolution[-1] in ("_", " "):
                resolution = resolution[:-1]
            # Replace consecutive underscores with a single underscore
            experiment = '_'.join(filter(None, experiment.split('_')))
            resolution = '_'.join(filter(None, resolution.split('_')))

            # Should separate the string formatting from the file handling?
            bedfile = val[0]
            matrixfile = val[1]
            bedpe = Pipeline.make_bedpe(bedfile, matrixfile)

            with open(experiment + "_" + resolution + ".bedpe", "w") as f:
                f.writelines(bedpe)
                f.close()

    @staticmethod
    def remove_blacklisted_regions(bedpe_file):
        """
        Removes blacklisted regions from bedpe file
        """
        try:
            reference_dir = SetDirectories.get_reference_dir()
            blacklisted_regions_path = os.path.join(reference_dir, "hg19/hg19-blacklist.v2.bed")
            with open(blacklisted_regions_path, "r") as f:
                blacklisted_regions = f.readlines()

            # Filter out mitochondrial DNA and Y chromosome
            with open(bedpe_file, "r") as f:
                bedpe_data = f.readlines()
            filtered_bedpe_data = [line for line in bedpe_data if not any(chrom in line for chrom in ["chrM", "chrY"])]

            blacklisted_pbt = pbt.BedTool(blacklisted_regions)
            blacklisted_bedpe = pbt.BedTool(filtered_bedpe_data)

            window_size = int(re.search(r"(\d+)[^/\d]*$", bedpe_file).group(1))
            if not isinstance(window_size, int):
                raise ValueError(f"Window size must be an integer, {window_size} is not an integer.")

            no_overlap_bedpe = blacklisted_bedpe.window(blacklisted_pbt, w=int(window_size), r=False, v=True)
            custom_print(f"Finished processing file: {bedpe_file}, PID: {os.getpid()}, TID: {threading.get_ident()}")

            return no_overlap_bedpe

        except Exception as e:
            tid = threading.get_ident()
            print(f"\nError processing file: {bedpe_file}, {e}, PID: {os.getpid()}, TID: {tid}")
            raise

    @staticmethod
    def input_to_remove_blacklist():
        """
        Calls the remove_blacklisted_regions function on each file in the BEDPE directory
        """

        global first_print
        first_print = True

        bedpe_dir = SetDirectories.get_temp_dir() + "/bedpe"
        bedpe_files = [os.path.join(bedpe_dir, file) for file in os.listdir(bedpe_dir)]

        # Create the output directory if it doesn't exist
        output_dir = os.path.join(SetDirectories.get_temp_dir(), "blacklisted")
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
        else:
            shutil.rmtree(output_dir)
            os.mkdir(output_dir)

        try:
            pickle.dumps(bedpe_files)
        except Exception as e:
            print(f"Pickling error: {e}")

        for file in bedpe_files:
            if not os.path.isfile(file):
                print(f"File {file} does not exist.")

        with concurrent.futures.ThreadPoolExecutor(max_workers=SetDirectories.get_threads()) as executor:
            futures = list(executor.map(Pipeline.remove_blacklisted_regions, bedpe_files))
            for bedpe_file, future in zip(bedpe_files, futures):
                try:
                    no_blacklist_bedpe = future
                    output_filename = os.path.basename(bedpe_file)[:-len('.bedpe')] + '_no_blacklist.bedpe'
                    converted_nb_bedpe = no_blacklist_bedpe.to_dataframe().to_csv(sep="\t", index=False, header=False)
                    output_filepath = os.path.join(output_dir, output_filename)
                    with open(output_filepath, "w") as f:
                        f.writelines(converted_nb_bedpe)

                except Exception as e:
                    tid = threading.get_ident()
                    print(f"Error processing {bedpe_file}: {e}, PID: {os.getpid()}, TID: {tid}")

    @staticmethod
    def remove_cytobands(blacklisted_bedpe_file):
        """
        Cytoband locations are determined in this case by Giemsa staining
        and are located and removed from the BEDPE file.
        """

        try:
            os.chdir(SetDirectories.get_reference_dir())
            cytobands = os.path.join(SetDirectories.get_reference_dir(), "hg19/cytoBand_hg19.txt")
            centromeric_regions = []
            with open(cytobands, "r") as f:
                cytobands = f.readlines()
                for line in cytobands:
                    line = line.strip().split("\t")
                    if line[4] == "acen":
                        centromeric_regions.append(line[0:5])

            os.chdir(SetDirectories.get_temp_dir() + "/blacklisted")
            blacklisted_pbt = pbt.BedTool(blacklisted_bedpe_file)

            window_size = int(re.search(r"(\d+)[^/\d]*$", blacklisted_bedpe_file).group(1))
            if not isinstance(window_size, int):
                raise ValueError(f"Window size must be an integer, {window_size} is not an integer.")

            no_cytobands = blacklisted_pbt.window(centromeric_regions, w=int(window_size), r=False, v=True)
            custom_print(f"Finished processing file: {blacklisted_bedpe_file}, PID: {os.getpid()}, TID: {threading.get_ident()}")
            return no_cytobands

        except Exception as e:
            tid = threading.get_ident()
            print(f"Error processing {blacklisted_bedpe_file}: {e}, PID: {os.getpid()}, TID: {tid}")
            raise

    @staticmethod
    def input_to_remove_cytobands():
        """
        Calls the remove_cytobands function on each file in the blacklisted directory
        """

        global first_print
        first_print = True

        blacklisted_dir_path = SetDirectories.get_temp_dir() + "/blacklisted"
        blacklisted_dir = os.listdir(blacklisted_dir_path)

        # Create the output directory if it doesn't exist
        output_dir = os.path.join(SetDirectories.get_temp_dir(), "no_cytobands")
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
        else:
            shutil.rmtree(output_dir)
            os.mkdir(output_dir)

        for file in blacklisted_dir:
            full_path = os.path.join(blacklisted_dir_path, file)
            if not os.path.isfile(full_path):
                print(f"File {file} does not exist.")

        with concurrent.futures.ThreadPoolExecutor(max_workers=SetDirectories.get_threads()) as executor:
            full_paths = [os.path.join(blacklisted_dir_path, file) for file in blacklisted_dir]
            futures = list(executor.map(Pipeline.remove_cytobands, full_paths))
            for bedpe_file, future in zip(blacklisted_dir, futures):
                try:
                    no_cytoband_bedpe = future
                    output_filename = f"{bedpe_file[:-len('.bedpe')]}_no_cytobands.bedpe"
                    converted_nc_bedpe = no_cytoband_bedpe.to_dataframe().to_csv(sep="\t", index=False, header=False)
                    output_filepath = os.path.join(output_dir, output_filename)
                    with open(output_filepath, "w") as f:
                        f.writelines(converted_nc_bedpe)

                except Exception as e:
                    tid = threading.get_ident()
                    print(f"Error processing {bedpe_file}: {e}, PID: {os.getpid()}, TID: {tid}")

    @staticmethod
    def find_significant_interactions(bedpe_file):
        """
        NCHG script to calculate the significance of interactions:
        m = minimum interaction length in bp, should be same as window size used to make the bedpe file (resolution)
        p = print counts, needed for FDR correction for padj with BH method.
        i = Use interchromosomal interactions as well when calculating p-values (w in args for this script)
        """

        try:
            # Setting min interactions length to same as bin size from HiC-Pro
            window_size = int(re.search(r"_(\d+)_", bedpe_file).group(1))  # Any number in file name with underscores on both sides
            if not isinstance(window_size, int):
                raise ValueError(f"Window size must be an integer, {window_size} is not an integer.")

            # Run NCHG
            nchg_flags = []
            if SetDirectories.inter_interactions:
                nchg_flags.append("-i")
            elif SetDirectories.mixed_interactions:
                intrachromosomal = True
                with open(bedpe_file, "r") as f:
                    for line in f:
                        columns = line.strip().split("\t")
                        if columns[0] != columns[3]:
                            intrachromosomal = False
                            break

                if not intrachromosomal:
                    nchg_flags.append("-i")

            if not SetDirectories.inter_interactions and not SetDirectories.intra_interactions and not SetDirectories.mixed_interactions:
                raise ValueError("Must select at least one type of interaction type for significance testing: Inter, intra, or mixed.")

            nchg_command = [SetDirectories.get_nchg_path(), bedpe_file, "-m", str(window_size), "-p"] + nchg_flags
            nchg_run = sp.run(nchg_command, capture_output=True, check=True)
            print(f"Finished processing file: {bedpe_file}, PID: {os.getpid()}, TID: {threading.get_ident()}")
            return nchg_run.stdout.decode("utf-8").split("\t")

        except Exception as e:
            tid = threading.get_ident()
            print(f"Error in {bedpe_file}: {e}, PID: {os.getpid()}, TID: {tid}")
            raise

    @staticmethod
    def split_bedpe_by_chromosome(bedpe_file, output_dir):
        try:
            with open(bedpe_file, "r") as f:
                data = f.readlines()

            # Group data by chromosome
            chromosomes = {}
            for line in data:
                chr_name = line.split("\t")[0]
                if chr_name not in chromosomes:
                    chromosomes[chr_name] = []
                chromosomes[chr_name].append(line.strip())  # Strip whitespace from the line

            # Write each chromosome to a separate file
            chr_files = []
            for chr_name, chr_data in chromosomes.items():
                output_filename = f"{os.path.basename(bedpe_file)[:-len('.bedpe')]}_{chr_name}_split.bedpe"
                output_filepath = os.path.join(output_dir, output_filename)
                with open(output_filepath, "w") as f:
                    f.write("\n".join(chr_data) + "\n")
                chr_files.append(output_filepath)

            return chr_files

        except Exception as e:
            tid = threading.get_ident()
            print(f"Error in {bedpe_file}: {e}, PID: {os.getpid()}, TID: {tid}")
            raise

    @staticmethod
    def input_to_nchg():
        """
        Calls the NCHG script on all files in the no_cytobands directory to find significant interactions
        """

        no_cytobands_dir_path = SetDirectories.get_temp_dir() + "/no_cytobands"
        no_cytobands_dir = os.listdir(no_cytobands_dir_path)

        # Create the output directory if it doesn't exist
        output_dir = os.path.join(SetDirectories.get_temp_dir(), "NCHG_output")
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
        else:
            shutil.rmtree(output_dir)
            os.mkdir(output_dir)

        # Create a directory to store split chromosome files
        chr_split_base_dir = os.path.join(SetDirectories.get_temp_dir(), "input_to_nchg")
        if not os.path.exists(chr_split_base_dir):
            os.mkdir(chr_split_base_dir)
        else:
            shutil.rmtree(chr_split_base_dir)
            os.mkdir(chr_split_base_dir)

        # Split input files by chromosome
        all_chr_files = []
        # keep original input file name for each chromosome file
        input_file_map = {}

        for file in no_cytobands_dir:
            full_path = os.path.join(no_cytobands_dir_path, file)
            if not os.path.isfile(full_path):
                print(f"File {file} does not exist.")

            # Create a subdirectory for each input file inside the chr_split_base_dir
            file_split_dir = os.path.join(chr_split_base_dir, f"{file[:-len('_no_blacklist_no_cytobands.bedpe')]}_split")
            if not os.path.exists(file_split_dir):
                os.mkdir(file_split_dir)
            else:
                shutil.rmtree(file_split_dir)
                os.mkdir(file_split_dir)

            chr_files = Pipeline.split_bedpe_by_chromosome(full_path, file_split_dir)

            for chr_file in chr_files:
                input_file_map[chr_file] = file
            all_chr_files.extend(chr_files)

        # Run find_significant_interactions on chromosome-specific files in parallel
        output_file_data = defaultdict(list)
        with concurrent.futures.ProcessPoolExecutor(max_workers=SetDirectories.get_threads()) as executor:
            futures = list(executor.map(Pipeline.find_significant_interactions, all_chr_files))
            for bedpe_file, future in zip(all_chr_files, futures):
                try:
                    nchg_output = future
                    input_file = input_file_map[bedpe_file]
                    output_file_data[input_file].append(nchg_output)
                except Exception as e:
                    tid = threading.get_ident()
                    print(f"Error processing {bedpe_file}: {e}, PID: {os.getpid()}, TID: {tid}")

        # Merge output files back together
        for input_file, nchg_outputs in output_file_data.items():
            output_filename = f"{os.path.basename(input_file)[:-len('_no_blacklist_no_cytobands.bedpe')]}_nchg_output.txt"
            output_filepath = os.path.join(output_dir, output_filename)
            with open(output_filepath, "w") as f:
                for nchg_output in nchg_outputs:
                    f.writelines(nchg_output)

    @staticmethod
    def adjust_pvalues(nchg_file, fdr_thresh=SetDirectories.get_fdr_threshold(), log_ratio_threshold=2, method="fdr_bh"):
        """
        Adjusts the p-values using the Benjamini-Hochberg method
        """

        try:
            # Finds the p-values and log ratios of the interactions
            p_values = []
            processed_lines = []
            with open(nchg_file, "r") as nchg_f:
                for line in nchg_f:
                    col = line.split()
                    if len(col) < 11:  # Expect 11 columns in the input file
                        print(f"Error: {nchg_file} does not have the correct number of columns.")
                    if col[7] == "0":
                        continue  # Skips the line if interactions/edges is 0

                    p_values.append(float(col[6]))

                    log_ratio = 0.0
                    if float(col[9]) != 0 and float(col[10]) != 0:
                        log_ratio = math.log(float(col[9]), 2) - math.log(float(col[10]), 2)

                    processed_line = ' '.join(col) + ' ' + str(log_ratio)
                    processed_lines.append(processed_line)

            padj = list(multicomp.multipletests(p_values, method=method))
            padj_out = []

            # Filters the interactions based on the log ratio and FDR thresholds
            for i, processed_line in enumerate(processed_lines):
                col = processed_line.split()
                col.append(str(padj[1][i]))
                if float(col[11]) >= log_ratio_threshold and float(col[12]) <= fdr_thresh:
                    padj_out.append("\t".join(col[:6] + [col[12]]))

            print(f"Finished processing file: {nchg_file}, PID: {os.getpid()}, TID: {threading.get_ident()}")
            return padj_out

        except Exception as e:
            tid = threading.get_ident()
            print(f"Error in {nchg_file}: {e}, PID: {os.getpid()}, TID: {tid}")
            raise

    @staticmethod
    def input_to_adjust_pvalues():
        """
        Calls the adjust_pvalues function on all files in the NCHG_output directory
        """

        nchg_dir_path = SetDirectories.get_temp_dir() + "/NCHG_output"
        nchg_dir = os.listdir(nchg_dir_path)

        # Create the output directory if it doesn't exist
        output_dir = os.path.join(SetDirectories.get_temp_dir(), "padj")
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
        else:
            shutil.rmtree(output_dir)
            os.mkdir(output_dir)

        for file in nchg_dir:
            full_path = os.path.join(nchg_dir_path, file)
            if not os.path.isfile(full_path):
                print(f"File {file} does not exist.")
                continue

        with concurrent.futures.ProcessPoolExecutor(max_workers=SetDirectories.get_threads()) as executor:
            full_paths = [os.path.join(nchg_dir_path, file) for file in nchg_dir]
            futures = list(executor.map(Pipeline.adjust_pvalues, full_paths))
            for nchg_file, future in zip(nchg_dir, futures):
                try:
                    padj = future
                    output_filename = f"{nchg_file[:-len('_nchg_output.txt')]}_padj.txt"
                    output_filepath = os.path.join(output_dir, output_filename)
                    with open(output_filepath, "w") as f:
                        for line in padj:
                            f.write(line + "\n")
                except Exception as e:
                    tid = threading.get_ident()
                    print(f"Error processing {nchg_file}: {e}, PID: {os.getpid()}, TID: {tid}")

    @staticmethod
    def make_weighted_edgelist(padj_file_path):
        """
        makes a weighted edgelist from padj file, padj values are weights
        """

        try:
            edge_list = []
            with open(padj_file_path, "r") as padj_file:
                for line in padj_file:
                    line = line.split()
                    edge_list.append(line[0] + ":" + line[1] + "-" + line[2] + " " + line[3] + ":" + line[4] + "-" + line[5] + " " + line[6])

            print(f"Finished processing file: {padj_file}, PID: {os.getpid()}, TID: {threading.get_ident()}")
            return edge_list

        except Exception as e:
            tid = threading.get_ident()
            print(f"Error in {padj_file}: {e}, PID: {os.getpid()}, TID: {tid}")
            raise

    @staticmethod
    def input_to_make_weighted_edgelist():
        """
        calls make_weighted_edgelist on all padj files
        """

        global first_print
        first_print = True

        padj_dir_path = SetDirectories.get_temp_dir() + "/padj"
        padj_dir = os.listdir(padj_dir_path)

        # Create the output directory if it doesn't exist
        output_dir = os.path.join(SetDirectories.get_temp_dir(), "weighted_edgelists")
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
        else:
            shutil.rmtree(output_dir)
            os.mkdir(output_dir)

        for file in padj_dir:
            full_path = os.path.join(padj_dir_path, file)
            if not os.path.isfile(full_path):
                print(f"File {file} does not exist.")

        with concurrent.futures.ThreadPoolExecutor(max_workers=SetDirectories.get_threads()) as executor:
            full_paths = [os.path.join(padj_dir_path, file) for file in padj_dir]
            futures = list(executor.map(Pipeline.make_weighted_edgelist, full_paths))
            for padj_file_path, future in zip(padj_dir, futures):
                try:
                    weighted_edgelist = future
                    output_filename = f"{os.path.basename(padj_file_path)[:-len('_padj.txt')]}_edgelist.txt"
                    output_filepath = os.path.join(output_dir, output_filename)
                    with open(output_filepath, "w") as f:
                        for line in weighted_edgelist:
                            f.write(line + "\n")
                except Exception as e:
                    tid = threading.get_ident()
                    print(f"Error processing {padj_file_path}: {e}, PID: {os.getpid()}, TID: {tid}")

    @staticmethod
    def make_edgelist(padj_file_path):

        """Makes edge list from padj file"""

        try:
            edge_list = []
            with open(padj_file_path, "r") as padj_file:
                for line in padj_file:
                    line = line.split()
                    edge_list.append(line[0] + ":" + line[1] + "-" + line[2] + "  " + line[3] + ":" + line[4] + "-" + line[5])

            custom_print(f"Finished processing file: {padj_file}, PID: {os.getpid()}, TID: {threading.get_ident()}")
            return edge_list

        except Exception as e:
            tid = threading.get_ident()
            print(f"Error in {padj_file}: {e}, PID: {os.getpid()}, tid: {tid}")
            raise

    @staticmethod
    def input_to_make_edgelist():
        """
        Calls make_edgelist on all padj files
        """

        global first_print
        first_print = True

        padj_dir_path = SetDirectories.get_temp_dir() + "/padj"
        padj_dir = os.listdir(padj_dir_path)

        # Create the output directory if it doesn't exist
        output_dir = os.path.join(SetDirectories.get_output_dir(), "edgelists")
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
        else:
            shutil.rmtree(output_dir)
            os.mkdir(output_dir)

        full_paths = [os.path.join(padj_dir_path, file) for file in padj_dir]

        with concurrent.futures.ThreadPoolExecutor(max_workers=SetDirectories.get_threads()) as executor:
            futures = list(executor.map(Pipeline.make_edgelist, full_paths))
            for padj_file_path, future in zip(full_paths, futures):
                try:
                    edgelist = future
                    output_filename = f"{os.path.basename(padj_file_path)[:-len('_padj.txt')]}_edgelist.txt"
                    output_filepath = os.path.join(output_dir, output_filename)
                    with open(output_filepath, "w") as f:
                        for line in edgelist:
                            f.write(line + "\n")
                except Exception as e:
                    tid = threading.get_ident()
                    print(f"Error processing {padj_file_path}: {e}, PID: {os.getpid()}, TID: {tid}")


def run_pipeline():
    """
    Call selected methods of the Pipeline, in the order specified
    """

    global first_print
    first_print = True

    start_time = time.time()

    # List of static method names to call
    method_names = [
        (lambda: Pipeline.input_to_make_bedpe(Pipeline_Input.group_files(SetDirectories.get_input_dir()))),
        "input_to_remove_blacklist",
        "input_to_remove_cytobands",
        "input_to_nchg",
        "input_to_adjust_pvalues",
        "input_to_make_edgelist",
        "input_to_make_weighted_edgelist"
    ]

    # Call each method once
    for method in tqdm(method_names):
        if type(method) == str:  # Check if method name is a string
            method = getattr(Pipeline, method)  # Get the method reference
        method()  # Call the method

    # Print runtime on completion
    end_time = time.time()
    print(f"Pipeline completed in {end_time - start_time:.2f} seconds. ({(end_time - start_time) / 60:.2f} minutes). ")


if __name__ == "__main__":
    run_pipeline()