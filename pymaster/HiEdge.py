# Import modules
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


###########################################################################
# Pre-processing pipeline for Hi-C data from HiC-Pro (command line version)
###########################################################################

help_message = "Pipeline for processing Hi-C data from HiC-Pro to statistically significant edge-lists for HG19. Data is run on one output folder from Hi-C pro at a time. \n\n" \
               "INPUT DIR:\n" \
               "Directory containing HiC-Pro output folders (bed and matrix files) is set as input directory. Any folder can be the input, as long as it contains the HiC-Pro output folders (raw, matrix) for one HiC-Pro run. \n\n" \
               "OUTPUT DIR: \n"\
               "Any directory to output processed data is set as output directory. The output directory will contain extremenly large temp files if running whole genome analysis (At least 60 GB needed). \n\n" \
               "REFERENCE DIR:\n" \
               "Directory containing reference genome files is set as reference directory. Put the reference files in a folder called hg19. \n" \
               "This directory should contain the following files: \n" \
               "    cytoBand_hg19.txt (USCS cytoband reference file: https://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/ \n" \
               "    hg19-blacklist.v2.bed (Encode hg19 blacklisted regions: https://github.com/Boyle-Lab/Blacklist/tree/master/lists \n\n"\
               "NCHG PATH: \n"\
               "Path to NCHG executable is set as NCHG path. NCHG is a C++ tool for calculating p-values for Hi-C interactions using the NCHG distribution. \n" \
               "It can be found here: https://github.com/Chrom3D/preprocess_scripts/blob/master/NCHG_hic.zip \n\n"\
               "NORM OPTION: \n"\
               "Specifies if normalized data or raw data is processed to edge lists. Options: raw, iced, norm, normalized. \n"\
               "If raw is selected, the script will look for raw data in the HiC-Pro output folder. If iced is selected, the script will look for ICE normalized data in the HiC-Pro output folder. \n"\
               "\n\n"\
               "If no arguments are given, the script will run with the hardcoded paths set in the SetDirectories class. Meaning, it's possible to run without providing arguments. \n"\
               "For instance, set the NCHG and reference paths hardcoded, and provide input and output directories for each run. \n"

parser = argparse.ArgumentParser(description=help_message, formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument("-i", "--input_dir", help="Directory containing HiC-Pro output folders (bed and matrix files)", required=False)
parser.add_argument("-o", "--output_dir", help="Any directory to output processed data", required=False)
parser.add_argument("-r", "--reference_dir", help="Directory containing reference genome files.", required=False)
parser.add_argument("-n", "--nchg_path", help="Path to NCHG executable", required=False)
parser.add_argument("-m", "--norm_option", help="Normalization option", choices=["raw", "iced", "norm", "normalized"], required=False)

args = parser.parse_args()

input_directory = args.input_dir if args.input_dir is not None else None
output_directory = args.output_dir if args.output_dir is not None else None
reference_directory = args.reference_dir if args.reference_dir is not None else None
nchg_executable_path = args.nchg_path if args.nchg_path is not None else None
norm_option = args.norm_option if args.norm_option is not None else None


class SetDirectories:
    """
    SET INPUT-, OUTPUT- AND REFERENCE DIRS AND FULLPATH TO NCHG HERE IF USING HARD CODED PATHS

    For each run, change the input and output directories to the appropriate directories
    Set input dir to root dir containing HiC-Pro output folders (raw, matrix).
    Set normalized data = True to process ICE matrices, False to process raw data.
    """

    input_dir = os.path.abspath("")
    output_dir = os.path.abspath("")
    reference_dir = os.path.abspath("")
    nchg_path = os.path.abspath("")
    normalized_data = True  # Checks for ICE normalized data in matrix folder

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

if input_directory is not None:
    SetDirectories.set_input_dir(input_directory)

if output_directory is not None:
    SetDirectories.set_output_dir(output_directory)

if reference_directory is not None:
    SetDirectories.set_reference_dir(reference_directory)

if nchg_executable_path is not None:
    SetDirectories.set_NCHG_path(nchg_executable_path)

if norm_option is not None:
    if norm_option in ["raw"]:
        SetDirectories.set_normalized_data(False)
    elif norm_option in ["iced", "norm", "normalized"]:
        SetDirectories.set_normalized_data(True)

# Set PyBedtools tmp dir to temp dir
pbt.set_tempdir(SetDirectories.get_temp_dir())

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
    def group_files(*args):
        """
        Groups bed and matrix files by resolution and experiment.
        :param args: one or more root directories containing raw data from HiC-Pro
        :return: dict of file paths for each BED and matrix file found, grouped by resolution and experiment
        """

        bedfiles = Pipeline_Input.find_files(*args)[0]
        matrixfiles = Pipeline_Input.find_files(*args)[1]
        iced_matrixfiles = Pipeline_Input.find_files(*args)[2]
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


class Pipeline:
    @staticmethod
    def make_bedpe(bed_file, matrix_file):
        """
        Makes bedpe file from HiC-Pro output
        """

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

    @staticmethod
    def input_to_make_bedpe(grouped_files):
        """
        Makes bedpe files from HiC-Pro output and saves them to the temp directory
        :param grouped_files: Output of Pipeline_Input.group_files()
        :return: BEDPE files saved to temp directory
        """

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

        os.chdir(SetDirectories.get_reference_dir())
        blacklisted_regions = open("hg19/hg19-blacklist.v2.bed", "r").readlines()
        blacklisted_pbt = pbt.BedTool(blacklisted_regions)

        os.chdir(SetDirectories.get_temp_dir() + "/bedpe")
        blacklisted_bedpe = pbt.BedTool(bedpe_file)

        window_size = int(re.search(r"(\d+)[^/\d]*$", bedpe_file).group(1))
        if not isinstance(window_size, int):
            raise ValueError(f"Window size must be an integer, {window_size} is not an integer.")

        no_overlap_bedpe = blacklisted_bedpe.window(blacklisted_pbt, w=int(window_size), r=False, v=True)

        return no_overlap_bedpe

    @staticmethod
    def input_to_remove_blacklist():
        """
        Calls remove_blacklisted_regions() on all bedpe files in the bedpe directory
        """
        bedpe_dir = os.listdir(SetDirectories.get_temp_dir() + "/bedpe")

        # Create the output directory if it doesn't exist
        output_dir = os.path.join(SetDirectories.get_temp_dir(), "blacklisted")
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
        else:
            shutil.rmtree(output_dir)
            os.mkdir(output_dir)

        # Iterate over input files and process/save them
        for bedpe_file in bedpe_dir:
            os.chdir(SetDirectories.get_temp_dir() + "/bedpe")
            no_blacklist_bedpe = Pipeline.remove_blacklisted_regions(bedpe_file)

            os.chdir(output_dir)
            output_filename = f"{bedpe_file[:-len('.bedpe')]}_no_blacklist.bedpe"
            converted_nbl_bedpe = no_blacklist_bedpe.to_dataframe().to_csv(sep="\t", index=False, header=False)
            with open(os.path.join(output_dir, output_filename), "w") as f:
                f.writelines(converted_nbl_bedpe)

    @staticmethod
    def remove_cytobands(blacklisted_bedpe_file):
        """
        Cytoband locations are determined in this case by Giemsa staining
        and are located and removed from the BEDPE file.
        """

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

        return no_cytobands

    @staticmethod
    def input_to_remove_cytobands():
        """
        Calls the remove_cytobands function on each file in the blacklisted directory
        """

        # Create the output directory if it doesn't exist
        output_dir = os.path.join(SetDirectories.get_temp_dir(), "no_cytobands")
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
        else:
            shutil.rmtree(output_dir)
            os.mkdir(output_dir)

        # Iterate over input files and process/save them
        blacklisted_dir = os.listdir(SetDirectories.get_temp_dir() + "/blacklisted")
        for blacklisted_file in blacklisted_dir:
            os.chdir(SetDirectories.get_temp_dir() + "/blacklisted")
            no_cytoband_bedpe = Pipeline.remove_cytobands(blacklisted_file)

            os.chdir(output_dir)
            output_filename = f"{blacklisted_file[:-len('.bedpe')]}_no_cytobands.bedpe"
            converted_nc_bedpe = no_cytoband_bedpe.to_dataframe().to_csv(sep="\t", index=False, header=False)
            with open(os.path.join(output_dir, output_filename), "w") as f:
                f.writelines(converted_nc_bedpe)

    @staticmethod
    def find_siginificant_interactions(bedpe_file):
        """
        NCHG script to calculate the significance of interactions:
        m = minimum interaction length in bp, should be same as window size used to make the bedpe file (resolution)
        p = input file, which is the output of the remove_cytobands function but reformatted to be compatible with NCHG
        """

        # Setting window size to be the same as the resolution
        window_size = int(re.search(r"(\d+)[^/\d]*$", bedpe_file).group(1))
        if not isinstance(window_size, int):
            raise ValueError(f"Window size must be an integer, {window_size} is not an integer.")

        # Run NCHG
        nchg_run = sp.run([SetDirectories.get_NCHG_path(), "-m", str(window_size), "-p", bedpe_file], capture_output=True)

        return nchg_run.stdout.decode("utf-8").split("\t")

    @staticmethod
    def input_to_nchg():
        """
        Calls the NCHG script on all files in the no_cytobands directory to find significant interactions
        """

        no_cytobands_dir = os.listdir(SetDirectories.get_temp_dir() + "/no_cytobands")

        # Create the output directory if it doesn't exist
        output_dir = os.path.join(SetDirectories.get_temp_dir(), "NCHG_output")
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
        else:
            shutil.rmtree(output_dir)
            os.mkdir(output_dir)

        # Iterate over input files and process/save them
        for no_cytobands_file in no_cytobands_dir:

            # run the NCHG script on the input bedpe file
            os.chdir(SetDirectories.get_temp_dir() + "/no_cytobands")
            nchg_output = Pipeline.find_siginificant_interactions(no_cytobands_file)

            # save the NCHG output to a file
            os.chdir(output_dir)
            output_filename = f"{no_cytobands_file[:-len('.bedpe')]}_nchg_output.txt"
            with open(os.path.join(output_dir, output_filename), "w") as f:
                f.writelines(nchg_output)

    @staticmethod
    def adjust_pvalues(nchg_file, fdr_threshold=0.01, log_ratio_threshold=2, method="fdr_bh"):
        """
        Adjusts the p-values using the Benjamini-Hochberg method
        """

        # Finds the p-values and log ratios of the interactions
        pval = []
        processed = []
        with open(nchg_file, "r") as nchg_file:
            for line in nchg_file:
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

        # Filters the interactions based on the log ratio and FDR thresholds
        for i in range(len(processed)):
            line = processed[i] + " " + str(padj[1][i])
            line = line.split()
            if float(line[11]) >= log_ratio_threshold and float(line[12]) <= fdr_threshold:
                padj_out.append(line[0] + "\t" + line[1] + "\t" + line[2] + "\t" + line[3] + "\t" + line[4] + "\t" + line[5] + "\t" + line[12])

        return padj_out

    @staticmethod
    def input_to_adjust_pvalues():
        """
        calls the adjust_pvalues function on all files in the NCHG_output directory
        """

        nchg_dir = os.listdir(SetDirectories.get_temp_dir() + "/NCHG_output")
        # Create the output directory if it doesn't exist
        output_dir = os.path.join(SetDirectories.get_temp_dir(), "padj")
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
        else:
            shutil.rmtree(output_dir)
            os.mkdir(output_dir)

        for nchg_file in nchg_dir:
            os.chdir(SetDirectories.get_temp_dir() + "/NCHG_output")
            padj = Pipeline.adjust_pvalues(nchg_file)

            os.chdir(output_dir)
            output_filename = f"{nchg_file[:-len('_nchg_output.txt')]}_padj.txt"
            with open(os.path.join(output_dir, output_filename), "w") as f:
                for line in padj:
                    f.write(line + "\n")

    @staticmethod
    def make_weighted_edgelist(padj_file):
        """
        makes a weighted edgelist from padj file, padj values are weights
        """

        edge_list = []
        with open(padj_file, "r") as padj_file:
            for line in padj_file:
                line = line.split()
                edge_list.append(line[0] + ":" + line[1] + "-" + line[2] + " " + line[3] + "-" + line[4] + ":" + line[5] + " " + line[6])

        return edge_list

    @staticmethod
    def input_to_make_weighted_edgelist():
        """
        calls make_weighted_edgelist on all padj files
        """

        os.chdir(SetDirectories.get_temp_dir() + "/padj")
        padj_dir = os.listdir(SetDirectories.get_temp_dir() + "/padj")
        # Create the output directory if it doesn't exist
        output_dir = os.path.join(SetDirectories.get_temp_dir(), "weighted_edgelists")
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
        else:
            shutil.rmtree(output_dir)
            os.mkdir(output_dir)

        for padj_file in padj_dir:
            os.chdir(SetDirectories.get_temp_dir() + "/padj")
            weighted_edgelist = Pipeline.make_weighted_edgelist(padj_file)

            os.chdir(output_dir)
            output_filename = f"{padj_file[:-len('_no_blacklist_no_cytobands_padj.txt')]}_weighted_edgelist.txt"
            with open(os.path.join(output_dir, output_filename), "w") as f:
                for line in weighted_edgelist:
                    f.write(line + "\n")

    @staticmethod
    def make_edgelist(padj_file):
        """Makes edge list from padj file"""

        edge_list = []
        with open(padj_file, "r") as padj_file:
            for line in padj_file:
                line = line.split()
                edge_list.append(line[0] + ":" + line[1] + "-" + line[2] + "  " + line[3] + "-" + line[4] + ":" + line[5])

        return edge_list

    @staticmethod
    def input_to_make_edgelist():
        """
        Calls make_edgelist on all padj files
        """

        os.chdir(SetDirectories.get_temp_dir() + "/padj")
        padj_dir = os.listdir(SetDirectories.get_temp_dir() + "/padj")
        # Create the output directory if it doesn't exist
        output_dir = os.path.join(SetDirectories.get_output_dir(), "edgelists")
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
        else:
            shutil.rmtree(output_dir)
            os.mkdir(output_dir)

        for padj_file in padj_dir:
            os.chdir(SetDirectories.get_temp_dir() + "/padj")
            edgelist = Pipeline.make_edgelist(padj_file)

            os.chdir(output_dir)
            output_filename = f"{padj_file[:-len('_no_blacklist_no_cytobands_padj.txt')]}_edgelist.txt"
            with open(os.path.join(output_dir, output_filename), "w") as f:
                for line in edgelist:
                    f.write(line + "\n")


def run_pipeline():
    """
    Call selected methods of the Pipeline, in the order specified
    """

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
    print(f"Pipeline completed in {end_time - start_time:.2f} seconds.")

if __name__ == "__main__":
    run_pipeline()
