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


# TODO: Enable reading in of different file structures (eg only raw folder, or only matrix folder)

####################################################
# Pre-processing pipeline for Hi-C data from HiC-Pro
####################################################

class SetDirectories:
    """
    SET INPUT-, OUTPUT- AND REFERENCE DIRS AND FULLPATH TO NCHG HERE

    For each run, change the input and output directories to the appropriate directories
    Set input dir to root dir containing HiC-Pro output folders (raw, matrix).
    Set normalized data = True to process ICE matrices, False to process raw data.
    """

    input_dir = os.path.abspath("/Users/GBS/Master/HiC-Data/HiC-Pro_out/chr18_inc/chr18_raw")
    output_dir = os.path.abspath("/Users/GBS/Master/HiC-Data/Pipeline_out/chr18_INC/chr18_raw_icedpipe")
    reference_dir = os.path.abspath("/Users/GBS/Master/Reference")
    nchg_path = os.path.abspath("/Users/GBS/Master/Scripts/NCHG_hic/NCHG")
    normalized_data = False

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
        for iced_matrixfile in iced_matrixfiles:
            resolution = int(iced_matrixfile.split("/")[-2])
            experiment = iced_matrixfile.split("/")[-4]
            key_prefix = "/".join(iced_matrixfile.split("/")[:-3])  # common part of file path excluding /raw/ and /iced/
            key_suffix = iced_matrixfile.split("/")[-1][:-len("_iced.matrix")]  # common part of file name excluding _iced.matrix
            key = f"{experiment, resolution}"

            # Group ICE-normalized matrix file
            for bedfile in bedfiles:
                bedfile_prefix = "/".join(bedfile.split("/")[:-3])
                bedfile_suffix = bedfile.split("/")[-1][:-len("_abs.bed")]
                if bedfile_prefix == key_prefix and bedfile_suffix == key_suffix:
                    if key not in grouped_iced_files:
                        grouped_iced_files[key] = (bedfile, iced_matrixfile)
                    else:
                        grouped_iced_files[key] += (bedfile, iced_matrixfile)

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

    # window_size = bedpe_file.strip(".bedpe").split("_")[2]

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

        # TODO: Fix this, there is something wrong with types in the NCHG script, it doesn't like the window size.
        # This could be because

        # Setting window size to be the same as the resolution
        window_size = int(re.search(r"(\d+)[^/\d]*$", bedpe_file).group(1))
        if not isinstance(window_size, int):
            raise ValueError(f"Window size must be an integer, {window_size} is not an integer.")

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
            # check if "iced" is in the file path (if the data is normalized)

            if "iced" in no_cytobands_file:
                # open the file and round the interaction values
                bedpe_inted = []
                with open(os.path.join(SetDirectories.get_temp_dir(), "no_cytobands", no_cytobands_file), "r") as file:
                    for line in file:
                        fields = line.split()
                        if len(fields) >= 7 and "." in fields[6]:  # rounds float iced matrix
                            fields[6] = str(int(round(float(fields[6]))))
                        bedpe_inted.append("\t".join(fields))

                # save the rounded bedpe file
                bedpe_inted_file = no_cytobands_file.replace(".bedpe", "_inted.bedpe")
                with open(os.path.join(SetDirectories.get_temp_dir(), "no_cytobands", bedpe_inted_file), "w") as file:
                    file.write("\n".join(bedpe_inted))

                # run the NCHG script on the rounded bedpe file
                os.chdir(SetDirectories.get_temp_dir() + "/no_cytobands")
                nchg_output = Pipeline.find_siginificant_interactions(bedpe_inted_file)

                # save the NCHG output to a file
                os.chdir(output_dir)
                output_filename = f"{no_cytobands_file[:-len('.bedpe')]}_nchg_output.txt"
                with open(os.path.join(output_dir, output_filename), "w") as f:
                    f.writelines(nchg_output)

            else:
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

        # TODO: Check if float points in pval of 0.0 fix yields correct results:
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

    @staticmethod
    def interactions_per_resolution(edge_list_file):
        """
        :input: edgelist files
        :output: Dictionary with resolution as key, and number of bins as value
        """

        resolution = int(re.search(r"(\d+)[^/\d]*$", edge_list_file).group(1))
        if not isinstance(resolution, int):
            raise ValueError(f"Resolution must be an integer, {resolution} is not an integer.")

        interaction_count = []
        if os.path.exists(edge_list_file):
            with open(edge_list_file, "r") as f:
                for line in f:
                    interaction_count.append(line)

            interactions_per_resolution = {resolution: len(interaction_count)}
            return interactions_per_resolution
        else:
            raise FileNotFoundError(f"File not found: {edge_list_file}")

    @staticmethod
    def call_interactions_per_resolution():

        # Create the output directory if it doesn't exist
        output_dir = os.path.join(SetDirectories.get_temp_dir(), "interactions_per_resolution")
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
        else:
            shutil.rmtree(output_dir)
            os.mkdir(output_dir)

        output_file = os.path.join(output_dir, "interactions_per_resolution.txt")

        # Get the list of edgelist files and sort them based on resolution
        edge_list_files = os.listdir(SetDirectories.get_output_dir() + "/edgelists")
        edge_list_files.sort(key=lambda x: int(re.search(r"(\d+)[^/\d]*$", edge_list_file).group(1)))

        # Iterate over the sorted list of edgelist files and process/save them
        with open(output_file, "w") as f:
            os.chdir(SetDirectories.get_output_dir() + "/edgelists")
            for edge_list_file in edge_list_files:
                interactions_per_resolution = Pipeline.interactions_per_resolution(os.path.join(SetDirectories.get_output_dir(), "edgelists", edge_list_file))
                for key, val in interactions_per_resolution.items():
                    f.write("{:<30} {:<10}{:<10}\n".format(edge_list_file, key, str(val)))

                # Check if the number of lines in the file matches the value of the dictionary key
                with open(edge_list_file) as edgelist:
                    line_count = sum(1 for _ in edgelist)
                    if line_count != val:
                        print(f"Warning: {edge_list_file} has {line_count} lines, expected {val}.")


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
        "input_to_make_weighted_edgelist",
        "call_interactions_per_resolution"
    ]

    # Call each method once
    for method in tqdm(method_names):
        if type(method) == str:  # Check if method name is a string
            method = getattr(Pipeline, method)  # Get the method reference
        method()  # Call the method

    # Print runtime on completion
    end_time = time.time()
    print(f"Pipeline completed in {end_time - start_time:.2f} seconds.")


run_pipeline()

