import os as os



# Place to read in files from custom dir, containing all HiC-Pro output files
# This code will read in all these files, split them on Celline/Chromosome and resolution, and store this information.
# Files will be handed 1 by 1 to the Pipeline, which after its run will create two directories,
# one containing the edge list for that specific cell line/chromosome and resolution, and one containing the
# output from pipeline for each step like BEDPE, NCHG out etc

# Abstract class for file handling


class FiletoPipeline:

    @staticmethod
    def read_files(input_directory):
        """
        This method should be handed an input directory (root dir containing HiC-Pro output files)
        :param input_directory: OS path to input directory
        :return: OS path to files
        """

        bedfiles = []

        for dirs, files in os.walk(input_directory):
            for file in files:
                if file.endswith(".bed"):
                    yield os.path.join(dirs, file)
                    bedfiles.append(os.path.join(dirs, file))

        return bedfiles

    @staticmethod
    def send_to_pipeline(input_directory):

        """
        This method sends the files to the pipeline one by one
        :param input_directory:
        :return:
        """

        for file in FiletoPipeline.read_files(input_directory):
            Pipeline(file)

        return "Done"


class FilefromPipeline:

    @staticmethod
    def pipeline_return():
        """
        This method returns the output from the pipeline
        :return: os path to output files
        """

        return Pipeline.output_files

    @staticmethod
    def files_to_processing():
        """
        This method sends the files to the processing class
        :return:
        """

        for file in FilefromPipeline.pipeline_return():
            Processing(file)

        return "Done"






class FileHandler:

    def __init__(self, input_directory, output_directory):
        self.input_directory = input_directory
        self.output_directory = output_directory

    def read_files(self):
        for root, dirs, files in os.walk(self.input_directory):
            for file in files:
                if file.endswith(".bed"):
                    FileHandler(os.path.join(root, file))

    def read_files_as_list(self):
        file_list = []
        for root, dirs, files in os.walk(self.input_directory):
            for file in files:
                if file.endswith(".bed"):
                    file_list.append(os.path.join(root, file))
        return file_list

    def write_files(self):
        for root, dirs, files in os.walk(self.output_directory):
            for file in files:
                if file.endswith(".bed"):
                    FileHandler(os.path.join(root, file))

    def write_files_as_list(self):
        file_list = []
        for root, dirs, files in os.walk(self.output_directory):
            for file in files:
                if file.endswith(".bed"):
                    file_list.append(os.path.join(root, file))
        return file_list


class DataReader:

    """
    This class should be handed an input file (root dir containing HiC-Pro output files),
    and then read in all files in this directory, and split them on cell line, chromosome and resolution.
    As outout it should, one at a time, give the paths of the files to the Pipeline class, which will
    then run the pipeline on each file, and write the output to a new directory.

    Use OS module walk and split to get the paths to the files.
    Then filter the files, rename them? or store them as objects? or in a dict?
    Then give one file path at a time to Pipeline class, refactor it so that this can work.
    We need to keep the resolution in the file name, so that we can refactor the Processing class
    to handle different resolutions and filter on them.

    Må vel også sette output directory, since after the Pipeline is done, we have a lot of output files.
    So we only need to supply input dir(s) and output dir.
    It should for any dir traverse the structure and grab .BED files and store them.
    It needs to be able to handle multiple input dirs.
    """

    def __init__(self, input_directory, output_directory):
        self.input_directory = input_directory
        self.output_directory = output_directory

    def read_files(self):
        for root, dirs, files in os.walk(self.input_directory):
            for file in files:
                if file.endswith(".bed"):
                    DataReader(os.path.join(root, file))

    def read_files_as_list(self):
        file_list = []
        for root, dirs, files in os.walk(self.input_directory):
            for file in files:
                if file.endswith(".bed"):
                    file_list.append(os.path.join(root, file))
        return file_list




# Example of HIC-Pro diretory:
# /Users/GBS/Master/Pipeline/python_pipe_test/diff_bedpe_res/HIC_res

# Example of full path to HiC-Pro output file:
# /Users/GBS/Master/Pipeline/python_pipe_test/diff_bedpe_res/HIC_res/hic_results/matrix/sample1/raw/50000/sample1_50000_abs.bed


class TargetDirectories:

    hicpro_output_directory = None
    pipeline_output_directory = None

    @staticmethod
    def output_directory():
        if TargetDirectories.hicpro_output_directory is not None:
            return TargetDirectories.hicpro_output_directory
        while True:
            path = input("Please enter the full path to the directory containing HiC-Pro output files: ")
            if os.path.exists(path):
                TargetDirectories.hicpro_output_directory = path
                return path
            print("Invalid path. Please enter a valid path.")

    @staticmethod
    def pipeline_directory():
        if TargetDirectories.pipeline_output_directory is not None:
            return TargetDirectories.pipeline_output_directory
        while True:
            path = input("Please enter the full path to the directory you want the output files to be written to: ")
            if os.path.exists(path):
                TargetDirectories.pipeline_output_directory = path
                return path
            print("Invalid path. Please enter a valid path.")

class FileStruct:

    @staticmethod
    def generate_file_struct():
        for dirnames in os.walk(TargetDirectories.output_directory()):
            if dirnames in TargetDirectories.output_directory() == "hic_results":
                print(os.listdir())


class FileHandler:

    """This class needs to find the hic pro output files, and split them on name (celline, chromosome)
    then it needs to open these files, and split again on resolution, and store this information.
    so to end up we get cell line x and cell line y, with chromosome 1 and 2, and resolution 50000 and 100000.
    as a result we get 4 files, 2 for each cell line, and 2 for each resolution.
    these files will be handed 1 by 1 to the pipeline, which will create two directories,
    the directiories for the output should be handled here, not in the pipeline I think?
    """

    @staticmethod
    def write_lines_to_file(file_path, lines, batch_size=1000):
        with open(file_path, "w") as f:
            for i in range(0, len(lines), batch_size):
                f.writelines(lines[i:i + batch_size])
        return file_path

    @staticmethod
    def locate_root_dir():
        return TargetDirectories.output_directory()

    @staticmethod
    def locate_data():
        for dirnames in os.walk(FileHandler.locate_root_dir()):
            if dirnames in FileHandler.locate_root_dir() == "hic_results":
                print(os.listdir())



# Might work, idk if this is a good approach?
class TempFileManager:

    # stores temp files and the class name that generated them
    temp_files = {}

    @classmethod
    # Call this method to register temp files: TempFileManager.register_temp_file(file_path, FileHandler.__name__)
    def register_temp_file(cls, file_path, class_name):
        cls.temp_files[file_path] = class_name

    @classmethod
    # Call this method to delete temp files: TempFileManager.delete_temp_files(file_path)
    def delete_temp_files(cls):
        for file_path in cls.temp_files:
            os.unlink(file_path)
        cls.temp_files.clear()

    @classmethod
    def get_temp_files(cls):
        return cls.temp_files.copy()

# TargetDirectories._hicpro_output_directory() # Noe med calling som er feil, kanskje feil usecase, må vi ha instanser?

