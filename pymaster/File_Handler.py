
# Place to read in files from custom dir, containing all HiC-Pro output files
# This code will read in all these files, split them on Celline/Chromosome and resolution, and store this information.
# Files will be handed 1 by 1 to the Pipeline, which after its run will create two directories,
# one containing the edge list for that specific cell line/chromosome and resolution, and one containing the
# output from pipeline for each step like BEDPE, NCHG out etc

# We need to wrap the NCHG script, and include the Gtrack and FDR script in the code or wrap them as well.

# Look at diagram for dir structure.


import os as os

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

class FileHandler:

    @staticmethod
    def read_hicpro_output():
        for dirnames in os.walk(TargetDirectories._hicpro_output_directory()):
            if dirnames in TargetDirectories._hicpro_output_directory() == "hic_results":
                print(os.listdir())




TargetDirectories._hicpro_output_directory() # Noe med calling som er feil, kanskje feil usecase, må vi ha instanser?
FileHandler.read_hicpro_output()



# class FileHandler:
#
#     def __init__(self, path):
#         self.path = path
#         self.cell_lines = []
#         self.chromosomes = []
#         self.resolutions = []
#         self.files = []
#         self.files_per_cell_line = {}
#         self.files_per_chromosome = {}
#         self.files_per_resolution = {}
#
#     def read_files(self):
#         for file in os.listdir(self.path):
#             if file.endswith(".bed") or file.endswith(".matrix"):
#                 self.files.append(file)
#
#     def split_files(self):
#         for file in self.files:
#             cell_line = file.split("_")[0]
#             chromosome = file.split("_")[1]
#             resolution = file.split("_")[2]
#
#             if cell_line not in self.cell_lines:
#                 self.cell_lines.append(cell_line)
#                 self.files_per_cell_line[cell_line] = []
#             if chromosome not in self.chromosomes:
#                 self.chromosomes.append(chromosome)
#                 self.files_per_chromosome[chromosome] = []
#             if resolution not in self.resolutions:
#                 self.resolutions.append(resolution)
#                 self.files_per_resolution[resolution] = []
#
#             self.files_per_cell_line[cell_line].append(file)
#             self.files_per_chromosome[chromosome].append(file)
#             self.files_per_resolution[resolution].append(file)
#
#     def run_pipeline(self):
#         for cell_line in self.cell_lines:
#             for chromosome in self.chromosomes:
#                 for resolution in self.resolutions:
#                     for file in self.files_per_cell_line[cell_line]:
#                         if file in self.files_per_chromosome[chromosome] and file in self.files_per_resolution[resolution]:
#                             # Run pipeline on file
#                             # Create dir for edge list and dir for output from pipeline
#                             # Write edge list to dir for edge list
#                             # Write output from pipeline to dir for output from pipeline
#                             # Move on to next file
#                             pass
#
#     def get_cell_lines(self):
#         return self.cell_lines
#
#     def get_chromosomes(self):
#         return self.chromosomes
#
#     def get_resolutions(self):
#         return self.resolutions
#
#     def get_files(self):
#         return self.files
#
#     def get_files_per_cell_line(self):
#         return self.files_per_cell_line
#
#     def get_files_per_chromosome(self):
#         return self.files_per_chromosome
#
#     def get_files_per_resolution(self):
#         return self.files_per_resolution
#
#
