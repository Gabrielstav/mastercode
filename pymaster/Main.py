# Use fit hic instead of NCHG?
# Find papaer and read it.

# Set the path to the Fit-Hi-C executable
# fit_hic_executable = "path/to/fitHiC.py"
#
# # Define the input files and parameters
# hicpro_contact_matrix = "path/to/hicPro_contact_matrix.txt"
# fit_hic_output = "path/to/fit_hic_output.txt"
# resolution = 20000  # Set the desired resolution, e.g., 20 kb

# Run Fit-Hi-C with the specified parameters'
# subprocess.run([
#     "python",
#     fit_hic_executable,
#     "-f", hicpro_contact_matrix,
#     "-o", fit_hic_output,
#     "-r", str(resolution),
#     # Add any other required parameters
# ])

# Chromosome parallelization for NCHG:

# with open(file in input_dir) as input_files:

# chr_dict = {}
# or list with comprehension? tuple? set? queue? stack? tree? graph? array? Think dict.
# for file in input_files:
# find fast way to filter to chromosomes, regex? What is the best way to grep for changes in chr string, if col[0] == chr1?
# Can also filter on inter/intra here but whY?
# Then find way to split the file, use temp files? The input to nchg has to be files infortunately.
# The big thing again though is concatenating the output files into one file (per resolution), before passing the whole genome file to padj.

# Research first: Split files into temp files based on chromosome (fastest), base temp file name on chromosome + file anme (res and exp).
# The new files can be output to a new directory, dict or filter etc to make them per chromosome, and these files are then passed to the NCHG method by just changin the input dir.
# The NCHG method doesn't need to change, jus the input to nchg, since we need to take the input files, and then split the, by chrom, then write to new dir, and then save the output of nchg to new dir,
# and then concatenate the output files into one file (per resolution).

# Pass these files to the NCHG method
# Concatenate the output files into one file (per resolution), in the correct order.

# How does this work with the parallelization? Can we have thread and one process pool for the method? The NCHG method needs to be called with process pool.
# Try to implement process pool for "find siginificant interactions" method.
