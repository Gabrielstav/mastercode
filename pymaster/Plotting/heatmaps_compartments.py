# Import modules
import os as os
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import subprocess as sp
import pathlib as path
from sklearn.decomposition import PCA
import h5py
import time
from scipy.stats import pearsonr
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
import cooler
from scipy.stats import spearmanr


# To generate .cool files from BEDPE count files using cooler cload pairs in command line:
# cooler cload pairs -c1 1 -p1 2 -c2 4 -p2 5 --field count=7:dtype=int --zero-based
# /Users/GBS/Master/Reference/hg19/chrom_hg19.sizes:1000000 /Users/GBS/Master/HiC-Data/testing_chrom_parallelization/output/temp_dir/bedpe/mcf10_1000000.bedpe matrix_mcf10.cool


class Dirs:

    root_path = path.Path("/Users/GBS/Master/HiC-Data")
    base_path_figures = path.Path("/Users/GBS/Master/Figures")

    chrom_sizes_file = "/Users/GBS/Master/Reference/hg19/chrom_hg19.sizes"
    hg19_cytoband_file = "/Users/GBS/Master/Reference/hg19/cytoBand_hg19.txt"

    compartments_path = root_path / "compartments"
    bedpe_path_intra = root_path / "bedpe_files/bedpe_intra"
    bedpe_path_inter = root_path / "bedpe_files/bedpe_inter"
    cool_path_intra = root_path / "coolfiles/intra"
    cool_path_inter = root_path / "coolfiles/inter"

    heatmap_dir_intra = base_path_figures / "heatmaps/intra"
    distance_dir_intra= base_path_figures / "decay_distances"
    heatmap_dir_inter = base_path_figures / "heatmaps/inter"
    distance_dir_inter = base_path_figures / "decay_distances/inter"


def bedpe_to_cool(bedpe_file, chrom_sizes_file, bin_size, cool_file):
    """
    Convert a BEDPE file to a COOL file using the Cooler command-line tool.

    Parameters:
    bedpe_file (str): The path to the BEDPE file.
    chrom_sizes_file (str): The path to the chrom_sizes file.
    bin_size (int): The bin size in base pairs.
    cool_file (str): The path to the output COOL file.
    """
    command = [
        "cooler",
        "cload",
        "pairs",
        "-c1",
        "1",
        "-p1",
        "2",
        "-c2",
        "4",
        "-p2",
        "5",
        "--field",
        "count=7:dtype=int",
        "--zero-based",
        f"{str(chrom_sizes_file)}:{str(bin_size)}",
        str(bedpe_file),
        str(cool_file)
    ]

    sp.run(command, check=True)



# bedpe_to_cool("/Users/GBS/Master/HiC-Data/bedpe_files/bedpe_intra/imr90_120000.bedpe", "/Users/GBS/Master/Reference/hg19/chrom_hg19.sizes", 1000000, "/Users/GBS/Master/HiC-Data/coolfiles/intra/imr90_1000000.cool")





def plot_contact_decay(cool_path, log_bins=True, bin_num=100):
    """
    Generate contact decay plot for Hi-C data.

    Parameters:
    cool_path: string
        Path to the .cool file.
    log_bins: Boolean
        - If True, use logarithmic bins.
        - If False, use linear bins.
    bin_num: The number of bins to use

    Returns:
    Matplotlib figure
    """

    # Load the Cooler file
    c = cooler.Cooler(cool_path)

    # Get the genomic positions of the bins
    bin_positions = c.bins()[:]["start"].values

    # Get the pairs DataFrame
    pairs = c.pixels()[:]

    # Calculate the genomic distances for each interaction
    pairs["genomic_distance"] = np.abs(bin_positions[pairs["bin2_id"]] - bin_positions[pairs["bin1_id"]])

    # If log_bins is True, calculate log bin edges
    if log_bins:
        bin_edges = np.logspace(np.log10(pairs["genomic_distance"].min()), np.log10(pairs["genomic_distance"].max()), num=bin_num)
        pairs['bin'] = pd.cut(pairs["genomic_distance"], bins=bin_edges)
        pairs.dropna(subset=["bin"], inplace=True)  # Drop any pairs not falling into a bin
        pairs["bin_center"] = pairs["bin"].apply(lambda x: (x.right + x.left) / 2)

    # Otherwise, use linear bins
    else:
        pairs['bin_center'] = pd.cut(pairs['genomic_distance'], bins=bin_num, labels=False)

    # Compute average interaction frequency for each genomic distance bin
    average_interactions = pairs.groupby("bin_center")["count"].mean()

    # Plot
    plt.figure(figsize=(10, 7))
    plt.loglog(average_interactions.index.values, average_interactions.values)
    plt.title("Contact Decay")
    plt.xlabel("Genomic Distance (bp)")
    plt.ylabel("Average Interaction Frequency")
    plt.show()

# plot_contact_decay("/Users/GBS/matrix_mcf10.cool", log_bins=True, bin_num=100)


def balance_cooler_file(cool_file):
    """Load and balance a cooler file."""
    c = cooler.Cooler(str(cool_file))
    if "weight" in c.bins().columns:
        print(f"The cooler file {cool_file} is already balanced.")
        return c
    else:
        weights, _ = cooler.balance.balance_cooler(c)
        bins = c.bins()[:]
        bins['weight'] = weights
        pixels = c.pixels()[:]
        pixels = pixels[['bin1_id', 'bin2_id', 'count']]
        pixels = pixels.join(bins['weight'], on='bin1_id', how='left').rename(columns={'weight': 'weight1'})
        pixels = pixels.join(bins['weight'], on='bin2_id', how='left').rename(columns={'weight': 'weight2'})
        pixels['weight'] = pixels['weight1'] * pixels['weight2']
        return c, bins, pixels

def validate_and_save_cooler(bins, pixels, balanced_file):
    """Save a balanced cooler file and validate it."""
    cooler.create_cooler(cool_uri=str(balanced_file), bins=bins, pixels=pixels)
    time.sleep(5)
    if os.path.exists(balanced_file):
        if h5py.is_hdf5(str(balanced_file)):
            print("Balanced cooler file was successfully created and is a valid HDF5 file.")
        else:
            print("File exists but is not a valid HDF5 file.")
    else:
        print("File does not exist.")
    print(f"The cooler file is now balanced and saved to {balanced_file}.")

def create_cooler_file():
    bedpe_file = Dirs.bedpe_path_intra / "imr90_120000.bedpe"
    chrom_file = Dirs.chrom_sizes_file
    resolution = 120000
    output_file = Dirs.cool_path_intra / "imr90_120000.cool"

    bedpe_to_cool(bedpe_file, chrom_file, resolution, output_file)

# create_cooler_file()

def balance_and_save_cooler(cool_file, balanced_file):
    """Balance and save a cooler file."""
    c, bins, pixels = balance_cooler_file(cool_file)
    validate_and_save_cooler(bins, pixels, balanced_file)


# balance_and_save_cooler("/Users/GBS/Master/HiC-Data/coolfiles/intra/mcf7_50000.cool", "/Users/GBS/Master/HiC-Data/coolfiles/intra/mcf7_1000000_balanced.cool")

def process_chromosome(chrom, cooler_obj, bins_df):
    """
    Process a single chromosome: compute its correlation matrix,
    perform PCA, and assign compartments
    """

    # Fetch the contact matrix for the chromosome
    chrom_matrix = cooler_obj.matrix(balance=True).fetch(chrom)

    # Skip chromosome if it only contains one bin
    if chrom_matrix.size == 1:
        print(f"Skipping chromosome {chrom} because it only contains one bin.")
        return None

    # Compute the correlation matrix and replace NaN/infinite values
    correlation_matrix = np.corrcoef(chrom_matrix)
    correlation_matrix = np.nan_to_num(correlation_matrix)

    # Compute the first principal component
    pca = PCA(n_components=1)
    principal_components = pca.fit_transform(correlation_matrix)

    # Assign compartments based on the sign of PC1
    compartments = np.sign(principal_components)

    # Create and return a DataFrame for this chromosome
    chrom_bins = bins_df.fetch(chrom)
    chrom_df = pd.DataFrame({
        "chrom": chrom,
        "start": chrom_bins["start"],
        "end": chrom_bins["end"],
        "compartment": ['A' if x > 0 else 'B' for x in compartments.flatten()],
        "eigenvector": principal_components.flatten()  # save the eigenvectors
    })

    return chrom_df


def call_compartments(cool_file, balanced_file, output_file):
    """
    Call chromatin compartments using Hi-C data in a cooler file
    """

    # Convert Path objects to string
    cool_file = str(cool_file)
    output_file = str(output_file)

    # Balance the cooler file
    c, bins, pixels = balance_cooler_file(cool_file)
    validate_and_save_cooler(bins, pixels, balanced_file)

    # Load the balanced cooler file
    c_balanced = cooler.Cooler(str(balanced_file))

    # Create a DataFrame to hold the results
    result = pd.DataFrame(columns=["chrom", "start", "end", "compartment", "eigenvector"])

    # Loop over each chromosome
    for chrom in c_balanced.chromnames:
        if chrom == 'chrM':
            print("Skipping mitochondrial chromosome (chrM).")
            continue
        chrom_df = process_chromosome(chrom, c_balanced, c_balanced.bins())
        if chrom_df is not None:
            result = result.append(chrom_df)

    # Save the results to a BED file
    result.to_csv(output_file, sep='\t', header=False, index=False)





def compartment_calling():
    input_cool_file = Dirs.cool_path_intra / "imr90_1000000.cool"
    balanced_cool_file = Dirs.cool_path_intra / "imr90_1000000_balanced_2.cool"
    output_file = Dirs.compartments_path / "imr90_1Mb_pca1_compartments.bed"

    cooler_obj, bins, pixels = balance_cooler_file(input_cool_file)
    validate_and_save_cooler(bins, pixels, balanced_cool_file)

    call_compartments(input_cool_file, balanced_cool_file, output_file)

# compartment_calling()

def compute_pearson_and_plot(bed_file1, bed_file2, cell_line1, cell_line2):
    # Load BED files into pandas DataFrames
    df1 = pd.read_csv(bed_file1, sep='\t', header=None)
    df2 = pd.read_csv(bed_file2, sep='\t', header=None)

    # Extract the eigenvector values and compartment calls
    eigenvector1 = df1.iloc[:, 4]
    eigenvector2 = df2.iloc[:, 4]
    compartments1 = df1.iloc[:, 3]
    compartments2 = df2.iloc[:, 3]

    # Compute the Pearson correlation coefficient
    correlation = pearsonr(eigenvector1, eigenvector2)[0]

    # Create a new column in the DataFrame indicating whether the compartment call changed
    compartment_change = compartments1 != compartments2

    # Create a scatter plot with points colored based on the compartment call
    plt.figure(figsize=(8, 6))
    plt.scatter(eigenvector1[compartment_change], eigenvector2[compartment_change], color='red', alpha=0.5, label='Compartment change')
    plt.scatter(eigenvector1[~compartment_change], eigenvector2[~compartment_change], color='blue', alpha=0.5, label='No compartment change')
    plt.xlabel(f'Eigenvector - {cell_line1}')
    plt.ylabel(f'Eigenvector - {cell_line2}')
    plt.title(f'Correlation between {cell_line1} and {cell_line2}: r = {correlation:.2f}')
    plt.legend()
    plt.show()

    # Return the correlation
    return correlation

def compare_compartments():
    # Load the BED files
    bed_file1 = Dirs.compartments_path / "mcf7_1mb_compartments.bed"
    bed_file2 = Dirs.compartments_path / "mcf10_1mb_compartments.bed"

    # Compute the Pearson correlation coefficient
    correlation = compute_pearson_and_plot(bed_file1, bed_file2, 'MCF7', 'MCF10')

# compare_compartments()

def merge_bed_files(bed_file1, bed_file2):
    # Load BED files into pandas DataFrames
    df1 = pd.read_csv(bed_file1, sep='\t', header=None, names=['chrom', 'start', 'end', 'compartment', 'eigenvector'])
    df2 = pd.read_csv(bed_file2, sep='\t', header=None, names=['chrom', 'start', 'end', 'compartment', 'eigenvector'])

    # Merge the two dataframes on chrom, start and end
    merged_df = pd.merge(df1, df2, on=['chrom', 'start', 'end'], suffixes=('_1', '_2'))

    return merged_df



def plot_compartments_on_heatmap_sns(cool_path, compartments_bed_path, chrom=None, saveas=None):
    clr = cooler.Cooler(str(cool_path))
    if chrom == "all":
        matrix = clr.matrix()[:]
    else:
        matrix = clr.matrix(balance=True).fetch(chrom)

    diag_sums = [np.nansum(np.diagonal(matrix, offset=i)) for i in range(matrix.shape[0])]
    diag_counts = [np.sum(~np.isnan(np.diagonal(matrix, offset=i))) for i in range(matrix.shape[0])]
    exp_const = 1e-10
    expected = np.array(diag_sums) - np.array(diag_counts) + exp_const
    print(matrix)
    print(expected)
    observed_over_expected = matrix / expected[:, None]
    print(observed_over_expected)

    corr = np.corrcoef(observed_over_expected, rowvar=False)
    compartments = pd.read_csv(compartments_bed_path, sep="\t", header=None, names=["chrom", "start", "end", "compartment", "score"])

    # Convert the genomic coordinates to bin coordinates
    binsize = clr.binsize
    if chrom != "all":
        compartments = compartments[compartments["chrom"] == chrom]
        compartments["start"] = compartments["start"] // binsize
        compartments["end"] = compartments["end"] // binsize

    # Generate the plot
    # fig, (ax1, ax2) = plt.subplots(nrows=2, figsize=(10, 10), gridspec_kw={'height_ratios': [1, 4]}, sharex="all")
    fig, ax2 = plt.subplots(figsize=(10, 10))

    # Plot PCa1
    # ax1.plot(compartments["start"], compartments["score"], color='k')

    # Plot the contact map
    vmax = np.percentile(corr[~np.isnan(corr)], 99)
    vmin = np.percentile(corr[~np.isnan(corr)], 1)

    # Make colors not blend:
    colormap = ["blue", "blue", "black", "red", "red"]
    gradient = [0.0, 0.3, 0.5, 0.7, 1.0]
    cmap = mcolors.LinearSegmentedColormap.from_list("cmap", list(zip(gradient, colormap)))
    sns.heatmap(corr, cmap=cmap, ax=ax2)

    # Create ticks at every 25 Mb
    tick_frequency = 25
    ticks = np.arange(0, matrix.shape[0], tick_frequency * 10 ** 6 // binsize)

    # Create tick labels with 'Mb' added
    tick_labels = [f"{tick * binsize // 10 ** 6} Mb" for tick in ticks]


    ax2.set_xticks(ticks[1:])
    ax2.set_xticklabels(tick_labels[1:])
    ax2.set_yticks(ticks[1:])
    ax2.set_yticklabels(tick_labels[1:])

    ax2.set_title(f"Contact map for {chrom}")
    # ax1.set_title(f"Compartments for {chrom}")
    plt.tight_layout()

    if saveas:
        plt.savefig(saveas, dpi=1000, format="png")

    plt.show()

# plot_compartments_on_heatmap_sns(Dirs.cool_path_intra / "mcf7_1000000_balanced.cool", Dirs.compartments_path / "mcf10_250kb_compartments.bed", chrom="chr2")  # :0-95000000")  # , saveas=Dirs.heatmap_dir_intra / "mcf10_250kb_chr2_compartments_plot.png")

def plot_heatmap2(cool_file, chrom=None, saveas=None, log=False):
    # Load the Cooler file
    c = cooler.Cooler(str(cool_file))

    # Define the region to plot
    if chrom == "all":
        matrix = c.matrix(balance=False)[:]
    else:
        matrix = c.matrix(balance=False).fetch(chrom)

    # Copy matrix for processing
    matrix_copy = matrix.copy()

    if log:
        matrix_copy = np.log1p(matrix_copy)
        matrix_copy[matrix_copy == -np.inf] = 0

    # Calculate the limits for the color scale
    vmin, vmax = np.percentile(matrix_copy, [0, 100])

    # Define colormap
    cmap = mpl.cm.get_cmap("YlOrRd")
    cmap.set_bad("white")

    # Plot the heatmap
    plt.matshow(matrix_copy, cmap=cmap, norm=mcolors.PowerNorm(gamma=0.2, vmin=vmin, vmax=vmax))
    plt.title(f"Heatmap of {chrom}")

    # Get current axes
    ax = plt.gca()

    # Create ticks at every 25 Mb
    tick_frequency = 25
    ticks = np.arange(0, matrix.shape[0], tick_frequency * 10 ** 6 // c.binsize)

    # Create tick labels with 'Mb' added
    tick_labels = [f"{tick * c.binsize // 10 ** 6} Mb" for tick in ticks]


    ax.set_xticks(ticks[1:])
    ax.set_xticklabels(tick_labels[1:])
    ax.set_yticks(ticks[1:])
    ax.set_yticklabels(tick_labels[1:])
    ax.xaxis.tick_bottom()
    plt.xticks(rotation='vertical')
    ax.tick_params(axis='x', which='major', pad=10)
    ax.tick_params(axis='x', labelsize=7)
    ax.tick_params(axis='y', labelsize=7)

    ax.set_title(f"Contact map for chromosome 2")
    plt.tight_layout()

    if saveas:
        plt.savefig(saveas, dpi=1000, format="png", bbox_inches='tight')


    plt.show()

# plot_heatmap2(Dirs.cool_path_intra / "imr90_250000_balanced.cool", chrom="chr2", saveas=Dirs.heatmap_dir_intra / "imr90_example250kb_chr2.png")