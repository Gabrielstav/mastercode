# Import modules
import os as os
import cooler
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import subprocess as sp
import pathlib as path
from sklearn.decomposition import PCA
import h5py
import time
from scipy.stats import pearsonr





# To generate .cool files from BEDPE count files using cooler cload pairs in command line:
# cooler cload pairs -c1 1 -p1 2 -c2 4 -p2 5 --field count=7:dtype=int --zero-based
# /Users/GBS/Master/Reference/hg19/chrom_hg19.sizes:1000000 /Users/GBS/Master/HiC-Data/testing_chrom_parallelization/output/temp_dir/bedpe/mcf10_1000000.bedpe matrix_mcf10.cool


class Dirs:

    root_path = path.Path("/Users/GBS/Master/HiC-Data")
    base_path_figures = path.Path("/Users/GBS/Master/Figures")

    chrom_sizes_file = "/Users/GBS/Master/Reference/hg19/chrom_hg19.sizes"

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
        f"{chrom_sizes_file}:{bin_size}",
        str(bedpe_file),
        str(cool_file)
    ]

    sp.run(command, check=True)



# bedpe_to_cool("/Users/GBS/Master/HiC-Data/bedpe_files/bedpe_intra/imr90_1000000.bedpe", "/Users/GBS/Master/Reference/hg19/chrom_hg19.sizes", 1000000, "/Users/GBS/Master/HiC-Data/coolfiles/intra/imr90_1000000.cool")


def plot_heatmap(cool_file, region=None, saveas=None):
    # Load the Cooler file
    c = cooler.Cooler(cool_file)

    # Define the region to plot
    if region is not None:
        matrix = c.matrix(balance=False).fetch(region)
    else:
        matrix = c.matrix(balance=False)

    # Calculate the limits for the color scale
    vmin, vmax = np.percentile(matrix.data, [1, 99])

    # Define colormap
    cmap = mpl.cm.get_cmap("YlOrRd")
    cmap.set_bad("white")

    # Plot the heatmap
    plt.matshow(matrix, cmap=cmap, norm=colors.PowerNorm(gamma=0.3, vmin=vmin, vmax=vmax))
    plt.title("Heatmap of chromosome 2")

    plt.xlabel("Position (Mb)")
    plt.ylabel("Position (Mb)")

    ax = plt.gca()  # Get current axes

    # Major and minor ticks
    binsize_Mb = c.binsize / 1_000_000  # convert binsize to Mb
    num_bins = matrix.shape[0]
    bin_range_Mb = num_bins * binsize_Mb

    # calculate the tick positions in bin units
    major_ticks = np.arange(0, num_bins, 50 / binsize_Mb)  # every 50 Mb
    minor_ticks = np.arange(0, num_bins, 10 / binsize_Mb)  # every 10 Mb

    # remove the first major tick
    major_ticks = major_ticks[1:]

    # Set the tick locations
    ax.set_xticks(major_ticks, minor=False)
    ax.set_xticks(minor_ticks, minor=True)
    ax.set_yticks(major_ticks, minor=False)
    ax.set_yticks(minor_ticks, minor=True)

    # Label the major ticks
    ax.set_xticklabels([f"{x * binsize_Mb:.0f}" for x in major_ticks])
    ax.set_yticklabels([f"{x * binsize_Mb:.0f}" for x in major_ticks])

    if saveas:
        plt.savefig(saveas, dpi=1000, format="png")

    plt.show()


# plot_heatmap("/Users/GBS/Master/HiC-Data/coolfiles/intra/imr90_120000.cool", region="chr2", saveas=Dirs.heatmap_dir_intra / "imr_120000_chr1_0:120Mb.png")


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

def create_cooler_file():
    bedpe_file = Dirs.bedpe_path_intra / "gsm2824367_1000000.bedpe"
    chrom_file = Dirs.chrom_sizes_file
    resolution = 1000000
    output_file = Dirs.cool_path_intra / "gsm2824367_1000000.cool"

    bedpe_to_cool(bedpe_file, chrom_file, resolution, output_file)

# create_cooler_file()

def compartment_calling():
    input_cool_file = Dirs.cool_path_intra / "gsm_1000000.cool"
    balanced_cool_file = Dirs.cool_path_intra / "gsm_1000000_balanced.cool"
    output_file = Dirs.compartments_path / "gsm_1mb_compartments.bed"

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

compare_compartments()