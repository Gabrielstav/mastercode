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




# To generate .cool files from BEDPE count files using cooler cload pairs in command line:
# cooler cload pairs -c1 1 -p1 2 -c2 4 -p2 5 --field count=7:dtype=int --zero-based
# /Users/GBS/Master/Reference/hg19/chrom_hg19.sizes:1000000 /Users/GBS/Master/HiC-Data/testing_chrom_parallelization/output/temp_dir/bedpe/mcf10_1000000.bedpe matrix_mcf10.cool


class Dirs:

    root_path = path.Path("/Users/GBS/Master/HiC-Data")
    base_path_figures = path.Path("/Users/GBS/Master/Figures")

    compartments_path = root_path / "compartments"
    bedpe_path_intra = root_path / "bedpe_files/bedpe/intra"
    bedpe_path_inter = root_path / "bedpe_files/bedpe/inter"
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
        bedpe_file,
        cool_file
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


# Compartment calling ??

# def balance_cooler_file(cool_file):
#     # Load the .cool file
#     c = cooler.Cooler(cool_file)
#
#     # Check if the cooler file is already balanced
#     if "weight" in c.bins().columns:
#         print(f"The cooler file {cool_file} is already balanced.")
#     else:
#         # Balance the cool file
#         cooler.balance.balance_cooler(c)
#
#     return c

def balance_cooler_file(cool_file, balanced_file):
    # Load the .cool file
    c = cooler.Cooler(cool_file)

    # Check if the cooler file is already balanced
    if "weight" in c.bins().columns:
        print(f"The cooler file {cool_file} is already balanced.")
        return c
    else:
        # Balance the cool file
        weights, balanced_c = cooler.balance.balance_cooler(c)

        # Create a copy of bins dataframe with new weights
        bins = c.bins()[:]
        bins['weight'] = weights

        # Fetch the necessary columns explicitly from PixelTable
        pixels = c.pixels()[:]
        pixels = pixels[['bin1_id', 'bin2_id', 'count']]

        # Map the weights from bins to pixels
        pixels = pixels.join(bins['weight'], on='bin1_id', how='left')
        pixels.rename(columns={'weight': 'weight1'}, inplace=True)
        pixels = pixels.join(bins['weight'], on='bin2_id', how='left')
        pixels.rename(columns={'weight': 'weight2'}, inplace=True)
        pixels['weight'] = pixels['weight1'] * pixels['weight2']

        print(pixels.head())
        print(bins.head())

        # Ensure the correect cols are present
        assert "bin1_id" in pixels.columns
        assert "bin2_id" in pixels.columns
        assert "count" in pixels.columns
        assert "weight" in pixels.columns

        pixels_copy = pixels.copy()
        print(pixels_copy.head())


        # Create a new cooler file with balanced weights
        with h5py.File(balanced_file, 'w') as h5:
            cooler.create_cooler(h5, "/", bins, pixels)

        print(f"The cooler file {cool_file} is now balanced and saved to {balanced_file}.")

    # Load the balanced cooler file
    return cooler.Cooler(balanced_file)

def call_compartments(cool_file, balanced_file, output_file):
    # Convert Path objects to string
    cool_file = str(cool_file)
    output_file = str(output_file)

    # Balance the cooler file
    c = balance_cooler_file(cool_file, balanced_file)

    # Create a DataFrame to hold the results
    result = pd.DataFrame(columns=["chrom", "start", "end", "compartment", "eigenvector"])

    # Loop over each chromosome
    for chrom in c.chromnames:
        # Fetch the contact matrix for the chromosome
        matrix = c.matrix(balance=True).fetch(chrom)

        # Compute the correlation matrix
        correlation_matrix = np.corrcoef(matrix)

        # Compute the first principal component
        pca = PCA(n_components=1)
        principalComponents = pca.fit_transform(correlation_matrix)

        # Assign compartments based on the sign of PC1
        compartments = np.sign(principalComponents)

        # Create a DataFrame for this chromosome
        df = pd.DataFrame({
            "chrom": chrom,
            "start": c.bins().fetch(chrom)["start"],
            "end": c.bins().fetch(chrom)["end"],
            "compartment": ['A' if x > 0 else 'B' for x in compartments.flatten()],
            "eigenvector": principalComponents.flatten()  # save the eigenvectors
        })

        # Append the results to the main DataFrame
        result = result.append(df)

    # Save the results to a BED file
    result.to_csv(output_file, sep='\t', header=False, index=False)

call_compartments(Dirs.cool_path_intra / "imr90_1000000.cool", Dirs.cool_path_intra / "imr90_1000000_balanced.cool", Dirs.compartments_path / "imr90_1mb_compartments.bed")
