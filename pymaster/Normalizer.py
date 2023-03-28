# Import modules
import matplotlib.pyplot as plt
from matplotlib import colors
from iced import normalization
from iced import datasets
# import HiNT-Package as hnt, does not work in IDE



# TODO: ICE normalization of CNV for cancer cell lines (CAIC or LOIC), unsure if I should do it here or on the HPC?
#       Probably try to do it on the HPC, the iced package is deprecated and throws errors, it'll take a while to process data locally, and I need to download the datasets as well.
#       So for now this is just a backup plan. Might now work on HPC either tho.


# Trying test data from iced package
counts, lengths, cnv = datasets.load_sample_cancer()

loic_normed = normalization.ICE_normalization(counts, counts_profile=cnv)
print(loic_normed)
block_biases = normalization.ICE_normalization(counts, lengths, cnv)
caic_normed = loic_normed / block_biases
print(caic_normed)


# Testing plotting
chromosomes = ["I", "II", "III", "IV", "V", "VI"]

fig, axes = plt.subplots(ncols=3, figsize=(14, 3))

axes[0].imshow(counts, cmap="RdBu_r", norm=colors.SymLogNorm(1),
               origin="bottom",
               extent=(0, len(counts), 0, len(counts)))

[axes[0].axhline(i, linewidth=1, color="#000000") for i in lengths.cumsum()]
[axes[0].axvline(i, linewidth=1, color="#000000") for i in lengths.cumsum()]
axes[0].set_title("Raw contact counts", fontweight="bold")

m = axes[1].imshow(block_biases, cmap="RdBu_r", norm=colors.SymLogNorm(1),
                   origin="bottom",
                   extent=(0, len(counts), 0, len(counts)))
[axes[1].axhline(i, linewidth=1, color="#000000") for i in lengths.cumsum()]
[axes[1].axvline(i, linewidth=1, color="#000000") for i in lengths.cumsum()]
axes[1].set_title("Estimated block biases", fontweight="bold")

m = axes[2].imshow(caic_normed,
                   cmap="RdBu_r", norm=colors.SymLogNorm(1),
                   origin="bottom",
                   extent=(0, len(counts), 0, len(counts)))
[axes[2].axhline(i, linewidth=1, color="#000000") for i in lengths.cumsum()]
[axes[2].axvline(i, linewidth=1, color="#000000") for i in lengths.cumsum()]
cb = fig.colorbar(m)
axes[2].set_title("Normalized contact counts with CAIC", fontweight="bold")

