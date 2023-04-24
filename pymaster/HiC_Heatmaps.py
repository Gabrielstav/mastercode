# Import modules
import os.path

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import hicexplorer as he
import re as re
from HiEdge import SetDirectories as Sd
import fithic as fh


a.load_data("/Users/GBS/Master/HiC-Data/Pipeline_out/chr18_INC/chr18_norm/temp_dir/bedpe/chr18_norm_250000.bedpe")
a.fit()
a.plot()



# TODO: Make plots of Hi-C heat maps to compare with LCC plots.

# Read in BEDPE file
def bedpe_directory():
    norm_chr18_inc = os.path.abspath("/Users/GBS/Master/HiC-Data/Pipeline_out/chr18_INC/chr18_norm/temp_dir/bedpe")
    raw_chr18_inc = os.path.abspath("/Users/GBS/Master/HiC-Data/Pipeline_out/chr18_INC/chr18_raw/temp_dir/bedpe")
    return norm_chr18_inc, raw_chr18_inc

def padj_edgelist_directory():
    pvals_chr18_raw = os.path.abspath("/Users/GBS/Master/HiC-Data/Pipeline_out/chr18_INC/chr18_raw/temp_dir/weighted_edgelists")
    pvals_chr18_iced = os.path.abspath("/Users/GBS/Master/HiC-Data/Pipeline_out/chr18_INC/chr18_norm/temp_dir/weighted_edgelists")
    return pvals_chr18_raw, pvals_chr18_iced

def pval_directory():
    pvals_chr18_raw = os.path.abspath("/Users/GBS/Master/HiC-Data/Pipeline_out/chr18_INC/chr18_raw/temp_dir/NCHG_output")
    pvals_chr18_iced = os.path.abspath("/Users/GBS/Master/HiC-Data/Pipeline_out/chr18_INC/chr18_norm/temp_dir/NCHG_output")
    return pvals_chr18_raw, pvals_chr18_iced

class import_bedpe:
    @staticmethod
    def get_bedpe_files():
        bedpe_files = []
        for files in bedpe_directory():
            if files.endswith(".bedpe"):
                bedpe_files.append(files)
        return bedpe_files

    @staticmethod
    def make_bedpe_df():
        for bedpe_file in import_bedpe.get_bedpe_files():
            bedpe_df = pd.read_csv(bedpe_file, sep="\t", header=None)
            return bedpe_df


class convert_bedpe:
    @staticmethod
    def to_hic_table(bedpe_file):
        hic_table = he.HiCTable.from_dataframe(
            import_bedpe.make_bedpe_df(),
            mode="full",
            resolution=int(re.search(r"(\d+)[^/\d]*$", bedpe_file).group(1)),
            chromosome1=18,  # need to automate this on file names?
            chromosome2=18,  # same as above, or set as flag?
            position1=1,  # same as below...
            position2=4,  # what col has position 2? But bedpe has 4 position cols
            count=6,  # what col has the IF vals
        )
        return hic_table

    @staticmethod
    def call_convert_bedpe():
        tables = []
        for bedpe_file in import_bedpe.get_bedpe_files():
            tables.append(convert_bedpe.to_hic_table(bedpe_file))  # how best to store as dfs?
        return tables


class plot_hic:
    pass

class plot_pvals:

    @staticmethod
    def get_pvals_files():
        pvals_files = []
        os.chdir(pval_directory()[0])
        for files in pval_directory():
            if files.endswith(".txt"):
                pvals_files.append(files)
        print(pvals_files)
        return pvals_files

    @staticmethod
    def make_pvals_df():
        for pvals_file in plot_pvals.get_pvals_files():
            pvals_df = pd.read_csv(pvals_file, sep="\t", header=None)
            print(pvals_df)
            return pvals_df

    @staticmethod
    def plot_pvals():
        pvals = []
        for pval in plot_pvals.make_pvals_df():
            pvals.append(pval[6])
        return pvals

plot_pvals.make_pvals_df()