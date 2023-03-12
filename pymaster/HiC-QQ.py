# Import modules
import os.path

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import hicexplorer as he
import re as re
from Pipeline import SetDirectories


# TODO: Kinda low priority, but might be nice later?

# Read in BEDPE file
def bedpe_directory():
    norm_chr18_inc = os.path.abspath("/Users/GBS/Master/HiC-Data/Pipeline_out/chr18_INC/chr18_norm/temp_dir/bedpe")
    raw_chr18_inc = os.path.abspath("/Users/GBS/Master/HiC-Data/Pipeline_out/chr18_INC/chr18_raw/temp_dir/bedpe")
    return norm_chr18_inc, raw_chr18_inc

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

