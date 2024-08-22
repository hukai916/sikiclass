#!/usr/bin/env python

"""
Generate a table summarizing indel size distribution: min, max, and mean.

Usage:
    python parse_indel_distribution.py raw_csv out_tsv
"""

#%%
from genericpath import exists
import numpy as np 
import os
import pandas as pd
import sys

csv = sys.argv[1]
tsv = sys.argv[2]
sample = str(sys.argv[3])

# csv = "/Users/kaihu/Projects/Nathan/work/sikiproject/new_strategy/h3f3d_umi_correct_cutoff_5/workflow/res/stats/indel_ref_pos/without_tag/SpCas9_gspPCR_SIKI2_1_1.csv"

#%%
df = pd.read_csv(csv)
df_del = df[df["indel"] < 0]
df_insertion = df[df["indel"] > 0]

#%%
del_min = min(abs(df_del["indel"]))
del_max = max(abs(df_del["indel"]))
del_mean = np.mean(abs(df_del["indel"]))

insertion_min = min(abs(df_insertion["indel"]))
insertion_max = max(abs(df_insertion["indel"]))
insertion_mean = np.mean(abs(df_insertion["indel"]))

os.makedirs(os.path.dirname(tsv), exist_ok = True)
with open(tsv, "w") as f:
    f.write("\t".join([sample, str(del_max), str(del_min), str(del_mean), str(insertion_max), str(insertion_min), str(insertion_mean)]))
    f.write("\n")