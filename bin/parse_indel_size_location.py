#!/usr/bin/env python

"""
Generate a table summarizing indel size and relative location to PAM.

Usage:
    python parse_indel_size_location.py raw_tsv 1113 out_tsv
"""

#%%
from genericpath import exists
import sys 
import os
import pandas as pd
import csv

input_tsv = sys.argv[1]
pam = int(sys.argv[2])
tsv = sys.argv[3]

# tsv = "/Users/kaihu/Projects/Nathan/work/sikiproject/new_strategy/h3f3d_umi_correct_cutoff_5/workflow/res/stats/indel_ref_pos/without_tag/SpCas9_gspPCR_SIKI2_1_1.tsv"
# pam = 1113

#%%
df = pd.read_csv(input_tsv, sep = "\t")
if not df.shape[0] == 0:
    df.loc[:, "pam_pos"] = pam -1
    df.loc[:, "distance_to_pam"] = df.loc[:, "ref_pos"] - df.loc[:, "pam_pos"]
else:
    df["pam_pos"] = pd.NA
    df["distance_to_pam"] = pd.NA
    
os.makedirs(os.path.dirname(tsv), exist_ok=True)
df.to_csv(tsv, sep = "\t", index = False, quoting = csv.QUOTE_NONE, quotechar = "",  escapechar = "\t")