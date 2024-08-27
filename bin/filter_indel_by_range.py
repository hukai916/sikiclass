#!/usr/bin/env python

"""
Filter out indels that appear in given range.
Usage:

python filter_indel_by_range.py raw.tsv "1:2000" output.tsv 
"""

#%%
import sys 
import pandas as pd 
import csv
import os

raw_tsv = sys.argv[1]
scan_range = sys.argv[2].strip()
out_tsv = sys.argv[3]

# %%
df_indel = pd.read_csv(raw_tsv, sep = "\t")
os.makedirs(os.path.dirname(out_tsv), exist_ok = True)

if scan_range == ":":
    df_indel.to_csv(out_tsv, sep = "\t", index = None, quoting = csv.QUOTE_NONE, quotechar = "",  escapechar = "\t")
else:
    start, end = [int(x) for x in scan_range.split(":")]
    df_indel = df_indel[(df_indel['ref_pos'] >= start) & (df_indel['ref_pos'] <= end)]
    df_indel.to_csv(out_tsv, sep = "\t", index = None, quoting = csv.QUOTE_NONE, quotechar = "",  escapechar = "\t")
