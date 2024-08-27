#!/usr/bin/env python

#%%
import os
import pandas as pd
import sys
import csv

control_samples = [x.strip() for x in sys.argv[1].split(",")]
out_tsv = sys.argv[2]

#%%
all_files = os.listdir(os.getcwd())
control_tsvs = [file for file in all_files if file.endswith('.tsv') and any(sample in file for sample in control_samples)]

df = pd.DataFrame()

for tsv in control_tsvs:
    df_tem = pd.read_csv(tsv, sep = ",")
    df = pd.concat([df, df_tem], ignore_index = True)

df.to_csv(out_tsv, sep = '\t', index = False, quoting = csv.QUOTE_NONE, quotechar = "",  escapechar = "\t") 