#!/usr/bin/env python

"""
Filter out indels that appear in uninjected samples.
Usage:

python filter_uninjected_indel.py uninjected_indel.csv input.csv output.csv 
"""

#%%
import sys 
import pandas as pd 
import csv

indel = sys.argv[1]
input = sys.argv[2]
output = sys.argv[3]
sample_id = sys.argv[4]
control_indel_freq_cutoff = int(sys.argv[5])
# indel = "../stats/summary_counts/uninjected_indel_with_any_insert.csv"
# input = "../stats/summary_counts/SpCas9_gspPCR_SIKI2_1_1/minimap2/fastq_with_any_insert.csv"
# output = "../stats/summary_counts/SpCas9_gspPCR_SIKI2_1_1/minimap2/fastq_with_any_insert_filter_uninjected.csv"
# sample_id = "test_sample"

# %%
df_indel = pd.read_csv(indel, sep = "\t")
# print(df_indel.head())
df_input = pd.read_csv(input, sep = "\t")
#%%
df_indel_grouped = df_indel.groupby(['ref_pos', 'ref_base', 'indel']).size().reset_index(name='count')
# filter uninjected indel with at least xxx occurrences:
df_indel_grouped = df_indel_grouped[df_indel_grouped['count'] >= control_indel_freq_cutoff]

# %%
df_set = set(zip(df_indel_grouped['ref_pos'], df_indel_grouped['ref_base'], df_indel_grouped['indel']))
mask = df_input.apply(lambda row: (row['ref_pos'], row['ref_base'], row['indel']) not in df_set, axis=1)

df_filtered = df_input[mask]
#%%
df_filtered.loc[:, "sample_id"] = [sample_id] * df_filtered.shape[0]
# %%
df_filtered.to_csv(output, sep = "\t", index = None, quoting = csv.QUOTE_NONE, quotechar = "",  escapechar = "\t")
# %%
