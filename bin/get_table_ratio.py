#!/usr/bin/env python

#%%
import sys
import pandas as pd  
import csv

out_tsv = sys.argv[1]
tsv_total = sys.argv[2]
class_total = sys.argv[3]
tsv_no_tag_indel = sys.argv[4]
class_no_tag_indel = sys.argv[5]
tsv_no_tag_deletion = sys.argv[6]
class_no_tag_deletion = sys.argv[7]
tsv_no_tag_insertion = sys.argv[8]
class_no_tag_insertion = sys.argv[9]
tsv_no_tag_complex = sys.argv[10]
class_no_tag_complex = sys.argv[11]
tsv_single_tag = sys.argv[12]
class_single_tag = sys.argv[13]
tsv_multiple_tag = sys.argv[14]
class_multiple_tag = sys.argv[15]
tsv_any_tag = sys.argv[16]
class_any_tag = sys.argv[17]
tsv_single_tag_precise_tag = sys.argv[18]
class_single_tag_precise_tag = sys.argv[19]
tsv_single_tag_5indel = sys.argv[20]
class_single_tag_5indel = sys.argv[21]
tsv_single_tag_3indel = sys.argv[22]
class_single_tag_3indel = sys.argv[23]
tsv_single_tag_anyindel = sys.argv[24]
class_single_tag_anyindel = sys.argv[25]


#%%
df_total = pd.read_csv(tsv_total, sep="\t", header = None)
df_total.columns = ["sample_id", class_total]
df_no_tag_indel = pd.read_csv(tsv_no_tag_indel, sep="\t", header = None)
df_no_tag_indel.columns = ["sample_id", class_no_tag_indel]
df_no_tag_deletion = pd.read_csv(tsv_no_tag_deletion, sep="\t", header = None)
df_no_tag_deletion.columns = ["sample_id", class_no_tag_deletion]
df_no_tag_insertion = pd.read_csv(tsv_no_tag_insertion, sep="\t", header = None)
df_no_tag_insertion.columns = ["sample_id", class_no_tag_insertion]
df_no_tag_complex = pd.read_csv(tsv_no_tag_complex, sep="\t", header = None)
df_no_tag_complex.columns = ["sample_id", class_no_tag_complex]
df_single_tag = pd.read_csv(tsv_single_tag, sep="\t", header = None)
df_single_tag.columns = ["sample_id", class_single_tag]
df_multiple_tag = pd.read_csv(tsv_multiple_tag, sep="\t", header = None)
df_multiple_tag.columns = ["sample_id", class_multiple_tag]
df_any_tag = pd.read_csv(tsv_any_tag, sep="\t", header = None)
df_any_tag.columns = ["sample_id", class_any_tag]
df_single_tag_precise_tag = pd.read_csv(tsv_single_tag_precise_tag, sep="\t", header = None)
df_single_tag_precise_tag.columns = ["sample_id", class_single_tag_precise_tag]
df_single_tag_5indel = pd.read_csv(tsv_single_tag_5indel, sep="\t", header = None)
df_single_tag_5indel.columns = ["sample_id", class_single_tag_5indel]
df_single_tag_3indel = pd.read_csv(tsv_single_tag_3indel, sep="\t", header = None)
df_single_tag_3indel.columns = ["sample_id", class_single_tag_3indel]
df_single_tag_anyindel = pd.read_csv(tsv_single_tag_anyindel, sep="\t", header = None)
df_single_tag_anyindel.columns = ["sample_id", class_single_tag_anyindel]

dfs = [df_total, df_no_tag_indel, df_no_tag_deletion, df_no_tag_insertion, df_no_tag_complex, df_single_tag, df_multiple_tag, df_any_tag, df_single_tag_precise_tag, df_single_tag_5indel, df_single_tag_3indel, df_single_tag_anyindel]

df_merge = dfs[0]
for df in dfs[1:]:
    df_merge = pd.merge(df_merge, df, on = "sample_id", how = "outer")

#%% 
df_merge[class_no_tag_indel + "_ratio"] = round(df_merge[class_no_tag_indel] / df_merge[class_total], 2)
df_merge[class_no_tag_deletion + "_ratio"] = round(df_merge[class_no_tag_deletion] / df_merge[class_total], 2)
df_merge[class_no_tag_insertion + "_ratio"] = round(df_merge[class_no_tag_insertion] / df_merge[class_total], 2)
df_merge[class_no_tag_complex + "_ratio"] = round(df_merge[class_no_tag_complex] / df_merge[class_total], 2)
df_merge[class_single_tag + "_ratio"] = round(df_merge[class_single_tag] / df_merge[class_total], 2)
df_merge[class_multiple_tag + "_ratio"] = round(df_merge[class_multiple_tag] / df_merge[class_total], 2)
df_merge[class_any_tag + "_ratio"] = round(df_merge[class_any_tag] / df_merge[class_total], 2)
df_merge[class_single_tag_precise_tag + "_ratio"] = round(df_merge[class_single_tag_precise_tag] / df_merge[class_total], 2)
df_merge[class_single_tag_5indel + "_ratio"] = round(df_merge[class_single_tag_5indel] / df_merge[class_total], 2)
df_merge[class_single_tag_3indel+ "_ratio"] = round(df_merge[class_single_tag_3indel] / df_merge[class_total], 2)
df_merge[class_single_tag_anyindel + "_ratio"] = round(df_merge[class_single_tag_anyindel] / df_merge[class_total], 2)

df_merge.to_csv(out_tsv, sep = "\t", index = None, quoting = csv.QUOTE_NONE, quotechar = "",  escapechar = "\t")
# %%
