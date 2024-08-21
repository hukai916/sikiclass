#!/usr/bin/env python

"""
Classify reads without tag into indel, del, insertion.

Usage:
    python scripts/classify_without_tag.py $csv $fastq 1113 $out_indel $out_del $out_insertion
"""

#%%
import sys 
from Bio import SeqIO 
import os
import pandas as pd
import gzip
import numpy as np

csv = sys.argv[1]
fastq = sys.argv[2]
pam = sys.argv[3] # 1-based PAM start site
out_indel = sys.argv[4]
out_del = sys.argv[5]
out_insertion = sys.argv[6]

# csv = "/Users/kaihu/Projects/Nathan/work/sikiproject/new_strategy/h3f3d_umi_correct_cutoff_5/workflow/res/stats/indel_ref_pos/without_tag/SpCas9_gspPCR_SIKI2_1_1.csv"

#%% 
df = pd.read_csv(csv)
df_del = df[df["indel"] < 0]
df_insertion = df[df["indel"] > 0]

#%%
# with indel:
df_indel = df.groupby("read_id").apply(lambda group: len(group))
if len(df_indel) == 0:
    df_indel = pd.DataFrame(columns=['indel_count'])
else:
    df_indel = df_indel.to_frame()
    df_indel.columns = ["indel_count"]
df_indel = df_indel.reset_index().rename(columns = {"index": "read_id"})
indel = df_indel[df_indel["indel_count"] > 0]["read_id"].tolist()

# with del: 
df_del = df_del.groupby("read_id").apply(lambda group: len(group))
if len(df_del) == 0:
    df_del = pd.DataFrame(columns=['indel_count'])
else:
    df_del = df_del.to_frame()
    df_del.columns = ["indel_count"]
df_del = df_del.reset_index().rename(columns = {"index": "read_id"})
deletion = df_del[df_del["indel_count"] > 0]["read_id"].tolist()

# with insertion: 
df_insertion = df_insertion.groupby("read_id").apply(lambda group: len(group))
if len(df_insertion) == 0:
    df_insertion = pd.DataFrame(columns=['indel_count'])
else:
    df_insertion = df_insertion.to_frame()
    df_insertion.columns = ["indel_count"]
df_insertion = df_insertion.reset_index().rename(columns = {"index": "read_id"})
insertion = df_insertion[df_insertion["indel_count"] > 0]["read_id"].tolist()


#%% write to corresponding fastq files 
os.makedirs(os.path.dirname(out_indel), exist_ok=True)
os.makedirs(os.path.dirname(out_del), exist_ok=True)
os.makedirs(os.path.dirname(out_insertion), exist_ok=True)

with gzip.open(fastq, "rt") as handle, gzip.open(out_indel, "wt") as out_indel, gzip.open(out_del, "wt") as out_del, gzip.open(out_insertion, "wt") as out_insertion:
    for record in SeqIO.parse(handle, "fastq"):
        if record.id in indel:
            SeqIO.write(record, out_indel, "fastq")
        if record.id in deletion:
            SeqIO.write(record, out_del, "fastq")
        if record.id in insertion:
            SeqIO.write(record, out_insertion, "fastq")
# %%
