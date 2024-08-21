#!/usr/bin/env python

"""
Classify reads with single tag into precise_insert, with_any_del, with_5_del, and with_3_del.

Usage:
    python classify_tag_single.py xxx.csv xxx.fastq.gz 1-based_tag_start_pos 1-based_tag_end_pos num_flanking_base_to_check out_precise out_any_del out_5del out_3del 
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
tag_start = int(sys.argv[3])
tag_end = int(sys.argv[4])
flanking = int(sys.argv[5])
out_precise = sys.argv[6]
out_5del = sys.argv[7]
out_3del = sys.argv[8]
out_any_del = sys.argv[9]

# csv = "/Users/kaihu/Projects/Nathan/work/sikiproject/new_strategy/h3f3d_umi_correct_cutoff_5/workflow/res/stats/indel_ref_pos/with_tag_single/SpCas9_gspPCR_SIKI2_1_1.csv"
# csv = "/Users/kaihu/Projects/Nathan/work/sikiproject/new_strategy/h3f3d_umi_correct_cutoff_5/workflow/res/stats/indel_ref_pos/with_tag_single/SpCas9_uninjected_SIKI2_3_10.csv"
# fastq = "/Users/kaihu/Projects/Nathan/work/sikiproject/new_strategy/h3f3d_umi_correct_cutoff_5/workflow/res/fastq/01_with_tag_single/SpCas9_uninjected_SIKI2_3_10.fastq.gz"
# fastq = "/Users/kaihu/Projects/Nathan/work/sikiproject/new_strategy/h3f3d_umi_correct_cutoff_5/workflow/res/fastq/01_with_tag_single/SpCas9_uninjected_SIKI2_1_10.fastq.gz"

# tag_start = 1003
# tag_end = 1089
# flanking = 60
# out_precise = "/Users/kaihu/Projects/Nathan/work/sikiproject/new_strategy/h3f3d_umi_correct_cutoff_5/workflow/res/fastq/02_with_tag_single_precise_insert/precise.fastq.gz"
# out_any_del = "/Users/kaihu/Projects/Nathan/work/sikiproject/new_strategy/h3f3d_umi_correct_cutoff_5/workflow/res/fastq/02_with_tag_any_del/any.fastq.g"
# out_5del = "/Users/kaihu/Projects/Nathan/work/sikiproject/new_strategy/h3f3d_umi_correct_cutoff_5/workflow/res/fastq/02_with_tag_5del/five_del.fastq.g"
# out_3del = "/Users/kaihu/Projects/Nathan/work/sikiproject/new_strategy/h3f3d_umi_correct_cutoff_5/workflow/res/fastq/02_with_tag_3del/three_del.fastq.g"

#%% 
df = pd.read_csv(csv)
# print("\n" + csv + "\n")
if len(df) == 0:
    df["del_range"] = None
else:    
    df.loc[:, "del_range"] = df.apply(lambda row: range(row["ref_pos"], row["ref_pos"] + abs(row["indel"])), axis = 1)

def group_has_overlap(group, check_range):
    # check for each read_id group to see if del_ranges have any overlap with check_range
    ol = 0
    for rng in group["del_range"]:
        if any(x in check_range for x in rng):
            ol += 1
    return ol

check_range_tag_plus_flank = range(tag_start - 1 - flanking, tag_end + flanking)
check_range_5_flank = range(tag_start - 1 - flanking, tag_start + int((tag_end - tag_start)/2))
check_range_3_flank = range(tag_end - int((tag_end - tag_start)/2) - 1, tag_end + flanking)

#%%
# precise insert:
df_precise_insert = df.groupby("read_id").apply(lambda group: group_has_overlap(group, check_range_tag_plus_flank), include_groups=False)
if len(df_precise_insert) == 0:
    df_precise_insert = pd.DataFrame(columns=['indel_count'])
else:
    df_precise_insert = df_precise_insert.to_frame()

df_precise_insert.columns = ["indel_count"]
df_precise_insert = df_precise_insert.reset_index().rename(columns = {"index": "read_id"})
precise_insert = df_precise_insert[df_precise_insert["indel_count"] == 0]["read_id"].tolist()

# 5_flank_del: 
df_5del = df[df["indel"] < 0].groupby("read_id").apply(lambda group: group_has_overlap(group, check_range_5_flank), include_groups=False)
if len(df_5del) == 0:
    df_5del = pd.DataFrame(columns=['indel_count'])
else:
    df_5del = df_5del.to_frame()

df_5del.columns = ["indel_count"]
df_5del = df_5del.reset_index().rename(columns = {"index": "read_id"})
five_del = df_5del[df_5del["indel_count"] != 0]["read_id"].tolist()

# 3_flank_del:
df_3del = df[df["indel"] < 0].groupby("read_id").apply(lambda group: group_has_overlap(group, check_range_3_flank), include_groups=False)
if len(df_3del) == 0:
    df_3del = pd.DataFrame(columns=['indel_count'])
else:
    df_3del = df_3del.to_frame()

df_3del.columns = ["indel_count"]
df_3del = df_3del.reset_index().rename(columns = {"index": "read_id"})
three_del = df_3del[df_3del["indel_count"] != 0]["read_id"].tolist()

# any_del:
df_any_del = df[df["indel"] < 0].groupby("read_id").apply(lambda group: group_has_overlap(group, range(0, 10 ** 100)), include_groups=False)
if len(df_any_del) == 0:
    df_any_del = pd.DataFrame(columns=['indel_count'])
else:
    df_any_del = df_any_del.to_frame()

df_any_del.columns = ["indel_count"]
df_any_del = df_any_del.reset_index().rename(columns = {"index": "read_id"})
anydel = df_any_del[df_any_del["indel_count"] != 0]["read_id"].tolist()

#%% write to corresponding fastq files 
os.makedirs(os.path.dirname(out_precise), exist_ok=True)
os.makedirs(os.path.dirname(out_any_del), exist_ok=True)
os.makedirs(os.path.dirname(out_5del), exist_ok=True)
os.makedirs(os.path.dirname(out_3del), exist_ok=True)

with gzip.open(fastq, "rt") as handle, gzip.open(out_precise, "wt") as out_precise, gzip.open(out_any_del, "wt") as out_any_del, gzip.open(out_5del, "wt") as out_5del, gzip.open(out_3del, "wt") as out_3del:
    for record in SeqIO.parse(handle, "fastq"):
        if record.id in precise_insert:
            SeqIO.write(record, out_precise, "fastq")
        if record.id in anydel:
            SeqIO.write(record, out_any_del, "fastq")
        if record.id in five_del:
            SeqIO.write(record, out_5del, "fastq")
        if record.id in three_del:
            SeqIO.write(record, out_3del, "fastq")
        if not record.id in df.groupby("read_id").groups.keys():
            SeqIO.write(record, out_precise, "fastq")

# %%
