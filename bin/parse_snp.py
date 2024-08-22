#!/usr/bin/env python

"""
Generate a table summarizing wt PAM and mut PAM for precise insert reads.

Usage:
    python scripts/parse_snp_v2.py $bam $fastq 999 A G (1-based SNP location, WT_SNP, MUT_SNP) sample 
"""

#%%
from genericpath import exists
import gzip
from Bio import SeqIO
import numpy as np 
import os
import pandas as pd
import sys
import pysam

#%%
bam = sys.argv[1]
fastq = sys.argv[2]
snp_pos = int(sys.argv[3]) - 1 # pileup uses 0-based coordinates
snp_wt = sys.argv[4]
snp_mut = sys.argv[5]
out_wt = sys.argv[6]
out_mut = sys.argv[7]
out_other = sys.argv[8]
out_tsv = sys.argv[9]
sample = sys.argv[10]

#%%
# bam = "/Users/kaihu/Projects/Nathan/work/sikiproject/new_strategy/h3f3d_umi_correct_cutoff_5/workflow/res/bam/with_tag_single/SpCas9_gspPCR_SIKI2_1_1.bam"
# snp_pos = 1116 # A vs G 1-based coordinate

precise_insert_read_id = []
with gzip.open(fastq, "rt") as f:
    for record in SeqIO.parse(f, "fastq"):
        precise_insert_read_id.append(record.id)

c_wt = 0
c_mut = 0
c_other = 0
c = len(precise_insert_read_id)

os.makedirs(os.path.dirname(out_wt), exist_ok = True)
os.makedirs(os.path.dirname(out_mut), exist_ok = True)
os.makedirs(os.path.dirname(out_tsv), exist_ok = True)
with pysam.AlignmentFile(bam, "rb") as samfile, gzip.open(out_wt, "wt") as out_wt, gzip.open(out_mut, "wt") as out_mut, gzip.open(out_other, "wt") as out_other:
    for pileupcolumn in samfile.pileup():
        if pileupcolumn.pos == snp_pos:
            tem = 0
            for pileupread in pileupcolumn.pileups:
                # if sample == "SpCas9_gspPCR_SIKI2_1_9":
                #     print(pileupcolumn.pos, pileupcolumn.n)
                    # print(pileupread.alignment.query_name, pileupread)
                    # if pileupread.alignment.query_name == "m54328U_231110_222942/6489078/ccs_GAAAATAG":
                    #     print(pileupread.alignment.query_sequence[pileupread.query_position].upper())
                    # if pileupread.alignment.query_name == "m54328U_231110_222942/5441967/ccs_TTACGGAA":
                    #     print(pileupread.alignment.query_sequence[pileupread.query_position].upper())
                    # if pileupread.alignment.query_name in precise_insert_read_id:
                    #     print(pileupread.alignment.query_name,
                    #           pileupread.alignment.query_sequence[pileupread.query_position].upper(),
                    #           pileupread.alignment.is_supplementary)
                    if not pileupread.alignment.is_supplementary: # only use primary read
                        if pileupread.alignment.query_name in precise_insert_read_id:
                            tem += 1
                            if not pileupread.is_del and not pileupread.is_refskip:
                                snp = pileupread.alignment.query_sequence[pileupread.query_position].upper()
                                if snp == snp_wt:
                                    c_wt += 1
                                    SeqIO.write(record, out_wt, "fastq")
                                elif snp == snp_mut:
                                    c_mut += 1
                                    SeqIO.write(record, out_mut, "fastq")
                                else:
                                    c_other += 1
                                    SeqIO.write(record, out_other, "fastq")
            # print(sample, c, tem)
if c != 0:
    with open(out_tsv, "w") as f:
        f.write("\t".join([sample, 
                        str(c), 
                        str(c_wt), 
                        str(round(c_wt/c, 2)),
                        str(c_mut), 
                        str(round(c_mut/c, 2)),
                        str(c_other), 
                        str(round(c_other/c, 2))
                        ]))
        f.write("\n")
else:
    with open(out_tsv, "w") as f:
        f.write("\t".join([sample, 
                        str(c), 
                        str(c_wt), 
                        str(0),
                        str(c_mut), 
                        str(0),
                        str(c_other), 
                        str(0)
                        ]))
        f.write("\n")