#!/usr/bin/env python

"""
Given a list of fastq files, count the number of reads per file.

Usage:
    python scripts/count_fastq.py $out_tsv $fq_total $fq_indel $fq_del $fq_insertion $fq_tag_single $fq_tag_multiple $fq_tag_single_precise_insert $fq_tag_single_any_del $fq_tag_single_3del $fq_tag_single_5del sample_name
"""

#%%
import sys 
from Bio import SeqIO 
import os
import gzip

fastq = sys.argv[1]
out_tsv = sys.argv[2]
sample_name = sys.argv[3]

def get_read_count(fastq_gz):
    c = 0
    with gzip.open(fastq_gz, "rt") as fastq:
        for record in SeqIO.parse(fastq, "fastq"):
            c += 1
    return c

c = get_read_count(fastq)

os.makedirs(os.path.dirname(out_tsv), exist_ok=True)
with open(out_tsv, "w") as f:
    f.write(sample_name + "\t" + str(c))
    f.write("\n")