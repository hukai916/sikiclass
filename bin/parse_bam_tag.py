#!/usr/bin/env python

"""
Parse bam file by scanning occurence of query read on the referecne.

Usage:
    python parse_bam_tag.py xxx.bam xxx.fastq.gz out_with_no_tag.fastq out_with_tag.single.fastq out_with_tag.multiple.fastq
"""

#%%
from genericpath import exists
import sys 
from Bio import SeqIO 
import os
import pysam
import gzip

bam = sys.argv[1]
fastq = sys.argv[2]
out_no_tag = sys.argv[3]
out_single_tag = sys.argv[4]
out_multiple_tag = sys.argv[5]
out_tag = sys.argv[6]

# bam = "/Users/kaihu/Projects/Nathan/work/sikiproject/new_strategy/h3f3d_umi_correct_cutoff_5/workflow/res/bam/tag/SpCas9_gspPCR_SIKI2_1_1.sam"
# fastq = "/Users/kaihu/Projects/Nathan/work/sikiproject/new_strategy/h3f3d_umi_correct_cutoff_5/workflow/fastq/SpCas9_gspPCR_SIKI2_1_1.fastq.gz"
# out_no_tag = "/Users/kaihu/Projects/Nathan/work/sikiproject/new_strategy/h3f3d_umi_correct_cutoff_5/workflow/res/fastq/01_without_tag/SpCas9_gspPCR_SIKI2_1_1.fastq.gz"
# out_single_tag = "/Users/kaihu/Projects/Nathan/work/sikiproject/new_strategy/h3f3d_umi_correct_cutoff_5/workflow/res/fastq/01_with_tag/SpCas9_gspPCR_SIKI2_1_1.single.fastq.gz"
# out_multiple_tag = "/Users/kaihu/Projects/Nathan/work/sikiproject/new_strategy/h3f3d_umi_correct_cutoff_5/workflow/res/fastq/01_with_tag/SpCas9_gspPCR_SIKI2_1_1.multiple.fastq.gz"

def cigar_consume(cigar_tuple):
    """_calculate the query and reference consumption given cigar_tuple, return tuple (query_consumed, ref_consumed)_
    
    Ref: https://samtools.github.io/hts-specs/SAMv1.pdf
    
    M 0 alignment match (can be a sequence match or mismatch) yes yes
    I 1 insertion to the reference yes no
    D 2 deletion from the reference no yes
    N 3 skipped region from the reference no yes
    S 4 soft clipping (clipped sequences present in SEQ) yes no
    H 5 hard clipping (clipped sequences NOT present in SEQ) no no
    P 6 padding (silent deletion from padded reference) no no
    = 7 sequence match yes yes
    X 8 sequence mismatch yes yes

    Args:
        cigar_tuple (_tuple_): _CIGAR tuple_
    """
    
    if cigar_tuple[0] == 0:
        return((cigar_tuple[1], cigar_tuple[1]))
    elif cigar_tuple[0] == 1: 
        return((cigar_tuple[1], 0))
    elif cigar_tuple[0] == 2: 
        return((0, cigar_tuple[1]))
    elif cigar_tuple[0] == 3: 
        return((0, cigar_tuple[1]))
    elif cigar_tuple[0] == 4: 
        return((cigar_tuple[1], 0))
    elif cigar_tuple[0] == 5: 
        return((0, 0))
    elif cigar_tuple[0] == 6: 
        return((0, 0))
    elif cigar_tuple[0] == 7: 
        return((cigar_tuple[1], cigar_tuple[1]))
    elif cigar_tuple[0] == 8:
        return((cigar_tuple[1], cigar_tuple[1]))
    else:
        print("wrong CIGAR!")
        exit()

#%% 
read_dict = {}

with pysam.AlignmentFile(bam, 'r') as bam_file:
    for i, alignment in enumerate(bam_file):
        ref_pos = alignment.reference_start
        # tag_pos = alignment.query_alignment_start
        tag_pos = 0 # query should start with 0 as soft clips will consume query; ref_pos should not start from 0 as CIGAR must NOT start with D or N.
        read_name = alignment.reference_name
        if not read_name in read_dict:
            read_dict[read_name] = [0] * alignment.query_length 
        for i, v in enumerate(alignment.cigartuples):
            if v[0] in [0, 7, 8]:
                for k in range(tag_pos, tag_pos + v[1]):
                    read_dict[read_name][k] += 1         
            query_consume, ref_consume = cigar_consume(v)
            tag_pos += query_consume

tag_dict = {}
for key in read_dict:
    tag_dict[key] = max(read_dict[key])

# print(read_dict["m54328U_231110_222942/1116482/ccs_CTACTTAG"])
# print(tag_dict["m54328U_231110_222942/1116482/ccs_CTACTTAG"])

os.makedirs(os.path.dirname(out_no_tag), exist_ok=True)
os.makedirs(os.path.dirname(out_single_tag), exist_ok=True)
os.makedirs(os.path.dirname(out_multiple_tag), exist_ok=True)
os.makedirs(os.path.dirname(out_tag), exist_ok=True)
# print(out_no_tag, out_single_tag, out_multiple_tag)

with gzip.open(fastq, "rt") as handle, gzip.open(out_no_tag, "wt") as out_no_tag, gzip.open(out_single_tag, "wt") as out_single_tag, gzip.open(out_multiple_tag, "wt") as out_multiple_tag, gzip.open(out_tag, "wt") as out_tag:
    for record in SeqIO.parse(handle, "fastq"):
        if not record.id in tag_dict:
            SeqIO.write(record, out_no_tag, "fastq")
        elif tag_dict[record.id] == 0:
            SeqIO.write(record, out_no_tag, "fastq")
        elif tag_dict[record.id] == 1:
            SeqIO.write(record, out_single_tag, "fastq")
            SeqIO.write(record, out_tag, "fastq")
        else:
            SeqIO.write(record, out_multiple_tag, "fastq") 
            SeqIO.write(record, out_tag, "fastq")
# %%
