#!/usr/bin/env python

"""
Given bam and ref file, parse out the locations of each indel occurences along the reference.

Usage:
python parse_bam_indel.py xxx.bam ref.fa outfile.csv    
"""
# %%
import sys 
import pysam
from Bio import SeqIO
import os

# bam = "../stats/summary_counts/SpCas9_uninjected_SIKI2_1_10/bwa/fastq_no_insert_loose.bam"
# ref = "../ref/h3f3d_wt.fa"
# outfile = "../test/test.csv"
bam = sys.argv[1]
ref = sys.argv[2]
outfile = sys.argv[3]

# %%
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

# %%
# Read in ref fasta file:
ref_record = [x for x in SeqIO.parse(ref, "fasta")]
assert len(ref_record) == 1, "More than one record in ref!"
ref_fa = ref_record[0].seq

# %%
# For each record in BAM, parse all I/D occurrences for each matched region
dict = {i: [] for i in range(len(ref_fa))} 
    # 0: [(read_id, I/D_length), (...)] # positives represent insertions; negatives represent deletions.
d_clip = {} # stores inferred deletions from clipped alignment regions

with pysam.AlignmentFile(bam, 'rb') as bam_file:
    len_ref = bam_file.get_reference_length(bam_file.references[0])
    for i, alignment in enumerate(bam_file):
        ref_pos  = alignment.reference_start
        ref_end  = alignment.reference_end
        if ref_pos == -1:
            continue
        
        if not ref_pos == 0 or not ref_end == len_ref:
            if not alignment.query_name in d_clip:
                d_clip[alignment.query_name] = [range(ref_pos, ref_end)]
            else:
                d_clip[alignment.query_name].append(range(ref_pos, ref_end))
                
        # if alignment.query_name == "m54328U_231110_222942/6489078/ccs_GAAAATAG":
        #     # covered_ref.append(range(alignment.reference_start, alignment.reference_end))
        #     print(len_ref, alignment.reference_start, alignment.reference_end)
        
        for i, v in enumerate(alignment.cigartuples):
            if v[0] == 1:
                indel = (alignment.query_name, v[1])
                dict[ref_pos].append(indel)
            elif v[0] == 2:
                indel = (alignment.query_name, -v[1])
                dict[ref_pos].append(indel)          
            query_consume, ref_consume = cigar_consume(v)
            ref_pos = ref_pos + ref_consume

# for each query, infer the deletions from clips and add them to dict
def merge_ranges(ranges):
    # Sort the ranges by their start value
    sorted_ranges = sorted(ranges, key=lambda r: r.start)
    
    merged_ranges = []
    for r in sorted_ranges:
        if not merged_ranges or merged_ranges[-1].stop < r.start:
            # No overlap, add a new range
            merged_ranges.append(r)
        else:
            # Overlap detected, merge with the last range
            merged_ranges[-1] = range(merged_ranges[-1].start, max(merged_ranges[-1].stop, r.stop))
    
    return merged_ranges

def find_gaps(merged_ranges):
    gaps = []
    for i in range(len(merged_ranges) - 1):
        end_current_range = merged_ranges[i].stop
        start_next_range = merged_ranges[i + 1].start
        
        if end_current_range < start_next_range:
            gap = range(end_current_range, start_next_range)
            gaps.append(gap)
    
    return gaps

for k in d_clip:
    merged_ranges = merge_ranges(d_clip[k])
    gaps = find_gaps(merged_ranges)
    for gap in gaps:
        dict[gap.start].append((k, gap.start - gap.stop + 1)) # pysam coordinates follow [)] rule.

os.makedirs(os.path.dirname(outfile), exist_ok=True)
with open(outfile, "w") as f:
    f.write("ref_pos\tref_base\tindel\tread_id\n")
    for pos in dict:
        if not len(dict[pos]) == 0:
            for hit in dict[pos]:
                f.write("\t".join([str(pos), ref_fa[pos], str(hit[1]), str(hit[0])]))
                f.write("\n")
                # print(",".join([str(pos), ref_fa[pos], str(hit[1]), str(hit[0])]))
                # print("\n")