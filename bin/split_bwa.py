#!/usr/bin/env python

"""
Split reference into individual referneces and run BWA separately, then merge resulting sam into a single sam.

Usage:
python split_bwa.py ref.fasta ref/tag.fastq bam/tag/split/${filename}.bam

tag.fastq is the reads to map against each of ref.fasta

Dev notes:
1. to capture shorter supplementary hits, use the following options:
-a -k 5 -T 15 -D 0.1
"""

from genericpath import exists
import sys 
from Bio import SeqIO 
import os
import subprocess

ref = sys.argv[1]
read = sys.argv[2]
outfile = sys.argv[3]

# read in ref seq:
ref_dict = SeqIO.to_dict(SeqIO.parse(ref, "fasta"))

# create result temp folder:
outdir = os.path.dirname(outfile)
if outdir == "":
    outtem_fa = outdir + "temp/fa/"
    outtem_bam = outdir + "temp/bam/"
else: 
    outtem_fa = outdir + "/temp/fa/"
    outtem_bam = outdir + "/temp/bam/"

os.makedirs(outtem_fa, exist_ok = True)
os.makedirs(outtem_bam, exist_ok = True)

# iterate through ref
for ref in ref_dict:
    print("processing ", ref, " ...")
    out_fa = outtem_fa + str(ref).replace("/", "_") + ".fa"
    SeqIO.write(ref_dict[ref], out_fa, format = "fasta")
    # bwa index
    bwa_index = "bwa index " + out_fa
    subprocess.call(bwa_index, shell = True, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)
    # print(bwa_index)
    
    # bwa mem
    out_bam = outtem_bam + str(ref).replace("/", "_") + ".bam"
    bwa_mem = "bwa mem -a -k 5 -T 15 -D 0.1 " + out_fa + " " + read + " | samtools view -o " + out_bam
    # print(bwa_mem)
    subprocess.call(bwa_mem, shell = True, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)

# merge bam files into a single bam
# bam_merge = "ulimit -n 99999999 && samtools merge -o " + outfile + ".tem.bam" + " " + outtem_bam + "*"
## Docker may not have permission to run ulimit, use -b instead
find_cmd = "find " + outtem_bam + " -name '*.bam' > " + outtem_bam + "bam_list.txt"
subprocess.call(find_cmd, shell = True, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)

bam_merge = "samtools merge -o " + outfile + ".tem.bam" + " -b " + outtem_bam + "bam_list.txt"
subprocess.call(bam_merge, shell = True, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)
print(bam_merge)

# filter out unmapped reads 
bam_filter = "samtools view -h -F 4 " + outfile + ".tem.bam" + " -o " + outfile + " && rm " + outfile + ".tem.bam"
subprocess.call(bam_filter, shell = True, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)

# remove temp folder:
if outdir == "":
    cmd_clean = "rm -r " + outdir + "temp"
else: 
    cmd_clean = "rm -r " + outdir + "/temp"
subprocess.call(cmd_clean, shell = True, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)