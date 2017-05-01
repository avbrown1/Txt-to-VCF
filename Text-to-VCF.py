#!/usr/bin/env python3.5
# NAME
#    VCF-Conversion.py - convert Hapamp data to VCF w.r.t. Wm82.a2
#
# SYNOPSIS
#     Hapmap to VCF file
#
# INPUT FILES
#     SEQLEN_FILE
#         Sequence lengths from JGI Wm82.a2 assembly (two columns: SEQID, LENGTH). The "mtDNA" length is needed.
#
#     ASSEMBLY_REPORT
#         ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000004515.3_Glycine_max_v2.0/GCA_000004515.3_Glycine_max_v2.0_assembly_report.txt
#
#     NCBI_TO_JGI_COORDINATE_FILE
#         Four tab-separated columns:
#             SEQID    START    END    OFFSET
#         Each entry in this file indicates an interval in the NCBI Wm82.a2 assembly
#         to which OFFSET should be added to map to the corresponding JGI Wm82.a2 assembly.
# 
#     SNPBATCH_FILE
#         dbSNP SNPBATCH file (snpBatch_BARC_1059007.gz)
#
#     SOYSNP50K_FILE
#         Korean_snp_working1.txt
#
# ENVIRONMENT VARIABLES
#     SOYSNP50K_ONLY_HOMOZYGOUS
#         If set, only output samples that contain all homozygous alleles (no heterozygous or missing).
#         Useful for input to Beagle.
#
# CHANGE HISTORY
#     2016-10-27    Use 
#     2016-09-28    Use sequence lengths from the JGI assembly. Translate coordinates from NCBI- to JGI-assembly space.
#     2016-06-03    Initial version
#
# AUTHOR
#     Anne Brown < anne.brown@ars.usda.gov> edited from Nathan Weeks <nathan.weeks@ars.usda.gov>

import datetime
import gzip
import operator
import os
import re
import sys


print('##fileformat=VCFv4.2')
print('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
print('##fileDate=' + datetime.date.today().strftime("%Y%m%d"))

# Get (JGI) sequence lengths
with open("/scratch/abrown1/Gmax_275_v2.0.softmasked.fa.seqlen") as seqlen:
    seqid_to_length = {seqid: int(length) for (seqid, length) in (line.split() for line in seqlen)}

for seqid in seqid_to_length:
    print('##contig=<ID={},length={},assembly=Wm82.a2,species="Glycine max",taxonomy=3847>'.format(seqid, seqid_to_length[seqid]))

VCF = []

with open("/scratch/abrown1/Song_new_data/New_Cultivars_50k") as f: 
    print("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", 
    	   "FORMAT", *f.readline().split()[6:],sep="\t")
    for line in f:
        A = line[:-1].split("\t")
        CHROm = A[4]
        POS = A[5]
        ID = A[0]
        REF = A[2]
        ALT = A[3]
        print("\t".join((CHROM, POS, ID, REF, ALT, ".", ".", ".", "GT",
                        ## SNP data has 2 nucleotides, not 1 so we need to make it so it recognizes both nucleotides as being REF/REF, ALT/ALT, etc. 
        *('0/0' if allele == REF else ('1/1' if allele == ALT else ('./.' if allele == 'N' else '0/1')) 
        for allele in A[6:]))))