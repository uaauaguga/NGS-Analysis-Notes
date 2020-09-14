#!/usr/bin/env python
import argparse
import pysam

parser = argparse.ArgumentParser(description='Get Reads Aligned to Viral Decoy')
parser.add_argument('--inbam', '-ib', type=str, required=True, help='input bam file')
parser.add_argument('--inbai','-ii',type=str,help="input index file")
parser.add_argument('--virus','-v',type=str,default="ref-data/decoy-virus-name.txt")
parser.add_argument('--outbam','-ob',type=str,help="Output viral decoy aligned reads")
args = parser.parse_args()

inBam = pysam.AlignmentFile(args.inbam,"rb",filepath_index=args.inbai)
outBam = pysam.AlignmentFile(args.outbam, "wb", template=inBam)
viruses = open(args.virus).read().strip().split("\n")

for virus in viruses:
    print("Processing {} ...".format(virus))
    for read in inBam.fetch(contig=virus):
        outBam.write(read)

inBam.close()
outBam.close()
        
        



