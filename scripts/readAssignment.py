#!/usr/bin/env python
from collections import defaultdict
import pysam
import pandas as pd
import argparse
from tqdm import tqdm
import HTSeq
from matplotlib import pyplot as plt

def main():
    parser = argparse.ArgumentParser(description='Assign reads to different genomic regions')
    parser.add_argument('--input','-i',type=str,required=True,help="Input bam file, should in hg38 coordinate, can either sorted by query name or coordiante")
    parser.add_argument('--strandness','-s',type=str,default="no",choices=["forward","reverse","no"])
    parser.add_argument('--output','-o',type=str,required=True,help="Output reads assignment")
    parser.add_argument('--beddir','-bd',type=str,default="genome/bed",help="Dir that contains bed files")
    parser.add_argument('--priority','-p',type=str,default="lncRNA,mRNA,snoRNA,snRNA,srpRNA,tRNA,tucpRNA,Y_RNA,pseudogene,exon,intron,antisense,promoter,enhancer,repeats")
    args = parser.parse_args()

    bam = HTSeq.BAM_Reader(args.input)
    regions = args.priority.strip().split(",")
    print("Load genomic regions ...")
    if args.strandness != "no":
        ga = HTSeq.GenomicArrayOfSets("auto", stranded=True)
    else:
        ga = HTSeq.GenomicArrayOfSets("auto", stranded=False)
    for region in regions:
        bed = args.beddir + "/" + region + ".bed"
        with open(bed,"r") as f:
            for line in f:
                line = line.strip()
                if line.startswith("#") or len(line) == 0:
                    continue
                fields = line.split("\t") 
                chrom,start,end,strand = fields[0],int(fields[1]),int(fields[2]),fields[5]
                iv = HTSeq.GenomicInterval(chrom,start,end,strand = strand) 
                ga[iv] += region
        print("{} loaded".format(region))
    print("Done .")
    stats = defaultdict(int)
    n_total_fragments = 0
    for read1,read2 in tqdm(HTSeq.pair_SAM_alignments_with_buffer(bam,max_buffer_size=5000000)):
        n_total_fragments += 1
        # ignore singletons
        if (read1 is None) or (read2 is None):
            stats['singleton'] += 1
            continue
        # ignore unmapped reads
        if not (read1.aligned and read2.aligned):
            stats['unmapped'] += 1
            continue
        if read1.iv.chrom != read2.iv.chrom:
            stats['diff_chrom'] += 1
            continue
        if args.strandness == 'forward':
            read2.iv.strand = read1.iv.strand
        elif args.strandness == 'reverse':
            read1.iv.strand = read2.iv.strand
        else:
            read1.iv.strand = "."
            read2.iv.strand = "."
        featureSet = set()
        for iv0, step_set in ga[read1.iv].steps():
            featureSet = featureSet.union(step_set)
        for iv0, step_set in ga[read2.iv].steps():
            featureSet = featureSet.union(step_set)
        for region in regions:
            if region in featureSet:
                stats[region] += 1
                break
    
    n_assigned = pd.Series(stats).sum()
    stats["unassigned"] = n_total_fragments - n_assigned
    stats["total"] = n_total_fragments
    with open(args.output,"w") as f:
        for region in regions:
            print(region,stats[region],sep="\t",file=f)
        for each in ['singleton','unmapped','diff_chrom','unassigned','total']:
            print(each,stats[each],sep="\t",file=f)
    

if __name__ == "__main__":
    main()
