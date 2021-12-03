#!/usr/bin/env python
from collections import defaultdict
import pysam
import argparse
from tqdm import tqdm
import HTSeq
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np


def checkStrandness(iv,mate,strandness):
    flipStrand = {"+":"-","-":"+"}
    if (mate == "1" and strandness == "reverse") or (mate == "2" and strandness == "forward"):
        iv.strand = flipStrand[iv.strand]
    if strandness == "no":
        iv.strand = "."
    return iv


def plotCoverage(df,path):
    fig,axes = plt.subplots(2,2,figsize=(7,7),sharey=True) #,sharex=True, 
    coordinate = np.arange(-50,50)
    axes[0,0].plot(coordinate,df["read1-5p"].values,color="black",lw=2.5)
    axes[0,0].set_title("Read1-Exon-5-Prime",fontsize=13,weight="bold")
    axes[0,0].set_ylim([0,df["read1-5p"].max()*1.02])
    axes[0,1].plot(coordinate,df["read1-3p"].values,color="black",lw=2.5)
    axes[0,1].set_title("Read1-Exon-3-Prime",fontsize=13,weight="bold")
    axes[1,0].plot(coordinate,df["read2-5p"].values,color="black",lw=2.5)
    axes[1,0].set_title("Read2-Exon-5-Prime",fontsize=13,weight="bold")
    axes[1,1].plot(coordinate,df["read2-3p"].values,color="black",lw=2.5)
    axes[1,1].set_title("Read2-Exon-3-Prime",fontsize=13,weight="bold")
    for i in [0,1]:
        for j in [0,1]:
            axes[i,j].tick_params(labelsize=13)
    _ = fig.text(0.5, 0.05, 'Position Relative to Boundary', ha='center',fontsize=17,weight="bold")
    _ = fig.text(0, 0.5, 'Average Coverage', va='center', rotation='vertical',fontsize=17,weight="bold")
    plt.savefig(path,bbox_inches="tight")


def main():
    parser = argparse.ArgumentParser(description='Assign reads to different genomic regions')
    parser.add_argument('--input','-i',type=str,required=True,help="Input bam file, should in hg38 coordinate")
    parser.add_argument('--strandness','-s',type=str,default="no",choices=["forward","reverse","no"])
    parser.add_argument('--filter','-f',type=int,default=3,help="Only consider exon with mean coverage higher than this value")
    parser.add_argument('--gtf','-a',type=str,default="genome/gtf/gencode.v27.annotation.gtf",help="gtf annotation")
    parser.add_argument('--coverage','-c',type=str,required=True,help="Output coverage")
    parser.add_argument('--pdf','-p',type=str,default=None,help="Output coverage plot")
    args = parser.parse_args()
    
    if args.strandness != "no":
        ga1 = HTSeq.GenomicArray("auto", stranded=True)
        ga2 = HTSeq.GenomicArray("auto", stranded=True)
    else:
        ga1 = HTSeq.GenomicArray("auto", stranded=False)
        ga2 = HTSeq.GenomicArray("auto", stranded=False)

    #chr1	HAVANA	exon	12613	12721 
           
             
    print("Load bam file ...")            
    bam = HTSeq.BAM_Reader(args.input)
    for read1,read2 in tqdm(HTSeq.pair_SAM_alignments_with_buffer(bam,max_buffer_size=5000000)):
        if (read1 is None) or (read2 is None):
            continue
        if not (read1.aligned and read2.aligned):
            continue
        if read1.iv.chrom != read2.iv.chrom:
            continue
        else:
            read1.iv.strand = "."
            read2.iv.strand = "."
        for cigop in read1.cigar:
            if cigop.type != "M":
                continue
            ga1[ checkStrandness(cigop.ref_iv,"1",args.strandness) ] += 1
        for cigop in read2.cigar:
            if cigop.type != "M":
                continue
            ga2[ checkStrandness(cigop.ref_iv,"2",args.strandness) ] += 1  
    print("Done.")
    

    exonNumber = 0
    print("Get coverage of exons in gtf annotation ...")
    fivePrime1 = np.zeros(100)
    fivePrime2 = np.zeros(100)
    threePrime1 = np.zeros(100)
    threePrime2 = np.zeros(100)
    with open(args.gtf) as f:
        for line in tqdm(f):
            line = line.strip()
            if len(line) == 0 or line.startswith("#"):
                continue
            fields = line.split("\t")
            if fields[2] != "exon":
                continue
            strand = fields[6]
            if strand == "+":
                fivePrimeBoundary = HTSeq.GenomicInterval(fields[0],int(fields[3])-1-50,int(fields[3])+49,strand)
                threePrimeBoundary = HTSeq.GenomicInterval(fields[0],int(fields[4])-1-50,int(fields[4])+49,strand)
                fivePrime1_ = np.fromiter(ga1[fivePrimeBoundary],dtype="i")
                fivePrime2_ = np.fromiter(ga2[fivePrimeBoundary],dtype="i")
                threePrime1_ = np.fromiter(ga1[threePrimeBoundary],dtype="i")
                threePrime2_ = np.fromiter(ga2[threePrimeBoundary],dtype="i")
            else:
                fivePrimeBoundary = HTSeq.GenomicInterval(fields[0],int(fields[4])-1-50,int(fields[4])+49,strand)
                threePrimeBoundary = HTSeq.GenomicInterval(fields[0],int(fields[3])-1-50,int(fields[3])+49,strand)
                fivePrime1_ = np.fromiter(ga1[fivePrimeBoundary],dtype="i")[::-1]
                fivePrime2_ = np.fromiter(ga2[fivePrimeBoundary],dtype="i")[::-1]
                threePrime1_ = np.fromiter(ga1[threePrimeBoundary],dtype="i")[::-1]
                threePrime2_ = np.fromiter(ga2[threePrimeBoundary],dtype="i")[::-1]
            if (fivePrime1_.mean() > args.filter) or (fivePrime2_.mean() > args.filter) or (threePrime1_.mean() > args.filter) or (threePrime2_.mean() > args.filter):
                exonNumber += 1
                fivePrime1 += fivePrime1_
                fivePrime2 += fivePrime2_
                threePrime1 += threePrime1_
                threePrime2 += threePrime2_
            
    print("Done .")
    df = pd.DataFrame({"read1-5p":fivePrime1,"read1-3p":threePrime1,"read2-5p":fivePrime2,"read2-3p":threePrime2})
    df = df/exonNumber
    df.to_csv(args.coverage,sep="\t")
    if args.pdf is not None:
        plotCoverage(df,args.pdf)
               

                
                            
    

if __name__ == "__main__":
    main()
