#!/usr/bin/env python
import argparse
import subprocess
import os
import math
import logging 
import sys

logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('build index')

def star(args):
    fai = args.fasta + ".fai" 
    if not os.path.exists(fai):
        logger.info("No fai file  detected, build one ...")
        subprocess.run(["samtools","faidx",args.fasta])
        logger.info("Done .")
    L = 0
    n = 0
    with open(fai) as f:
        for line in f:
            L += int(line.strip().split("\t")[1])   
            n += 1
    genomeSAindexNbases = min(14,int(math.log2(L)/2-1))
    avgL = max(int(L/n),100)
    genomeChrBinNbits = min(18,int(math.log2(avgL)))
    if not os.path.exists(args.prefix):
        os.makedirs(args.prefix)
    if  os.path.exists(args.tmp_dir):
        logger.info("The specified tmp dir already exists .")
        sys.exit(1)        
    tmp_pdir = os.path.dirname(args.tmp_dir)
    if not os.path.exists(tmp_pdir):
        if len(tmp_pdir) > 0:
            os.makedirs(os.path.dirname(tmp_pdir))
    if args.gtf is not None:
        logger.info("Genome annotation provided in gtf format .")
        logger.info("Build star sequence index ...")
        cmd = ["STAR","--runMode","genomeGenerate","--genomeDir",args.prefix,"--genomeFastaFiles",args.fasta,"--genomeSAindexNbases",str(genomeSAindexNbases),"--genomeChrBinNbits",str(genomeChrBinNbits),"--sjdbGTFfile",args.gtf,"--runThreadN",str(args.threads),"--outTmpDir",args.tmp_dir]
        logger.info("Run: "+" ".join(cmd))
        subprocess.run(cmd)
    else:
        logger.info("Genome annotation not provided .")
        logger.info("Build star sequence index ...")
        cmd = ["STAR","--runMode","genomeGenerate","--genomeDir",args.prefix,"--genomeFastaFiles",args.fasta,"--genomeSAindexNbases",str(genomeSAindexNbases),"--genomeChrBinNbits",str(genomeChrBinNbits),"--runThreadN",str(args.threads),"--outTmpDir",args.tmp_dir]
        logger.info("Run: "+" ".join(cmd))
        subprocess.run(cmd)
    logger.info("Done .")
        

def bowtie2(args):
    if not os.path.exists(os.path.dirname(args.prefix)):
        os.makedirs(os.path.dirname(args.prefix))
    cmd = ["bowtie2-build","--threads",str(args.threads),args.fasta,args.prefix]
    logger.info("Build bowtie2 sequence index ...")
    logger.info("Run: "+" ".join(cmd))
    subprocess.run(cmd)
    logger.info("Done .")

def main():
    parser = argparse.ArgumentParser(description='Build index for STAR/bowtie2')
    parser.add_argument('--fasta','-f',required=True, help="Path of input fasta file of the reference sequence")
    parser.add_argument('--aligner','-a',required=True,choices=["star","bowtie2"],help="Aligner for index building")
    parser.add_argument('--gtf','-g',help="Path of genome annotation, only work for STAR genome index building")
    parser.add_argument('--prefix', '-o', required=True, help = "Prefix of the index")
    parser.add_argument('--threads','-p',type=int,default=1,help="Threads used for index building")
    parser.add_argument('--tmp-dir','-t',default="tmp",help="Tmp-dir for build STAR index")
    args = parser.parse_args()
    if args.aligner == "star":
        star(args)
    elif args.aligner == "bowtie2":
        bowtie2(args)

if __name__ == "__main__":
    main()


