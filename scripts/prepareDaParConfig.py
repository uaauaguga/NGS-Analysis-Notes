#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
parser = argparse.ArgumentParser(description='Prepare Configuration File for DaPar')
parser.add_argument('--pos','-p',type=str,required=True,help='Positive class')
parser.add_argument('--neg','-n',type=str,required=True,help='Negative class')
parser.add_argument('--indir','-i',type=str,required=True,help="Dir contains input wig files")
parser.add_argument('--outdir','-o',type=str,required=True,help="Output dir for dapar results")
parser.add_argument('--config','-c',type=str,required=True,help="Where to put the configuration file")
parser.add_argument('--coverage',type=str,default="30",help="Minimal coverage")
args = parser.parse_args()

coverage = args.coverage

posIds = open(args.pos).read().strip().split("\n")
negIds = open(args.neg).read().strip().split("\n")

print("{} postive samples".format(len(posIds)))
print("{} negative samples".format(len(negIds)))

posWigs=[ args.indir + "/" + sample_id + ".wig" for sample_id in posIds]
negWigs=[ args.indir + "/" + sample_id + ".wig" for sample_id in negIds]

posWigStr = ",".join(posWigs)
negWigStr = ",".join(negWigs)

template="""Annotated_3UTR=genome/bed/gencode.v27.annotation.UTR3p.bed
Group1_Tophat_aligned_Wig={}
Group2_Tophat_aligned_Wig={}
Output_directory={}
Output_result_file=result
Num_least_in_group1=1
Num_least_in_group2=1
Coverage_cutoff={}
FDR_cutoff=0.05
PDUI_cutoff=0.5
Fold_change_cutoff=0.59"""

with open(args.config,"w") as f:
    f.write(template.format(posWigStr,negWigStr,args.outdir,coverage))
