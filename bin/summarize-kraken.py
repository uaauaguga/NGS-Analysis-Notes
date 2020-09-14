#!/usr/bin/env python
import argparse
import pandas as pd
import numpy as np
import os
from tqdm import tqdm

levels = ["0","D","K","P","C","O","F","G","S"]
parser = argparse.ArgumentParser(description='summarize kraken2 report')
parser.add_argument('--indir', '-i', type=str, required=True, help='dir contains  input kraken2 report')
parser.add_argument('--level','-l',type=str,required=True,choices=levels,help='level of summarize')
parser.add_argument('--output','-o',type=str,required=True,help="Where to put the configuration file")
args = parser.parse_args()

"""
1. Percentage of fragments covered by the clade rooted at this taxon
2. Number of fragments covered by the clade rooted at this taxon
3. Number of fragments assigned directly to this taxon
4. A rank code, indicating (U)nclassified, (R)oot, (D)omain, (K)ingdom,
   (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies.
   Taxa that are not at any of these 10 ranks have a rank code that is
   formed by using the rank code of the closest ancestor rank with
   a number indicating the distance from that rank.  E.g., "G2" is a
   rank code indicating a taxon is between genus and species and the
   grandparent taxon is at the genus rank.
5. NCBI taxonomic ID number
6. Indented scientific name

"""

records = []
for report in tqdm(os.listdir(args.indir)):
    sample_id = report.split(".")[0]
    with open(args.indir+"/"+report) as f:
        if args.level == "0":
                unclassified = next(f).strip().split("\t")[1]
                classified = next(f).strip().split("\t")[1]
                records.append((sample_id,unclassified,"unclassified"))
                records.append((sample_id,classified,"classified"))
        else:
            for line in f:
                percent,count,_,l,ID,name = line.strip().split("\t")
                name = name.strip()
                if l.strip() == args.level:
                    records.append((sample_id,count,ID+"|"+name))
df = pd.DataFrame.from_records(records)
df.columns = ["sample_id","count","taxo"]
df = df.pivot(index="taxo",columns="sample_id",values="count")
df.to_csv(args.output,sep="\t")
            
