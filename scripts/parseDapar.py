#!/usr/bin/env python
import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='Parse Dapar Result')
parser.add_argument('--input', '-i', type=str, required=True, help='dapar results')
parser.add_argument('--config','-c',type=str,required=True,help='dapar configuration')
parser.add_argument('--long','-l',type=str,required=True,help='long isoform expression')
parser.add_argument('--short','-s',type=str,required=True,help="short isoform expression")
parser.add_argument('--pdui','-p',type=str,required=True,help="PDUI matrix")
args = parser.parse_args()


def sampleIdfromString(string):
    return [ substring.split("/")[-1].split(".")[0] for substring in string.split(",")]

f = open(args.config)

print("Parse config file ...")
for line in f:
    key,value = line.strip().split("=")
    if key == "Group1_Tophat_aligned_Wig":
        group1_ids = sampleIdfromString(value)
    elif key == "Group2_Tophat_aligned_Wig":
        group2_ids = sampleIdfromString(value)
print("Done .")


print("Parse dapar results ...")
df = pd.read_csv(args.input,sep="\t",index_col=0)
longAIndex = ["A_{}_long_exp".format(str(i+1)) for i in range(len(group1_ids))]
longBIndex = ["B_{}_long_exp".format(str(i+1)) for i in range(len(group2_ids))]
shortAIndex = ["A_{}_short_exp".format(str(i+1)) for i in range(len(group1_ids))]
shortBIndex = ["B_{}_short_exp".format(str(i+1)) for i in range(len(group2_ids))]
PDUIAIndex = ["A_{}_PDUI".format(str(i+1)) for i in range(len(group1_ids))]
PDUIBIndex = ["B_{}_PDUI".format(str(i+1)) for i in range(len(group2_ids))]
print("Done .")

print("Write results ...")
longMat = df.loc[:,list(longAIndex)+list(longBIndex)]
longMat.columns = group1_ids + group2_ids
longMat.to_csv(args.long,sep="\t")

shortMat = df.loc[:,list(shortAIndex)+list(shortBIndex)]
shortMat.columns =  group1_ids + group2_ids
shortMat.to_csv(args.short,sep="\t")

PDUIMat = df.loc[:,list(PDUIAIndex)+list(PDUIBIndex)]
PDUIMat.columns =  group1_ids + group2_ids
PDUIMat.to_csv(args.pdui,sep="\t")
print("Done .")
