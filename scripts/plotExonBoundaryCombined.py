#!/usr/bin/env python
from collections import defaultdict
import pysam
import argparse
from tqdm import tqdm
import HTSeq
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np



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
    parser = argparse.ArgumentParser(description='Combine Coverage of Multiple Sample and Generate Plot')
    parser.add_argument('--coverages','-cs',type=str,required=True,help="Paths of coverages of different samples, separated by comma")
    parser.add_argument('--pdf','-p',type=str,default=None,help="Output coverage plot")
    args = parser.parse_args()
    paths = args.coverages.split(",")
    n = len(paths)
    df = pd.read_csv(paths[0],sep="\t",index_col=0) 
    for path in paths[1:]:
        df += pd.read_csv(paths[0],sep="\t",index_col=0)
    df = df/n
    plotCoverage(df,args.pdf)
if __name__ == "__main__":
    main()
