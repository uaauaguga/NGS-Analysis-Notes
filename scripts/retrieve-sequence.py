#!/usr/bin/env python
from Bio import Entrez
import sys
import time
import argparse

parser = argparse.ArgumentParser(description='Download genebank sequence')
parser.add_argument('--query', '-q',help="Genebank id for download",required=True)
parser.add_argument('--fasta','-o',help="Downloaded sequences",required=True)
args = parser.parse_args()

Entrez.email = "jinyf16@mails.tsinghua.edu.cn"
gbIds = open(args.query).read().strip().split("\n")
print("{} query sequence".format(len(gbIds)))

fout = open(args.fasta,"a")

for gbId in gbIds:
    try:
        print("Start retriving {} from ncbi nucleotide...".format(gbId),file=sys.stderr)
        handle = Entrez.efetch(db="nucleotide", id=gbId, rettype="fasta", retmode="text")
        content = handle.read().strip()
        contents = content.split("\n")
        entry = contents[0] + "\n" + "".join(contents[1:])
        print(entry,file=fout)
        print("Done.")
    except:
        print("Error retriving {}, skip ...".format(gbId),file=sys.stderr)
    time.sleep(0.5)
fout.close()
