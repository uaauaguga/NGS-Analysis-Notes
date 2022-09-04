#!/usr/bin/env python
import argparse
import logging
import numpy as np
from tqdm import tqdm
import sys
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(name)s] %(message)s')
logger = logging.getLogger("bed to bedgraph")


def aggregate(records):
    gstart = records[0][0]
    gend = records[-1][1]
    L = gend - gstart
    try:
        scores = np.zeros(L)
        scalers = np.zeros(L)
    except:
        print(records)
        sys.exit(0)
    for start, end, score in records:
        start, end = start - gstart, end - gstart
        scores[start:end] += score
        scalers[start:end] += 1
    scores[scalers>0] = scores[scalers>0]/scalers[scalers>0]
    scores = np.round(scores,2)
    mask = np.concatenate([[True],scores[0:-1] != scores[1:]])
    start_offsets = np.where(mask)[0]
    end_offsets = np.concatenate([start_offsets[1:],[L]])
    records = []
    for start_offset, end_offset in zip(start_offsets,end_offsets):
        score = scores[start_offset]
        if score == 0:
            continue
        records.append((gstart+start_offset,gstart+end_offset,score))
    return records

def main():
    parser = argparse.ArgumentParser(description='combine score defined column 5 of bed file')
    parser.add_argument('--input','-i',required=True,help="input bed, should contain column 5")
    parser.add_argument('--output','-o',required=True,help="output bedgraph")
    args = parser.parse_args()

    
    fout = open(args.output,"w")
    last_chrom = ""
    cache = []
    with open(args.input) as fin:
        for line in tqdm(fin):
            fields = line.strip().split("\t")
            chrom, start, end, score = fields[0], int(fields[1]), int(fields[2]), float(fields[4])
            update = False
            if chrom != last_chrom:
                # processing a new chromsome
                last_end = -1
            if start > last_end:
                update = True
                # current interval does not overlap with the last one
            if update and len(cache) > 0:
                for start_, end_, score_ in aggregate(cache):
                    print(last_chrom, start_, end_, score_, file=fout,sep="\t")
                cache = []
            cache.append((start, end, score))
            last_chrom, last_start, last_end, last_score = chrom, start, end, score
    if len(cache) > 0:
        for start_, end_, score_ in aggregate(cache):
            print(start_, end_, score_, file=fout,sep="\t")
    fout.close()


if __name__ ==  "__main__":
    main()







             
