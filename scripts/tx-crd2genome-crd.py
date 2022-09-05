#!/usr/bin/env python
import argparse
import logging
import numpy as np
from tqdm import tqdm
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(name)s] %(message)s')
logger = logging.getLogger("tx crd to genome crd")


def main():
    parser = argparse.ArgumentParser(description='map transcript scores to genome coordinate')
    parser.add_argument('--input','-i',required=True,help="input bedgraph")
    parser.add_argument('--output','-o',required=True,help="output bedgraph")
    parser.add_argument('--bed','-b',required=True,help="transcript information in bed12 format")
    args = parser.parse_args()
    # 0           1       2             3                4      5         6      7      8       9             10             11
    #chr1	11868	14409	ENST00000456328.2	100	+	11868	11868	0	3	359,109,1189,	0,744,1352, 
    tx_info = {}
    logger.info("Load transcript model ...")
    with open(args.bed) as f:
        for line in f:
            fields = line.strip().split("\t")
            chrom, start, end, tx_id, strand, block_sizes, block_starts = fields[0], int(fields[1]), int(fields[2]), fields[3], fields[5], fields[10], fields[11]
            block_sizes = [int(i) for i in block_sizes.split(",") if len(i) > 0]
            block_starts = [int(i) for i in block_starts.split(",") if len(i) > 0]
            L = sum(block_sizes)
            tx_info[tx_id] = (chrom, start, end, strand, block_sizes, block_starts, L)
    tx_scores = {}
    scores = None
    scalers = None
    logger.info("load scores mapped to transcript coordinate ...")
    fout = open(args.output,"w")
    with open(args.input) as fin:
        last_tx_id = ""
        for line in fin:
            tx_id, start, end, score = line.strip().split("\t")
            start, end, score = int(start), int(end), float(score)
            if tx_id not in tx_info:
                continue
            if tx_id != last_tx_id:
                if len(last_tx_id)!= 0:
                    scores[scalers>0] = scores[scalers>0]/scalers[scalers>0]
                    tx_scores[last_tx_id] = np.round(scores,2)
                L = tx_info[tx_id][-1]
                scores = np.zeros(L)
                scalers = np.zeros(L)
            scores[start:end] += score
            scalers[start:end] += 1
            last_tx_id = tx_id
    if scores is not None and scalers is not None and tx_id in tx_info:
        scores[scalers>0] = scores[scalers>0]/scores[scalers>0]
        tx_scores[tx_id] = scores.copy()
    logger.info("map scores to genome coordinate ...")
    for tx_id in tqdm(tx_scores):
        chrom, start, end, strand, block_sizes, block_starts, L = tx_info[tx_id]
        t_pos = 0
        scores = tx_scores[tx_id]
        g_pos = start
        if strand == "-":
            scores = scores[::-1]
        for block_size, block_start in zip(block_sizes,block_starts):
            #print(block_size, block_start)
            subseq_scores = scores[t_pos:t_pos+block_size]
            mask = subseq_scores[:-1] != subseq_scores[1:]
            if mask.sum() == mask.shape[0]:
                start_offsets = [0]
                end_offsets = [block_size]
            else:
                start_offsets = np.where(np.concatenate([np.array([True]),mask]))[0]
                end_offsets = np.concatenate([start_offsets[1:],[block_size]])
            for start_offset, end_offset in zip(start_offsets,end_offsets):
                score = subseq_scores[start_offset]
                gstart = g_pos + block_start + start_offset
                gend = g_pos + block_start + end_offset
                print(chrom,gstart,gend,tx_id,score,strand,sep="\t",file=fout)
            t_pos += block_size
    fout.close()

if __name__ == "__main__":
    main()
