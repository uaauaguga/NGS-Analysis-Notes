#!/usr/bin/env python
import argparse
from scipy import stats
import gzip
import sys
import numpy as np
from collections import defaultdict
from tqdm import tqdm

def get_kmer_count(sequence, k = 3):
    counts = defaultdict(int)
    for i in range(len(sequence)-k):
        if b"N" in sequence[i:i+k]:
            continue
        counts[sequence[i:i+k]] += 1
    return counts

def get_N_count(sequence):
    n_N = 0
    for c in sequence:
        if c == 78:
            n_N += 1
    return n_N
        
def entropy(sequence, k=3):
    if len(sequence) < k:
        return 0
    n_N = get_N_count(sequence)
    if n_N >= 4:
        return 0
    counts = get_kmer_count(sequence, k=k)
    counts = np.array(list(counts.values()))
    ent = stats.entropy(counts)
    return 100*ent/entropy_scaler


def dust(sequence, k=3):
    if len(sequence) < k:
        return 100
    counts = get_kmer_count(sequence, k=k)
    dust_score = 0
    for count in counts.values():
        dust_score += count*(count-1)
    n_N = get_N_count(sequence)
    if n_N >= 4:
        return 100
    n = len(sequence) - n_N - 2 
    ## for polymer of one nucleotide, this score is
    ## n * (n - 1)
    return 100*dust_score/((n - 1)*n)


def main():
    parser = argparse.ArgumentParser(description='Analyze sequence complexity of fastq file')
    parser.add_argument('--in-fastq-1','-if1', help = "Pair 1 of input fastq file", required=True)
    parser.add_argument('--in-fastq-2','-if2', help = "Pair 2 of input fastq file", required=True)
    parser.add_argument('--method', '-m',  help="Method for sequence comlexicity calculation",default="dust",choices=["dust","entropy"])
    parser.add_argument('--word-size', '-k', help="k mer for complexity calculation",type=int,default=3)
    parser.add_argument('--out-fastq-1','-of1', help = "Pair 1 of output fastq file")
    parser.add_argument('--out-fastq-2','-of2', help = "Pair 2 of output fastq file")
    parser.add_argument('--filter','-f',type=float,help="If this paramter is setted, read pairs with entropy lower than this value or dust score higher then this value will be filtered in output fastq file")
    parser.add_argument('--statistics','-s',help="histogram of complexity distribution")
    args = parser.parse_args()
    input_1 =  gzip.open(args.in_fastq_1)
    input_2 =  gzip.open(args.in_fastq_2)
    count  = 0
    filtered = 0
    k = args.word_size

    if args.method == "entropy":
        global entropy_scaler
        entropy_scaler = stats.entropy(np.ones(4**args.word_size))

    if args.out_fastq_1 and args.out_fastq_2:
        output_1 = gzip.open(args.out_fastq_1,'w')
        output_2 = gzip.open(args.out_fastq_2,'w')

    scores = np.zeros(101)

    for one_1,one_2 in tqdm(zip(input_1,input_2)):
        # one_1 and one_2 are ids of a pair of reads
        count += 1
        sequence_1,sequence_2 = next(input_1), next(input_2)
        three_1, three_2 = next(input_1), next(input_2)
        qual_1, qual_2 = next(input_1), next(input_2)
        if args.method == "dust":
            score_1, score_2 = dust(sequence_1, k=k),dust(sequence_2, k=k)
        else:
            score_1, score_2 = entropy(sequence_1, k=k),entropy(sequence_2, k=k)
        score = (score_1 + score_2)/2
        scores[int(np.round(score))] += 1
        if score > 10:
            print("low complexcity",sequence_1.strip().decode())
        if score <= 1:
            print("high complexcity",sequence_1.strip().decode())
        if args.filter is not None:
            if (args.method == "entropy" and score < args.filter) or (args.method == "dust" and score > args.filter):
                filtered += 1
                continue
        # add complexity score to read ids
        score_1, score_2 = str(int(np.round(score_1))).encode(), str(int(np.round(score_2))).encode()
        one_1 = one_1.strip() + b":" + score_1 + b"\n"
        one_2 = one_2.strip() + b":" + score_2 + b"\n"
        if args.out_fastq_1 and args.out_fastq_2:
            # output to fastq 1
            output_1.write(one_1)
            output_1.write(sequence_1)
            output_1.write(three_1) 
            output_1.write(qual_1) 
            # output to fastq 2
            output_2.write(one_2)
            output_2.write(sequence_2)
            output_2.write(three_2)
            output_2.write(qual_2)    

    if args.out_fastq_1 and args.out_fastq_2:
        output_1.close()
        output_2.close()
    print(f"{count} reads were processed.")
    if args.filter is not None:
        print(f"{filtered} reads were filtered.")
    if args.statistics is not None:
        fstat = open(args.statistics,"w")
        for i in range(101):
            print(f"{i}\t{int(scores[i])}",file=fstat)
        fstat.close()

if __name__ == "__main__":
    main()
