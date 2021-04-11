#!/usr/bin/env python
import argparse
import sys
import gzip
import re
from tqdm import tqdm
import logging


logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')


def main():
    parser = argparse.ArgumentParser(description='Fetch reads from fastq file with reads id')
    parser.add_argument('--input', '-i', type=str, default="-",help='Input fastq file')
    parser.add_argument('--output','-o',type=str, default="-",help='Output fastq file')
    parser.add_argument('--read-ids','-r',type=str,required=True,help="Read ids to fetch")
    parser.add_argument('--fasta','-f',action="store_true",default=False,help="Whether output fasta. Output fastq by default")
    args = parser.parse_args()
    gzin = False
    gzout = False
    if args.input == "-":
        fin = sys.stdin
    elif args.input.endswith(".gz"):
        gzin = True
        fin = gzip.open(args.input,"rb")
    else:
        fin = open(args.input,"r")
    
    if args.output == "-":
        fout = sys.stdout
    elif args.output.endswith(".gz"):
        gzout = True
        fout = gzip.open(args.output,"wb")
    else:
        fout = open(args.output,"w")
    logging.info("Load read ids ...")
    read_ids = set(open(args.read_ids).read().strip().split("\n"))
    logging.info("Done .")

    logging.info("Fetch reads with their ids ...")
    for line in tqdm(fin):
        if gzin:
            line = line.decode()
        assert line.startswith("@"),"Malformated fastq input"
        read_id = re.split(r"\s",line[1:].strip())[0]
        if read_id in read_ids:
            if args.fasta:
                line = ">" + read_id + "\n"
            if gzout:
                line = line.encode()
            fout.write(line)             
            for i in range(3):
                line = next(fin)
                if i > 0 and args.fasta:
                    continue
                if gzin == gzout:
                    fout.write(line)
                elif gzin:
                    fout.write(line.decode())
                else:
                    fout.write(line.encode())
        else:
            for i in range(3):
                line = next(fin)
    logging.info("Done .")
    fin.close()
    fout.close()   
    
                
    
    
    


if __name__ == "__main__":
    main()
