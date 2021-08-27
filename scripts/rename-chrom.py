#!/usr/bin/env python
import argparse
from tqdm import tqdm

def main():
    parser = argparse.ArgumentParser(description='Rename sequence name in bed file')
    parser.add_argument('--input', '-i', type=str, required=True,help='Input bed file')
    parser.add_argument('--output','-o', type=str, required=True,help='Output bed file')
    parser.add_argument('--name-mapping','-n',type=str, required=True,help="Name mapping")
    args = parser.parse_args()
    name_mapping = {}
    with open(args.name_mapping) as f:
        for line in f:
            key, value = line.strip().split("\t") 
            name_mapping[key] = value
    fout = open(args.output,"w")
    with open(args.input) as f:
        for line in tqdm(f):
            fields = line.strip().split("\t")
            if fields[0] not in name_mapping:
                print(f"{fields[0]} not in name dict, skip")
                continue
            fields[0] = name_mapping[fields[0]]
            line = "\t".join(fields) + "\n"
            fout.write(line)
    fout.close()

if __name__ == "__main__":
    main()
