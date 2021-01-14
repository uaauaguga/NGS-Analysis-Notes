#!/usr/bin/env python
import gzip
import sys
import re
import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-i','--input_file',help='Path and prefix of  input *_[12].fastq.gz file',required=True)
args = parser.parse_args()


input_file_1 = args.input_file+'_1.fastq.gz'
input_file_2 = args.input_file+'_2.fastq.gz'

input_1 =  gzip.open(input_file_1)
input_2 =  gzip.open(input_file_2)

unpaired = False
for i,(label_1,label_2) in enumerate(zip(input_1,input_2)):
    if i%4 != 0:
        continue
    if label_1 != label_2:
        unpaired = True
        break

if unpaired:
    print("Unpaired reads detected")
else:
    print("No unpaired reads detected")
    


        
       
 
       

