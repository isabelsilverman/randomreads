#!/usr/bin/env python
# coding: utf-8

import sys, getopt
from collections import Counter
import csv

def main(argv):
    #get parameters from command line
    inFile = "reads.fq"
    outFile = "freq.tsv"

    opts, args = getopt.getopt(argv,"i:o:",["inFile=","outFile="])  
    for opt, arg in opts:
        if opt in ("-i", "--inFile"):
            inFile = arg
        elif opt in ("-o", "--outFile"):
            outFile = arg

    #reads sequences from inFile
    with open(inFile) as file:
        reads = []
        for count, line in enumerate(file):
            if count % 4 == 1:
                reads.append(line.strip())    

    #parse reads into codons 
    codons = [[] for i in range(28)]
    for read in reads:
        for i in range(0,28):
            codons[i].append(read[i:i+3])

    #calculate codn frequency at each position
    freq = []
    for row in codons:
        freq.append(Counter(row))

    #write frequencies into tsv
    with open(outFile, 'w') as tsv:
        writer = csv.writer(tsv, delimiter='\t')
        writer.writerow(["Position", "Codon", "Frequency"])
        for i in range(0,28):
            for item in freq[i]:
                writer.writerow([str(i+1), item, freq[i][item]])

if __name__ == "__main__":
   main(sys.argv[1:])