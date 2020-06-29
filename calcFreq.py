#!/usr/bin/env python
# coding: utf-8

import sys, getopt
from collections import Counter
import pandas as pd
import csv

def getCounts(inFile):
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

    #calculate codon counts at each position
    counts = []
    for row in codons:
        counts.append(Counter(row))
    
    return counts, len(reads)

def main(argv):
    #get parameters from command line
    inFile = []
    outFile = "freq.tsv"

    opts, args = getopt.getopt(argv,"i:o:",["inFile=","outFile="])  
    for opt, arg in opts:
        if opt in ("-i", "--inFile"):
            arg = arg.split(",")
            if type(arg) is str:
                inFile.append(arg)
            else:
                inFile.extend(arg)
        elif opt in ("-o", "--outFile"):
            outFile = arg
    
    #get counts for each input file
    counts = []
    for file in inFile:
        counts.append(getCounts(file))


    table = { 
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'} 

    codons = ['AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', 
            'AGA', 'AGC', 'AGG', 'AGT', 'ATA', 'ATC', 'ATG', 'ATT', 
            'CAA', 'CAC', 'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT', 
            'CGA', 'CGC', 'CGG', 'CGT', 'CTA', 'CTC', 'CTG', 'CTT', 
            'GAA', 'GAC', 'GAG', 'GAT', 'GCA', 'GCC', 'GCG', 'GCT', 
            'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT', 
            'TAA', 'TAC', 'TAG', 'TAT', 'TCA', 'TCC', 'TCG', 'TCT', 
            'TGA', 'TGC', 'TGG', 'TGT', 'TTA', 'TTC', 'TTG', 'TTT']

    #write counts and frequencies into tsv
    with open(outFile, 'w') as tsv:
        writer = csv.writer(tsv, delimiter='\t')
        colnames = ["Position", "Codon", "Amino Acid"]
        for file in inFile:
            colnames.append(file[:-3] + " Count")
        for file in inFile:
            colnames.append(file[:-3] + " Frequency")
        writer.writerow(colnames)
        for i in range(0,28):
            for j in range(0,64):
                codon = codons[j]
                row = [str(i+1), codon, table[codon]]
                for item, numReads in counts:
                    row.append(item[i][codon])
                for item, numReads in counts:
                    row.append(item[i][codon]/numReads)
                writer.writerow(row)

if __name__ == "__main__":
   main(sys.argv[1:])