#!/usr/bin/env python
# coding: utf-8

import sys, getopt
import random

def main(argv):
    #get parameters from command line
    numReads = 100000
    outFile = "reads.fq"

    opts, args = getopt.getopt(argv,"n:o:",["numReads=","outFile="])    
    for opt, arg in opts:
        if opt in ("-n", "--numReads"):
            numReads = arg
        elif opt in ("-o", "--outFile"):
            outFile = arg

    #generate random sequences
    dna = ["A","G","C","T"]

    with open(outFile, "w") as reads:
        for i in range(0, int(numReads)):
            reads.write("@SRR000000.1 read:" + str(i) + " length=30\n")
            for i in range(0,30):
                reads.write(random.choice(dna))
            reads.write("\n+\n||||||||||||||||||||||||||||||\n")

if __name__ == "__main__":
   main(sys.argv[1:])
