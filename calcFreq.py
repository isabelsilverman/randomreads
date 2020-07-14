#!/usr/bin/env python
# coding: utf-8

import sys, getopt
from collections import Counter

def getCounts(inFile, length):
    #reads sequences from inFile
    with open(inFile) as file:
        reads = []
        for count, line in enumerate(file):
            if count % 4 == 1:
                reads.append(line.strip())    

    #parse reads into codons 
    codons = [[] for i in range(length)]
    for read in reads:
        for i in range(0,length):
            codons[i].append(read[3*i:3*i+3])

    #calculate codon counts at each position
    counts = []
    for row in codons:
        counts.append(Counter(row))
    
    return counts, len(reads)

def main(argv):
    #get parameters from command line
    inFile = []
    outFile = ""
    wtFile = ""

    opts, args = getopt.getopt(argv,"i:o:w:",["inFile=","outFile=","wtSeq="])  
    for opt, arg in opts:
        if opt in ("-i", "--inFile"):
            arg = arg.split(",")
            if type(arg) is str:
                inFile.append(arg)
            else:
                inFile.extend(arg)
        elif opt in ("-o", "--outFile"):
            outFile = arg
        elif opt in ("-w", "--wtSeq"):
            wtFile = arg
    
    #read in wild-type amino acid sequence
    with open(wtFile) as file:
        wt_seq = list(file.read()) 
    length = len(wt_seq)

    #get counts for each input file
    counts = []
    for file in inFile:
        counts.append(getCounts(file, length))


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
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*', 
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'} 

    codons = list(table.keys())

    #write counts and frequencies into tsv
    import csv
    with open(outFile, 'w') as tsv:
        writer = csv.writer(tsv, delimiter='\t')
        colnames = ["Position", "Codon", "Amino Acid"]
        for file in inFile:
            colnames.append(file[:-3] + " Count")
        for file in inFile:
            colnames.append(file[:-3] + " Frequency")
        writer.writerow(colnames)
        for i in range(0,length):
            for j in range(0,64):
                codon = codons[j]
                row = [str(i+1), codon, table[codon]]
                for item, numReads in counts:
                    row.append(item[i][codon])
                for item, numReads in counts:
                    row.append(item[i][codon]/numReads)
                writer.writerow(row)

    

    #read tsv file into a data frame
    import pandas as pd
    df = pd.read_csv("freq.tsv", delimiter = '\t')

    #calculate A 
    import math
    import numpy as np
    for file in inFile[1:]:
        name = file[:-3]
        wt_freq = []
        for i in range(0, length):
            wt_freq.append(sum(df[(df['Amino Acid'] == wt_seq[i]) & (df['Position'] == i+1)][name + ' Frequency']))
        wt_freq = np.repeat(wt_freq, 64)
        df[name + ' A'] = df[name + ' Frequency']/wt_freq
        df[name + ' A'] = df[name + ' A'].apply(math.log, args = [2])

    #calculate F and s
    for file in inFile[2:]:
        name = file[:-3]
        df[name + ' F'] = df[name + ' A'] - df[inFile[1][:-3] + ' A']
        wt_syn = []
        for i in range(0, length):
            wt_syn.append(sum(df[(df['Amino Acid'] == wt_seq[i]) & (df['Position'] == i+1)][name + ' F']))
        wt_syn = np.repeat(wt_syn, 64)
        stop = []
        for i in range(0, length):
            stop.append(sum(df[(df['Amino Acid'] == '*') & (df['Position'] == i+1)][name + ' F']))
        stop = np.repeat(stop, 64)
        df[name + ' s'] = (df[name + ' F'] - wt_syn)/(wt_syn - stop)
        reps = []
        new_df = df[['Position','Amino Acid',name + ' s']]
        new_df.columns = ['Position','Amino Acid','fitness effect']
        reps.append(new_df)

    #create a new data frame with the fitness effect for every amino acid at each position
    new_df_total = pd.concat(reps)
    mat = pd.pivot_table(new_df_total, index = 'Amino Acid', columns = 'Position', values = 'fitness effect', aggfunc = np.median)

    #plot the heatmap
    from matplotlib import pyplot as plt
    from matplotlib.patches import Rectangle
    import seaborn as sns

    sns.set()
    ax = sns.heatmap(mat, yticklabels = 1, square = 'TRUE')
    for i in range(0,length):
        ax.add_patch(Rectangle((i, list(mat.index).index(wt_seq[i])), 1, 1, fill=False, edgecolor='blue', lw=3))
    ax.set_title("Fitness Effect")
    plt.savefig("heatmap.png", bbox_inches='tight')



if __name__ == "__main__":
   main(sys.argv[1:])