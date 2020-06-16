#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys, getopt

numReads = 100000
outFile = "reads.fq"

opts, args = getopt.getopt(sys.argv,"no",["numReads","outFile"])
for opt, arg in opts:
    if opt in ("-n", "--numReads"):
        numReads = arg
    elif opt in ("-o", "--outFile"):
        outFile = arg
   


# In[2]:


import random
dna = ["A","G","C","T"]

with open(outFile, "w") as reads:
    for i in range(0, numReads):
        reads.write("@SRR000000.1 read:" + str(i) + " length=30\n")
        for i in range(0,30):
            reads.write(random.choice(dna))
        reads.write("\n+\n||||||||||||||||||||||||||||||\n")


# In[ ]:




