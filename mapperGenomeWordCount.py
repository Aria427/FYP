#!/usr/bin/env python

import genomePatternTesting
import sys
from path import path

#genomeFile = path('Data\HumanGenome.fa.gz').abspath()
words = []

#input file comes from STDIN (standard input)
for piece in genomePatternTesting.readInChunks(sys.stdin):#genomeFile):
    nucleotides = filter(str.isalpha, piece)
    int1 = nucleotides[:-3] 
    int2 = nucleotides[:-2]
    int3 = nucleotides[:-1]
    int4 = nucleotides[1:]
    words.append(int1+int2+int3+int4)

for word in words:
    #write results to STDOUT (standard output)
    print '%s\t%s' % (word, 1) #tab-delimited; the trivial word count is 1
#The output here will be the input for the reduce step.

