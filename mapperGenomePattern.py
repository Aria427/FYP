#!/usr/bin/env python

import sys
from collections import Counter
#from path import path

#genomeFile = path('Data\HumanGenome.fa.gz').abspath()
def readInput(file):
    for line in file:
        if line and line[0] != '>':
            yield line.upper().replace('\n', '').replace('N', '')

#This function returns a list of k-mers, similar to the sliding window approach.
def kmerList(sequence, k):
    kmers = []
    for i in xrange(0, len(sequence) + 1 - k):
        kmers.append( sequence[i:i+k] )
    return kmers            
            
def main():
    #input file comes from STDIN (standard input)
    subseqs = readInput(sys.stdin)
    
    wordCount = Counter()
    last3 = '' #store the last 3 (4-1) bases of each line
    for s in subseqs:
        s = last3 + s #append to start of next line to handle patterns found between lines
        kmers = kmerList(s, 4)  #generate 4-mers of line
        wordCount.update(kmers)
        last3 = s[-3:]
    
    #write results to STDOUT (standard output)
    for key, value in wordCount.items():
        print '%s\t%s' % (key, value) #tab-delimited
#The output here will be the input for the reduce step.

if __name__ == '__main__':
    main()