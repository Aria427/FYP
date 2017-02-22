#!/usr/bin/env python
#This file contains some testing functions with regards to finding repetitive patterns in the original genome.

import gzip
from collections import Counter

#This function uses the sliding window approach with the aid of a generator
#in order to process the genome sequence bit-by-bit (handling overlapping regions). 
def slidingWindow(sequence, windowSize, stepSize):
    chunksCount = ( (len(sequence) - windowSize) / stepSize) + 1
    for i in xrange(0, chunksCount*stepSize, stepSize):
        yield sequence[i:i+windowSize]

#This function returns a list of k-mers, similar to the sliding window approach.
def kmerList(sequence, k):
    kmers = []
    for i in xrange(0, len(sequence) + 1 - k):
        kmers.append( sequence[i:i+k] )
    return kmers

#This function returns the frequency of each 4-letter (int) word found in the genome.
def countIntWords(genome):
    wordCount = Counter()
    with gzip.open(genome) as f:
        #read genome into line-by-line generator
        subseqs = (line.upper().replace('\n', '').replace('N', '')
                    for line in f if line and line[0] != '>') #ignore header line
        
        last3 = '' #store the last 3 (4-1) bases of each line
        for s in subseqs:
            s = last3 + s #append to start of next line to handle patterns found between lines
            kmers = kmerList(s, 4)  #generate 4-mers of line
            wordCount.update(kmers) #update counter
            last3 = s[-3:]
    return wordCount
   
#This function returns the frequency of each 8-letter (long) word found in the genome.
def countLongWords(genome):
    wordCount = Counter()
    with gzip.open(genome) as f:
        #read genome into line-by-line generator
        subseqs = (line.upper().replace('\n', '').replace('N', '')
                    for line in f if line and line[0] != '>') #ignore header line
        
        last7 = '' #store the last 7 (8-1) bases of each line
        for s in subseqs:
            s = last7 + s #append to start of next line to handle patterns found between lines
            kmers = kmerList(s, 8)  #generate 8-mers of line
            wordCount.update(kmers) #update counter
            last7 = s[-7:]
    return wordCount    
           