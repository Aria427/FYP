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

"""
#The below functions proved to be inefficient.

#This function reads a file in pieces/chunks using a lazy approach (generator).
def readInChunks(genomeFile, chunkSize=1024):
    with gzip.open(genomeFile, 'r') as f:
        while True:
            data = f.read(chunkSize).rstrip().upper().replace('N', '').replace(' ', '')
            if not data:
                break
            yield data
    #return iter(lambda: genomeFile.read(chunkSize), '')
    
#This function returns the frequency of each 4-letter (int) word found in the genome.
#The implementation is based on the slidingWindow() function defined above.
def countIntWords(genome):
    wordCount = Counter()
    for piece in readInChunks(genome, 4096):
        #window size = 4 as int (4-letter matches) is considered
        #step size = 1 to consider each nucleotide
        chunks = slidingWindow(piece, 4, 1)
        for c in chunks:
            wordCount += Counter(c)
    return wordCount

#This function returns the frequency of each 8-letter (long) word found in the genome.
#The implementation is based on the slidingWindow() function defined above.
def countLongWords(genome):
    wordCount = Counter()
    with gzip.open(genome) as f:
        for line in f:
            line = line.rstrip().upper().replace('N', '').replace('\n', '').replace(' ', '')
            #window size = 8 as long (8-letter matches) is considered
            #step size = 1 to consider each nucleotide
            chunks = slidingWindow(line, 8, 1) 
            for c in chunks:
                wordCount += Counter(c)
    return wordCount
"""              