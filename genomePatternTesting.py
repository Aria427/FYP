#!/usr/bin/env python
#This file contains some testing functions with regards to finding repetitive patterns in the original genome.

import gzip
import mmap
from collections import Counter

#This function reads a file in pieces/chunks using a lazy approach (generator).
def readInChunks(genomeFile, chunkSize=1024):
    with gzip.open(genomeFile, 'r') as f:
        while True:
            data = f.read(chunkSize).upper().replace('N', '').replace(' ', '')
            if not data:
                break
            yield data
    #return iter(lambda: genomeFile.read(chunkSize), '')
        
#This function returns the frequency of each 4-letter (int) word found in the genome.
#The implementation is based on the readInChunks() function defined above.
#This function lasts 2536.83899999 seconds when run on the human genome;
#i.e. the most efficient out of the constructed functions
def countIntWordsChunks(genome):
    wordCount = Counter()
    for piece in readInChunks(genome):
        nucleotides = filter(str.isalpha, piece)
        int1 = nucleotides[:-3] 
        int2 = nucleotides[:-2]
        int3 = nucleotides[:-1]
        int4 = nucleotides[1:]
        wordCount += Counter([i1+i2+i3+i4 for i1, i2, i3, i4 in zip(int1, int2, int3, int4)])
    print wordCount

#This function returns the frequency of each 4-letter (int) word found in the genome.
#The implementation is based on the readBytes() function defined within this one. 
#This function lasts 2614.78399992 seconds when run on the human genome;
#i.e. takes longer than the above function. 
def countIntWordsBytes(genome):
    wordCount = Counter()
    f = gzip.open(genome, 'r')
    def readBytes(): #This function reads a file 1024 bytes at a time.
        return f.read(1024)  
    for piece in iter(readBytes, ''):
        nucleotides = filter(str.isalpha, piece)
        int1 = nucleotides[:-3] 
        int2 = nucleotides[:-2]
        int3 = nucleotides[:-1]
        int4 = nucleotides[1:]
        wordCount += Counter([i1+i2+i3+i4 for i1, i2, i3, i4 in zip(int1, int2, int3, int4)])
    print wordCount 
    f.close()
  
#This function returns the frequency of each 4-letter (int) word found in the genome.
#This function crashes when run on the human genome.
def countIntWords(genome):
    with open(genome, 'r') as f:
        m = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ) #memory-map whole file (0->whole file)
        nucleotides = filter(str.isalpha, m.read(-1)) #filter->ignore new line characters; -1->return all bytes
        int1 = nucleotides[:-3] 
        int2 = nucleotides[:-2]
        int3 = nucleotides[:-1]
        int4 = nucleotides[1:]
        wordCount = Counter([i1+i2+i3+i4 for i1, i2, i3, i4 in zip(int1, int2, int3, int4)])
        print wordCount

"""
from collections import Counter
from itertools import chain

with open(genomeFile, 'r') as f:
    prev = f.read(1)
    c = Counter()
    for ch in filter(str.isalpha, chain.from_iterable(f)): #filter to remove new line characters
        next1 = f.read(1)
        next2 = f.read(1)
        c[prev + ch + next1 + next2] += 1
        prev = ch
print c 
"""
        