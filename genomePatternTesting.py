#!/usr/bin/env python
#This file contains some testing functions with regards to finding repetitive patterns in the original genome.

import gzip
import mmap
from collections import Counter
import re

#This function returns the frequency of each 4-letter (int) word found in the genome.
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
        