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
#This function lasts ~2310.65400004 seconds when run on the human genome;
#i.e. the most efficient out of the constructed functions for int
def countIntWordsChunks(genome):
    wordCount = Counter()
    for piece in readInChunks(genome):
        nucleotides = filter(str.isalpha, piece) #isalpha -> only alphabetic characters
        #nucleotides = [base in nucleotides for base in 'ACGT']
        int1 = nucleotides[:-3] 
        int2 = nucleotides[:-2]
        int3 = nucleotides[:-1]
        int4 = nucleotides[1:]
        wordCount += Counter([i1+i2+i3+i4 
                              for i1, i2, i3, i4 
                              in zip(int1, int2, int3, int4)])
    print wordCount
#Result:
#Counter({'TTTT': 298456206, 'AAAA': 295878805, 'AAAT': 233204521, 'TTTG': 221603116, 
#'CCCA': 220519406, 'AAAG': 213147885, 'CCCT': 212879264, 'TTTA': 196975835, 
#'GGGA': 182228315, 'TTTC': 181239017, 'GGGG': 158825594, 'CCCC': 157575620, 
#'GGGT': 153734106, 'AAAC': 153370668, 'GGGC': 129676917, 'CCCG': 30886900, 
#'CCCH': 433, 'HHHR': 430, 'IIIV': 325, 'KKKI': 325, 'VVVA': 291, 'AAAL': 259, 
#'LLLT': 259, 'RRRK': 212, 'RRRU': 125, 'UUUK': 112, 'GGGL': 88, 'LLLV': 88, 
#'RRRG': 79, 'VVVG': 45, 'AAAD': 42, 'VVVR': 42, 'OOOM': 42, 'DDDO': 42, 
#'RRRA': 42, 'VVVT': 26, 'MMMG': 25, 'VVVC': 21, 'UUUG': 9, 'MMMC': 7, 'JJJH': 6, 
#'HHHV': 6, 'MMMA': 6, 'RRRJ': 6, 'MMMT': 5, 'XXXK': 3, 'RRRX': 3, 'BBBV': 2, 
#'KKKB': 2, 'YYYK': 1, 'RRRY': 1, 'RRRM': 1}    

    
#This function returns the frequency of each 8-letter (long) word found in the genome.
#The implementation is based on the readInChunks() function defined above.
#This function lasts ~2922.11099982 seconds when run on the human genome.
def countLongWordsChunks(genome):
    wordCount = Counter()
    for piece in readInChunks(genome):
        nucleotides = filter(str.isalpha, piece) #isalpha -> only alphabetic characters
        int1 = nucleotides[:-7] 
        int2 = nucleotides[:-6]
        int3 = nucleotides[:-5]
        int4 = nucleotides[:-4]
        int5 = nucleotides[:-3] 
        int6 = nucleotides[:-2]
        int7 = nucleotides[:-1]
        int8 = nucleotides[1:]   
        wordCount += Counter([i1+i2+i3+i4+i5+i6+i7+i8 
                              for i1, i2, i3, i4, i5, i6, i7, i8 
                              in zip(int1, int2, int3, int4, int5, int6, int7, int8)])
    print wordCount   
#Result:
#Counter({'TTTTTTTT': 297263469, 'AAAAAAAA': 294698011, 'AAAAAAAT': 232272322, 
#'TTTTTTTG': 220718195, 'CCCCCCCA': 219638990, 'AAAAAAAG': 212295859, 
#'CCCCCCCT': 212028425, 'TTTTTTTA': 196188540, 'GGGGGGGA': 181499174, 
#'TTTTTTTC': 180514796, 'GGGGGGGG': 158190458, 'CCCCCCCC': 156944851, 
#'GGGGGGGT': 153119333, 'AAAAAAAC': 152756329, 'GGGGGGGC': 129158522, 
#'CCCCCCCG': 30763035, 'CCCCCCCH': 421, 'HHHHHHHR': 419, 'IIIIIIIV': 321, 
#'KKKKKKKI': 321, 'VVVVVVVA': 290, 'LLLLLLLT': 259, 'AAAAAAAL': 259, 
#'RRRRRRRK': 209, 'RRRRRRRU': 118, 'UUUUUUUK': 109, 'LLLLLLLV': 88, 
#'GGGGGGGL': 88, 'RRRRRRRG': 79, 'VVVVVVVG': 44, 'OOOOOOOM': 42, 'RRRRRRRA': 42, 
#'VVVVVVVR': 42, 'DDDDDDDO': 42, 'AAAAAAAD': 42, 'VVVVVVVT': 26, 'MMMMMMMG': 25, 
#'VVVVVVVC': 21, 'UUUUUUUG': 9, 'MMMMMMMC': 7, 'RRRRRRRJ': 6, 'MMMMMMMA': 6, 
#'HHHHHHHV': 6, 'JJJJJJJH': 6, 'MMMMMMMT': 5, 'RRRRRRRX': 3, 'XXXXXXXK': 3, 
#'KKKKKKKB': 2, 'BBBBBBBV': 2, 'RRRRRRRM': 1, 'RRRRRRRY': 1, 'YYYYYYYK': 1})
    
    
#This function returns the frequency of each 4-letter (int) word found in the genome.
#The implementation is based on the readBytes() function defined within this one. 
#This function takes longer than the 'chunks' function to run on the human genome. 
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
        