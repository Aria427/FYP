#!/usr/bin/env python
#This file contains the map step for the aligment MapReduce implementation.

import sys
import gzip
from itertools import chain, islice
from operator import itemgetter

ALPHABET = ['A', 'C', 'G', 'T']
O = {} #dictionary with keys = $, A, C, G, T and values = arrays of counts
C = {} #number of lexographically greater symbols in the reference genome
D = [] # estimated lower bounds of mismatches in the short read

#rewards/penalties => score
GAP = 1
MISMATCH = 1
MATCH = 0

prunes = 0

#This function reads the genome in chunks.
def readGenomeChunks(s3File, bytesNum=100):
    with gzip.open(s3File, 'r') as f:
        f.readline() #ignore header line with genome information
        for chunk in iter(lambda: f.read(bytesNum), ''):
            data = chunk.rstrip().upper().replace('N', '').replace('\n', '').replace(' ', '')
            yield data

#This function reads the whole genome from the file.
def readGenome(s3File):
    genome = ''
    with gzip.open(s3File, 'r') as f: 
        for line in f:
            if line and line[0] != '>': #ignore header line with genome information
                l = line.rstrip().upper().replace('N', '').replace(' ', '')
                genome += l
        return genome

#This function reads the sequencing reads after optimisation as input to the mapper.      
def readOptimisedReads(file):
    sequence, quality = '', ''
    while True: #runs until EOF
    	line = file.readline() 
        if not line: #reached EOF
             break
            
        line = line.split()
        sequence = line[0]
        quality = line[1]

        yield sequence, quality

def readInputPhiXReads(file):  
    readID, sequence, quality = '', '', ''    
    while True: #runs until EOF
        line = file.readline() 
        if not line: #reached EOF
            break

        if line.startswith('@'): #first line of read/record
            #reset to default values
            readID = line.rstrip()
            sequence = ''
            quality = ''   

        elif not readID: #if no previous line starts with @
            readID = line.rstrip() #get first ID
            continue
        
        elif not sequence or not quality:
            sequenceLines = [] 
            while not line.startswith('+'): #not placeholder line (third line)
                #rstrip() - removes leading/trailing whitespace
                #replace() - removes whitespace from within string
                N = [pos for pos, char in enumerate(sequence) if char == 'N'] #positions of N in read
                line = line.rstrip().upper().replace(' ', '')
                sequenceLines.append(line) #no whitespace in string sequence
                line = file.readline()
            sequence = ''.join(sequenceLines) #merge lines to form original sequence
            sequenceNoNs = sequence.replace('N', '') #remove Ns
            temp = sequenceNoNs
        
            qualityLines = []
            qualityNoNs = ''
            
            while True: #collect base qualities
                line = line.rstrip().replace(' ', '')
                qualityLines.append(line) 
                quality = ''.join(qualityLines) #merge lines to form quality
                if len(quality) >= len(sequence): #bases and qualities line up
                    break
                else:
                    line = file.readline()
            
            for i in range(len(quality)): 
                if i not in N: #remove indices corresponding to Ns in read
                    qualityNoNs = qualityNoNs + quality[i]
                
            yield temp, qualityNoNs

#This function generates the suffix array for some text.
#It implements the algorithm of Vladu and Negruseri:
#   http://web.stanford.edu/class/cs97si/suffix-array.pdf
def suffixArray(text):
    L = sorted((b, i) for i, b in enumerate(text)) #list of pairs = (base, index)
    n = len(text)+1
    count = 1
    
    while count < n:
        U = [0] * n
        for (r, i), (s, j) in zip(L, islice(L, 1, None)):
            U[j] = U[i] + (r != s)
        #Invariant: U[i] 
        #   index of text[i:i+count] in sorted list of the text's unique substrings of length count

        L = sorted(chain((((U[i],  U[i+count]), i) for i in range(n - count)),
                         (((U[i], -1), i) for i in range(n - count, n))))
        #Invariant: L[i][1]
        #   starting index in text of substring i of length count in sorted order
        
        count *= 2
        
    return [i for _, i in L]

#This function generates the Burrows-Wheeler Transform of some text using the suffix array.
def bwt(text):
    bwt = []

    for i in suffixArray(text):
        if i == 0:
            bwt.append('$')
        else:
            bwt.append(text[i-1])
    
    return ''.join(bwt)

#This function generates a list of the number of occurences of a character for each reference subsequence.
def rank(bwt):
    ranks, totals = {}, {}

    for char in ALPHABET:
        if (char not in totals) and (char != '$'):
            ranks[char] = []
            totals[char] = 0

    for char in bwt:
        if char != '$':
            totals[char] += 1
        for t in totals.iterkeys():
            ranks[t].append(totals[t])

    return ranks, totals

#This function calculates the number of lexographically greater symbols C in the reference genome.
def computeC(totals):
    C = {'A': 0, 'C': 0, 'G': 0, 'T': 0}

    for k in ALPHABET:
        for ref in ALPHABET:
            if ref < k: C[k] += totals[ref]

    return C

#This function calculates the estimated lower bounds of mismatches in the short read.
def computeD(pattern, C, Oprime, bw):
    k = 1
    l = len(bw)-2
    z = 0
    D = [0] * len(pattern)

    for i in range(0, len(pattern)):
        k = C[pattern[i]] + Oprime[pattern[i]][k-1] + 1
        l = C[pattern[i]] + Oprime[pattern[i]][l]
        if k > l:
            k = 1
            l = len(bw)-1
            z = z+1
        D[i] = z

    return D

#This function ensures that if D[i] = -1, its value will be considered as 0.
def getD(i):
    if i < 0:
        return 0
    else:
        return D[i]

#This function operates in a similar manner to getD().
def getO(char, index):
    if index < 0:
        return 0
    else:
        return O[char][index]

#This function implements the Burrows Wheeler approximate matching algorithm using recursion.
def burrowsWheelerApproximateRecursion(pattern, i, maxMismatches, k, l):    
    global prunes
    saIndex = set() #set of suffix array indices at which a match starts 
    
    #pruning based on estimated mistakes
    if maxMismatches < getD(i):
        prunes += 1
        return set()

    #end of query condition
    temp = set()
    if i < 0:
        for j in range(k,l+1):
            temp.add((j, maxMismatches))
        return temp
    
    #Insertion
    saIndex = saIndex.union(burrowsWheelerApproximateRecursion(
                            pattern, i-1, maxMismatches-GAP, k, l))

    for char in ALPHABET:
        tempk = C[char] + getO(char, k-1) + 1
        templ = C[char] + getO(char, l)

        if tempk <= templ:
            #Deletion
            saIndex = saIndex.union(burrowsWheelerApproximateRecursion(
                                    pattern, i, maxMismatches-GAP, tempk, templ))
            
            if char == pattern[i]: #match
                saIndex = saIndex.union(burrowsWheelerApproximateRecursion(
                                        pattern, i-1, maxMismatches+MATCH, tempk, templ))
                
            else: #mismatch
               saIndex = saIndex.union(burrowsWheelerApproximateRecursion(
                                        pattern, i-1, maxMismatches-MISMATCH, tempk, templ))

    return saIndex

#This function implements the Burrows Wheeler approximate matching algorithm.
#bw = bwt of genome, bwr = bwt of reverse genome, pattern = short read to match, maxMismatches = max number of differences 
def burrowsWheelerApproximate(bw, bwr, pattern, maxMismatches):    
    global prunes
    global O
    global C
    global D
    
    O, totals = rank(bw)
    Oprime, _ = rank(bwr) #reverse ranks
    C = computeC(totals) 
    D = computeD(pattern, C, Oprime, bw)

    saIndexSet = burrowsWheelerApproximateRecursion(pattern, len(pattern)-1, maxMismatches, 0, len(bw)-1)
    indexDict = {}

    for (i, j) in saIndexSet:
        if i in indexDict: #index already exists
            if indexDict[i] < j:
                indexDict[i] = j #pick higher maxMismatches value
                prunes += 1
        else:
            indexDict[i] = j

    return sorted(indexDict.items(), key=itemgetter(1), reverse=True) #sort by maxMismtaches in descending order                  

#This function aligns the reads to the genome using the Burrows Wheeler approximate algorithm.
def alignBurrowsWheeler(read, quality, bw, bwr):
    readQualityDictionary = {} #key:read, value:list of quality integers
    
    #maximum number of mismatches = 5
    matchOffset = burrowsWheelerApproximate(bw, bwr, read, 5) #check if read matches in forward/backward direction of genome
        
    if len(list(matchOffset)) > 0: #match - read aligned in at least one place
        readQualityDictionary[read] = quality
        
    return matchOffset, readQualityDictionary
    
def main():
    genomeFile = '/home/aria427/test/data/HumanGenome_Part100Update.gz' #-cacheArchive s3://fyp-input/HumanGenome.fa.gz#human  
    preprocGenome = '/home/aria427/test/preprocessed/BurrowsWheelerGenome'
    
    genome = readGenome(genomeFile)
    readSeq = readOptimisedReads(sys.stdin) #Human reads=28,094,847
    
    f = open(preprocGenome, 'r')
    bw, bwr = f.read().rstrip().split('\n', 1)
    
    for read, quality in readSeq:
        #Human genome=64,185,939 lines -> 3,273,481,150 bytes
        offset, rqDict = alignBurrowsWheeler(read, quality, bw, bwr)
         
        #write results to STDOUT (standard output)
        for o in offset: #to remove empty list and '[' ']' characters
            #tab-delimited, key:offset of match with reads, value:<default count of 1, genome subsequence matched, read matched, corresponding quality> 
            print '%s\t%s\t%s\t%s\t%s' % (o[0], 1, genome[o[0]:o[0]+len(read)], read, quality) #[0] as output is tuple (offset, maxMismatches)
	    #The output here will be the input for the reduce step
        
    f.close()
        
if __name__ == '__main__':
    main() 
    