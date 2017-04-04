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

#This function reads the genome in chunks (groups of bytes) from the file stored in S3.
def readGenomeChunks(s3File, bytesNum=100):
    with gzip.open(s3File, 'r') as f:
        f.readline()
        for chunk in iter(lambda: f.read(bytesNum), ''):
            data = chunk.rstrip().upper().replace('N', '').replace('\n', '').replace(' ', '')
            yield data

#This function reads the whole genome from the file stored in S3.
def readGenome(s3File):
    genome = ''
    with gzip.open(s3File, 'r') as f: 
        for line in f:
            if line and line[0] != '>': #ignore header line with genome information
                l = line.rstrip().upper().replace('N', '').replace(' ', '')
                genome += l
        return genome
                  
#This function reads the sequencing reads as input to the mapper.        
def readInputReads(file):
    flag, sequence, quality = '', '', ''
    while True: #runs until EOF
        line = file.readline() 
        if not line: #reached EOF
            break
        
        if line.startswith('#'): #read details
            line = file.readline()
            pass
        
        elif line.startswith('>'): #>flags reads scores
            line = file.readline()
            pass
        
        else:
            line = line.split()
            flag = line[0]
            sequence = line[1]
            quality = line[2]
            
            #Each read has length = 60
            if sequence == 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN':
                line = file.readline()
                pass
            
            else:
                qualityNoNs = ''
                N = [pos for pos, char in enumerate(sequence) if char == 'N'] #positions of N in read
                sequenceNoNs = sequence.rstrip().upper().replace('N', '').replace(' ', '') 
                
                for i in range(len(quality)):
                    if i not in N: #remove indices corresponding to Ns in read
                        qualityNoNs = qualityNoNs + quality[i]
                
                yield sequenceNoNs, qualityNoNs

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
                
#This function finds the reverse complement of a sequencing read.   
def reverseComplement(read):
    complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'} #each base is associated with its complementary base
    reverseRead = ''
    for base in read:
        reverseRead = complement[base] + reverseRead #complement added to beginning in order to reverse the read from end to start
    return reverseRead

#This function takes the quality value Q (rounded integer) and converts it to its respective character. 
def QtoPhred33(Q):
    return chr(Q + 33) #converts integer to character according to ASCII table

#This function takes the Phred-33 encoded character and converts it back to Q.  
def phred33ToQ(qual):
    return ord(qual) - 33 #converts character to integer according to ASCII table

#This function aligns the reads to the genome using the Burrows Wheeler approximate algorithm.
def alignBurrowsWheeler(read, quality, genome):
    readQualityDictionary = {} #key:read, value:list of quality integers
    
    #sa = bwa.suffixArray(genome)
    bw = bwt(genome) 
    bwr = bwt(genome[::-1]) 
    
    #maximum number of mismatches = 2
    reverseRead = reverseComplement(read)
    matchOffset = burrowsWheelerApproximate(bw, bwr, read, 2) #check if read matches in forward direction of genome
    matchOffset.extend(burrowsWheelerApproximate(bw, bwr, reverseRead, 2)) #add results of any matches in reverse complement of genome
        
    if len(list(matchOffset)) > 0: #match - read aligned in at least one place
        qualityQ = []
        for q in quality:
            qualityQ.append(phred33ToQ(q))
        readQualityDictionary[read] = qualityQ
        
    return matchOffset, readQualityDictionary
    
def main():
    #hard-coded reference genome stored in S3 via Amazon EMR
    #genomeFile = 's3://fyp-input-gen/HumanGenome_200000.fa.gz' #Frankfurt region doesn't work, Ireland does
    genomeFile = './human' #-cacheArchive s3://fyp-input/HumanGenome.fa.gz#human
    #g = readGenome(genomeFile)   
    readSeq = readInputReads(sys.stdin) #Human reads=28,094,847
    
    for read, quality in readSeq: 
        #Human genome=64,185,939 lines -> 3,273,481,150 bytes
        genome = readGenomeChunks(genomeFile, 250000) #250,000 bytes = 0.23842MB
        overlap = '' #size of read-1
        filesOffset = 0 #file is split in chunks so offset needs to change according to chunk
        
        for g in genome:
            g = overlap + g
            offset, rqDict = alignBurrowsWheeler(read, quality, g)
         
            #write results to STDOUT (standard output)
            for o in offset: #to remove empty list and '[' ']' characters
                #tab-delimited, key:offset of match with reads, value:<default count of 1, genome subsequence matched, read matched, corresponding quality> 
                print '%s\t%s\t%s\t%s\t%s' % (o[0]+filesOffset, 1, g[o[0]:o[0]+len(read)], read, quality) #[0] as output is tuple (offset, maxMismatches)
                #The output here will be the input for the reduce step  
            
            overlap = g[-99:] #100-1 for PhiX, 60-1 for Human read => -1 as last 60 have already been read
            filesOffset += (len(g)-100) #store offset according to overlap as file is read in chunks

        
if __name__ == '__main__':
    main()            
   