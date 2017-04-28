#!/usr/bin/env python
#This file contains the map step for the aligment MapReduce implementation.

import sys
import gzip
import bisect

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

#This class is an index object which is applied for k-mer indexing. 
class KmerIndex(object):
    #This function initialises the k-mer object to be used by pre-processing text T.
    def __init__(self, text, k): #initialise index from all k-mer length substrings
        self.k = k  #k-mer length
        self.index = []

        for i in xrange(len(text) - k + 1): #for each index such that k-mer doesn't run past end of T
            self.index.append((text[i:i+k], i)) #add (k-mer, offset) tuple
        self.index.sort()  #sort in ascending order according to k-mer
    
    #This function finds the number of index hits for first k-mer of pattern P
    def query(self, pattern): 
        kmer = pattern[:self.k]  #query with first k-mer
        i = bisect.bisect_left(self.index, (kmer, -1)) #binary search for 1st position in list where k-mer occurs
        indexHits = [] #all indices in T where the first k bases of P match
        
        while i < len(self.index): #collect matching index entries
            if self.index[i][0] != kmer:
                break
            indexHits.append(self.index[i][1]) #append 2nd value of tuple
            i += 1
            
        return indexHits
 
#This function implements approximate matching on k-mer indexing using the pigeonhole principle. 
def kmerIndexApproximate(pattern, text, n, indexObject): #n = max number of mismatches
    segmentLength = int(round(len(pattern) / (n+1))) #n+1 segments where at least 1 must match perfectly against T
    allMatches = set() #all the indices where matches where found
    indexHits = 0
    
    for i in range(n+1): #for each segment in P
        #bounds of P for segment being searched
        start = i*segmentLength
        end = min((i+1)*segmentLength, len(pattern)) #min() to not run past end of P
        
        matches = indexObject.query(pattern[start:end])
        indexHits += 1
        
        #step through each match position to make sure rest of P matches T with no more than n mismatches
        for m in matches:
            if m < start or m-start+len(pattern) > len(text):
                continue #P runs off the start or end of T
            
            mismatches = 0
            for j in range(0, start): #compare segment of P before start against corresponding positions in T
                if not pattern[j] == text[m-start+j]:
                    mismatches += 1
                    if mismatches > n: #exceeds maximum
                        break
                    
            for j in range(end, len(pattern)): #compare suffix after segment
                if not pattern[j] == text[m-start+j]:
                    mismatches += 1
                    if mismatches > n:
                        break
                    
            if mismatches <= n:
                allMatches.add(m - start) #add start position of match
                
    return list(allMatches)                                

#This function aligns the reads to the genome using the k-mer approximate indexing method.
def alignKmer(read, quality, genome):
    readQualityDictionary = {} #key:read, value:list of quality integers

    index = KmerIndex(genome, 10) #k-mer of length 10
    
    #maximum number of mismatches = 5
    matchOffset = kmerIndexApproximate(read, genome, 5, index) #check if read matches in forward/backward direction of genome
        
    if len(list(matchOffset)) > 0: #match - read aligned in at least one place
        readQualityDictionary[read] = quality
     
    return matchOffset, readQualityDictionary
    
def main():
    genomeFile = '/home/aria427/test/data/HumanGenome_Part100Update.gz' #-cacheArchive s3://fyp-input/HumanGenome.fa.gz#human  
    genome = readGenome(genomeFile)
    readSeq = readOptimisedReads(sys.stdin) #Human reads=28,094,847
    
    for read, quality in readSeq:
	#Human genome=64,185,939 lines -> 3,273,481,150 bytes
	offset, rqDict = alignKmer(read, quality, genome)
         
        #write results to STDOUT (standard output)
        for o in offset: #to remove empty list and '[' ']' characters
            #tab-delimited, key:offset of match with reads, value:<default count of 1, genome subsequence matched, read matched, corresponding quality> 
            print '%s\t%s\t%s\t%s\t%s' % (o, 1, genome[o:o+len(read)], read, quality) 
            #The output here will be the input for the reduce step 
 
if __name__ == '__main__':
    main() 
    