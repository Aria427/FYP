#!/usr/bin/env python
#This file contains the map step for the aligment MapReduce implementation.

import sys
import gzip

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

#This function calculates the edit distance of a read against the genome.
#Edit distance = minimum no of edits (substitutions, insertions, deletions) required to change one string into another.
def editDistance(pattern, text): #pattern=read, text=genome
    D = [] #distance matrix
    
    for i in xrange(len(pattern) + 1): #initialize D as x+1 by y+1 array of 0s 
        D.append([0] * (len(text) + 1))
    for i in xrange(len(pattern) + 1): #first row of ascending integers
        D[i][0] = i
    for i in xrange(len(text) + 1): #first column of ascending integers
        D[0][i] = i
    for i in xrange(1, len(pattern) + 1): #fill in rest of matrix row by row
        for j in xrange(1, len(text) + 1):
            distHorizontal = D[i][j-1] + 1
            distVertical = D[i-1][j] + 1
            if pattern[i-1] == text[j-1]:
                distDiagonal = D[i-1][j-1]
            else:
                distDiagonal = D[i-1][j-1] + 1
            D[i][j] = min(distHorizontal, distVertical, distDiagonal)
    
    return D[-1][-1] #edit distance = bottom-right                
                
#This function is a naive algorithm for approximate matching using the Edit distance.
#Takes too long => inefficient.
def naiveApproxEdit(pattern, text):
    matchOffsets = []
    editDist = editDistance(pattern, text)

    #loop through every position from where P could start without running past the end of T 
    for i in xrange(len(text) - len(pattern) + 1): #loop over all possible alignments of P in T from left to right	 
        mismatches = 0
        for j in xrange(len(pattern)): #loop over characters in P from left to right
            #i'th alignment and j'th character   
            if text[i+j] != pattern[j]: #mismatch
                mismatches += 1  
                if mismatches > editDist:   
                    break #exceeded maximum distance
                    
        if mismatches <= editDist: #approximate match
            matchOffsets.append(i)
            
    return matchOffsets                               

#This function aligns the reads to the genome using the Edit distance approximate matching method.
def alignEdit(read, quality, genome):
    readQualityDictionary = {} #key:read, value:list of quality integers

    matchOffset = naiveApproxEdit(read, genome) #check if read matches in forward/backward direction of genome
        
    if len(list(matchOffset)) > 0: #match - read aligned in at least one place
        readQualityDictionary[read] = quality

    return matchOffset, readQualityDictionary  
    
def main():
    genomeFile = '/home/aria427/test/data/HumanGenome_Part100Update.gz' #-cacheArchive s3://fyp-input/HumanGenome.fa.gz#human  
    genome = readGenome(genomeFile)
    readSeq = readOptimisedReads(sys.stdin) #Human reads=28,094,847
    
    for read, quality in readSeq:
	#Human genome=64,185,939 lines -> 3,273,481,150 bytes
	offset, rqDict = alignEdit(read, quality, genome)
         
        #write results to STDOUT (standard output)
        for o in offset: #to remove empty list and '[' ']' characters
            #tab-delimited, key:offset of match with reads, value:<default count of 1, genome subsequence matched, read matched, corresponding quality> 
            print '%s\t%s\t%s\t%s\t%s' % (o, 1, genome[o:o+len(read)], read, quality) 
            #The output here will be the input for the reduce step 
 
if __name__ == '__main__':
    main()  
    