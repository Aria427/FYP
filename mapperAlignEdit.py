#!/usr/bin/env python
#This file contains the map step for the aligment MapReduce implementation.

import sys
import gzip

#This function reads the genome in chunks (groups of bytes) from the file stored in S3.
def readGenomeChunks(s3File, bytesNum=100):
    with gzip.open(s3File, 'r') as f:
        f.readline()
        for chunk in iter(lambda: f.read(bytesNum), ''):
            data = chunk.rstrip().upper().replace('N', '').replace('\n', '').replace(' ', '')
            yield data

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

#This function aligns the reads to the genome using the Edit distance approximate matching method.
def alignEdit(read, quality, genome):
    readQualityDictionary = {} #key:read, value:list of quality integers

    reverseRead = reverseComplement(read)
    matchOffset = naiveApproxEdit(read, genome) #check if read matches in forward direction of genome
    matchOffset.extend(naiveApproxEdit(reverseRead, genome)) #add results of any matches in reverse complement of genome
        
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
            offset, rqDict = alignEdit(read, quality, g)
         
            #write results to STDOUT (standard output)
            for o in offset: #to remove empty list and '[' ']' characters
                #tab-delimited, key:offset of match with reads, value:<default count of 1, genome subsequence matched, read matched, corresponding quality> 
                print '%s\t%s\t%s\t%s\t%s' % (o+filesOffset, 1, g[o:o+len(read)], read, quality) 
                #The output here will be the input for the reduce step  
            
            overlap = g[-99:] #100-1 for PhiX, 60-1 for Human read => -1 as last 60 have already been read
            filesOffset += (len(g)-100) #store offset according to overlap as file is read in chunks  
 
if __name__ == '__main__':
    main()            
