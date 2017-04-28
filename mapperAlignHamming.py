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

#This function is also an online naive algorithm but for approximate matching using the Hamming distance.
def naiveApproxHamming(pattern, text, maxHammingDist=1):
    matchOffsets = []

    #loop through every position from where P could start without running past the end of T 
    for i in xrange(len(text) - len(pattern) + 1): #loop over all possible alignments of P in T from left to right	 
        mismatches = 0
        for j in xrange(len(pattern)): #loop over characters in P from left to right
            #i'th alignment and j'th character   
            if text[i+j] != pattern[j]: #mismatch
                mismatches += 1  
                if mismatches > maxHammingDist:   
                    break #exceeded maximum distance
                    
        if mismatches <= maxHammingDist: #approximate match
            matchOffsets.append(i)
            
    return matchOffsets                 

#This function aligns one read (with its corresponding quality) to the genome using the Hamming distance approximate matching method.
def alignHamming(read, quality, genome):
    readQualityDictionary = {} #key:read, value:list of quality integers

    #maximum Hamming distance = 5
    matchOffset = naiveApproxHamming(read, genome, 5) #check if read matches in forward/backward direction of genome

    if len(list(matchOffset)) > 0: #match - read aligned in at least one place
        readQualityDictionary[read] = quality

    return matchOffset, readQualityDictionary
    
def main():
    genomeFile = '/home/aria427/test/data/HumanGenome_Part100Update.gz' #-cacheArchive s3://fyp-input/HumanGenome.fa.gz#human  
    genome = readGenome(genomeFile)
    readSeq = readOptimisedReads(sys.stdin) #Human reads=28,094,847
    
    for read, quality in readSeq:
	#Human genome=64,185,939 lines -> 3,273,481,150 bytes
	offset, rqDict = alignHamming(read, quality, genome)
         
        #write results to STDOUT (standard output)
        for o in offset: #to remove empty list and '[' ']' characters
            #tab-delimited, key:offset of match with reads, value:<default count of 1, genome subsequence matched, read matched, corresponding quality> 
            print '%s\t%s\t%s\t%s\t%s' % (o, 1, genome[o:o+len(read)], read, quality) 
            #The output here will be the input for the reduce step  
    """        
    for read, quality in readSeq: 
        #Human genome=64,185,939 lines -> 3,273,481,150 bytes
        genome = readGenomeChunks(genomeFile, 250000) #250,000 bytes = 0.23842MB
        overlap = '' #size of read-1
        filesOffset = 0 #file is split in chunks so offset needs to change according to chunk
        
        for g in genome:
            g = overlap + g
            offset, rqDict = alignHamming(read, quality, g)
         
            #write results to STDOUT (standard output)
            for o in offset: #to remove empty list and '[' ']' characters
                #tab-delimited, key:offset of match with reads, value:<default count of 1, genome subsequence matched, read matched, corresponding quality> 
                print '%s\t%s\t%s\t%s\t%s' % (o+filesOffset, 1, g[o:o+len(read)], read, quality) 
                #The output here will be the input for the reduce step  
            
            overlap = g[-59:] #100-1 for PhiX, 60-1 for Human read => -1 as last 60 have already been read
            filesOffset += (len(g)-60) #store offset according to overlap as file is read in chunks
    """

if __name__ == '__main__':
    main()
    