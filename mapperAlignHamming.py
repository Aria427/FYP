#!/usr/bin/env python
#This file contains the map step for the aligment MapReduce implementation.

import sys
import smart_open #Add bootstrap action to s3://fyp-input/installPythonModules.sh

#This function reads the genome in chunks (groups of bytes) from the file stored in S3.
def readGenomeEMR(s3FileUrl, bytesNum=100):
    with smart_open.smart_open(s3FileUrl, 'r') as dataFile:
        firstLine = dataFile.readline()
        if firstLine.startswith('>'):
            firstLine = ''
            pass #ignore header information
        data = firstLine + dataFile.read(bytesNum).rstrip().upper().replace('N', '').replace(' ', '')
        yield data

def readGenome(s3File):
    genome = ''
    with smart_open.smart_open(s3File, 'r') as f: 
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

#This function aligns one read (with its corresponding quality) to the genome using the Hamming distance approximate matching method.
def alignHamming(read, quality, genome):
    readQualityDictionary = {} #key:read, value:list of quality integers

    #maximum Hamming distance = 2
    reverseRead = reverseComplement(read)
    matchOffset = naiveApproxHamming(read, genome, 2) #check if read matches in forward direction of genome
    matchOffset.extend(naiveApproxHamming(reverseRead, genome, 2)) #add results of any matches in reverse complement of genome

    if len(list(matchOffset)) > 0: #match - read aligned in at least one place
        qualityQ = []
        for q in quality:
            qualityQ.append(phred33ToQ(q))
        readQualityDictionary[read] = qualityQ
    return matchOffset, readQualityDictionary
    
def main():
    #hard-coded reference genome stored in S3 via Amazon EMR
    genomeFile = 's3://fyp-input.genome/HumanGenome.fa.gz' #Frankfurt region doesn't work, Ireland does
    g = readGenome(genomeFile)   
    readSeq = readInputReads(sys.stdin) #Human reads=28,094,847
    
    #readsMatched, readsMismatched = 0, 0
    for read, quality in readSeq: 
        #Human genome=64,185,939 lines -> 3,273,481,150 bytes
        #genome = readGenomeEMR(genomeFile, 500000000) #500,000,000 bytes = 0.5GB
        #overlap = ''
        
        #for g in genome:
        #    g = overlap + g
        offset, rqDict = alignHamming(read, quality, g)
            
            #write results to STDOUT (standard output)
        for o in offset: #to remove empty list and '[' ']' characters
                #if o == '':
                #    readsMismatched += 1
                #else:
                #tab-delimited, key:offset of match with reads, value:<default count of 1, genome subsequence matched, read matched, corresponding quality> 
                print '%s\t%s\t%s\t%s' % (o, 1, g[o:o+60], read, quality)
                    #readsMatched += 1
                #The output here will be the input for the reduce step  
            
        #    overlap = g[-100:] #100 for PhiX, 60 for Human   
 
if __name__ == '__main__':
    main()            
