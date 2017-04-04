#!/usr/bin/env python
#This file contains the map step for the aligment MapReduce implementation.

import sys
import gzip
import bisect

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

#This function aligns the reads to the genome using the k-mer approximate indexing method.
def alignKmer(read, quality, genome):
    readQualityDictionary = {} #key:read, value:list of quality integers

    index = KmerIndex(genome, 10) #k-mer of length 10
    
    #maximum number of mismatches = 2
    reverseRead = reverseComplement(read)
    matchOffset = kmerIndexApproximate(read, genome, 2, index) #check if read matches in forward direction of genome
    matchOffset.extend(kmerIndexApproximate(reverseRead, genome, 2, index)) #add results of any matches in reverse complement of genome
        
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
            offset, rqDict = alignKmer(read, quality, g)
         
            #write results to STDOUT (standard output)
            for o in offset: #to remove empty list and '[' ']' characters
                #tab-delimited, key:offset of match with reads, value:<default count of 1, genome subsequence matched, read matched, corresponding quality> 
                print '%s\t%s\t%s\t%s\t%s' % (o+filesOffset, 1, g[o:o+len(read)], read, quality) 
                #The output here will be the input for the reduce step  
            
            overlap = g[-99:] #100-1 for PhiX, 60-1 for Human read => -1 as last 60 have already been read
            filesOffset += (len(g)-100) #store offset according to overlap as file is read in chunks  
 
if __name__ == '__main__':
    main()            
