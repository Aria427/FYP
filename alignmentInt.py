#!/usr/bin/env python
#This file includes functions for aligning reads with the reference genome by using the different matching algorithms using integers. 
   
import matchingInt

#The genome is double stranded and so the reads can come from one strand or the other.    
#To match both the read and the reverse complement of the read to the genome: 
def reverseComplement(read):
    reverseRead = read[::-1]
    return reverseRead

#This aligns using the Hamming distance approximate matching method:
def alignHamming(reads, genome):
    readsMatched = 0
    readsCount = 0
    readsOffsets = []
    print '1'
    for read in reads: #loop over generator
        print '2'
        #reverseRead = reverseComplement(read)
        matchOffsets = matchingInt.naiveApproxHamming(read, genome) #check if read matches in forward direction of genome 
        #matchOffsets.extend(matchingInt.naiveApproxHamming(reverseRead, genome)) #add results of any matches in reverse complement of genome
        readsCount += 1 
        if (readsCount % 100) == 0:
            print "*"
        if len(list(matchOffsets)) > 0: #match - read aligned in at least one place
            readsMatched += 1
        readsOffsets.append(matchOffsets)  
    return readsMatched, readsCount, readsOffsets

#The edit distance is not implemented using integers as there is no accurate method of 
#finding insertions and deletions as opposed to substitutions with 1s and 0s.
     
#This aligns using the k-mer indexing method:
def alignKmer(reads, genome):
    readsMatched = 0
    readsCount = 0
    readsOffsets = []
    for read in reads: 
        reverseRead = reverseComplement(read)
        index = matchingInt.kmerIndex(genome, 10)
        matchOffsets = matchingInt.queryKmerIndex(read, genome, index) #check if read matches in forward direction of genome
        matchOffsets.extend(matchingInt.queryKmerIndex(reverseRead, genome, index)) #add results of any matches in reverse complement of genome
        readsCount += 1
        if (readsCount % 100) == 0:
            print "*"
        if len(list(matchOffsets)) > 0: #match - read aligned in at least one place
            readsMatched += 1
        readsOffsets.append(matchOffsets) 
    return readsMatched, readsCount, readsOffsets 
   
#This aligns using the FM indexing method:
def alignFM(reads, genome): 
    readsMatched = 0
    readsCount = 0
    readsOffsets = []
    for read in reads: 
        reverseRead = reverseComplement(read)
        fm = matchingInt.fmIndex(genome)
        matchOffsets = fm.occurrences(read) #check if read matches in forward direction of genome
        matchOffsets.extend(fm.occurrences(reverseRead)) #add results of any matches in reverse complement of genome
        readsCount += 1
        if (readsCount % 100) == 0:
            print "*"
        if len(list(matchOffsets)) > 0: #match - read aligned in at least one place
            readsMatched += 1
        readsOffsets.append(matchOffsets) 
    return readsMatched, readsCount, readsOffsets  