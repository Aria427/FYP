#!/usr/bin/env python
#This file includes functions for aligning reads with the reference genome by using the different matching algorithms using integers. 
    
import matchingInt

#The genome is double stranded and so the reads can come from one strand or the other.    
#To match both the read and the reverse complement of the read to the genome: 
def reverseComplement(read):
    reverseRead = ~read
    return reverseRead

#This aligns using the Hamming distance approximate matching method:
def alignHamming(reads, genome):
    readsMatched = 0
    readsCount = 0
    readsOffsets = []
    print '1'
    nextReads = next(reads) #program is halting on this statement
    print '2'
    for read in nextReads: #TypeError: 'long' object is not iterable 
        print '3'
        nextReverseReads = reverseComplement(nextReads)
        print '4'
        matchOffsets = matchingInt.naiveApproxHamming(nextReads, genome) #check if read matches in forward direction of genome
        print '5'    
        matchOffsets.extend(matchingInt.naiveApproxHamming(nextReverseReads, genome)) #add results of any matches in reverse complement of genome
        readsCount += 1 
        if (readsCount % 50) == 0:
            print "*"
        if len(list(matchOffsets)) > 0: #match - read aligned in at least one place
            readsMatched += 1
        readsOffsets.append(matchOffsets)     
        try:
            nextReads = next(reads)
        except StopIteration:
            break
    return readsMatched, readsCount, readsOffsets

#The edit distance is not implemented using integers as there is no accurate method of 
#finding insertions and deletions as opposed to substitutions with 1s and 0s.
     
#This aligns using the k-mer indexing method:
def alignKmer(reads, genome):
    readsMatched = 0
    readsCount = 0
    readsOffsets = []
    nextReads = next(reads)
    while nextReads: 
        nextReverseReads = reverseComplement(nextReads)
        index = matchingInt.kmerIndex(genome, 10)
        matchOffsets = matchingInt.queryKmerIndex(nextReads, genome, index) #check if read matches in forward direction of genome
        matchOffsets.extend(matchingInt.queryKmerIndex(nextReverseReads, genome, index)) #add results of any matches in reverse complement of genome
        readsCount += 1
        if (readsCount % 50) == 0:
            print "*"
        if len(list(matchOffsets)) > 0: #match - read aligned in at least one place
            readsMatched += 1
        readsOffsets.append(matchOffsets) 
        try:
            nextReads = next(reads)
        except StopIteration:
            break
    return readsMatched, readsCount, readsOffsets 
   
#This aligns using the FM indexing method:
def alignFM(reads, genome): 
    readsMatched = 0
    readsCount = 0
    readsOffsets = []
    nextReads = next(reads)
    while nextReads: 
        nextReverseReads = reverseComplement(nextReads)
        fm = matchingInt.fmIndex(genome)
        matchOffsets = fm.occurrences(nextReads) #check if read matches in forward direction of genome
        matchOffsets.extend(fm.occurrences(nextReverseReads)) #add results of any matches in reverse complement of genome
        readsCount += 1
        if (readsCount % 50) == 0:
            print "*"
        if len(list(matchOffsets)) > 0: #match - read aligned in at least one place
            readsMatched += 1
        readsOffsets.append(matchOffsets) 
        try:
            nextReads = next(reads)
        except StopIteration:
            break
    return readsMatched, readsCount, readsOffsets  