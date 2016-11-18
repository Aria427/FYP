#!/usr/bin/env python
#This file includes functions for aligning reads with the reference genome by using the different matching algorithms using integers. 
    
import matching

#The genome is double stranded and so the reads can come from one strand or the other.    
#To match both the read and the reverse complement of the read to the genome: 
def reverseComplement(read):
    complement = {0:1, 1:0}
    for zero, one in complement.items():
        read = read.replace(zero, one)
    return read

#This aligns using the Hamming distance approximate matching method:
def alignHamming(reads, genome):
    readsMatched = 0
    readsCount = 0
    nextReads = next(reads)
    while nextReads: 
        #nextReverseReads = reverseComplement(nextReads)
        readsMatched += matching.naiveApproxHammingInt(nextReads, genome) #check if read matches in forward direction of genome
        #add results of any matches in reverse complement of genome
        readsCount += 1 
        if (readsCount % 50) == 0:
            print "*"
        try:
            nextReads = next(reads)
        except StopIteration:
            break
    return readsMatched, readsCount
