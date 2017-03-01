#!/usr/bin/env python
#This file includes functions for aligning reads with the reference genome by using the different matching algorithms using strings. 
    
import matchingKmerIndex as kIdx

import matchingString

#This function finds the reverse complement of a sequencing read.
#The genome is double stranded, so the reads can come from one strand or the other.    
def reverseComplement(read):
    complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'} #each base is associated with its complementary base
    reverseRead = ''
    for base in read:
        reverseRead = complement[base] + reverseRead #complement added to beginning in order to reverse the read from end to start
    return reverseRead

#This function aligns the reads to the genome using the Hamming distance approximate matching method.
def alignHamming(reads, genome):
    readsMatched = 0
    readsCount = 0
    readsOffsets = []
    for read in reads: 
        reverseRead = reverseComplement(read)
        matchOffsets = matchingString.naiveApproxHamming(read, genome) #check if read matches in forward direction of genome
        matchOffsets.extend(matchingString.naiveApproxHamming(reverseRead, genome)) #add results of any matches in reverse complement of genome

        if len(list(matchOffsets)) > 0: #match - read aligned in at least one place
            readsMatched += 1
        readsOffsets.append(matchOffsets) 
        
        readsCount += 1
        if (readsCount % 1) == 0:
            print '*'
    return readsMatched, readsCount, readsOffsets      

#This function aligns the reads to the genome using the Edit distance approximate matching method.
def alignEdit(reads, genome):
    readsMatched = 0
    readsCount = 0
    readsOffsets = []
    for read in reads: 
        reverseRead = reverseComplement(read)
        matchOffsets = matchingString.approxEdit(read, genome) #check if read matches in forward direction of genome
        matchOffsets.extend(matchingString.approxEdit(reverseRead, genome)) #add results of any matches in reverse complement of genome
        
        if len(list(matchOffsets)) > 0: #match - read aligned in at least one place
            readsMatched += 1
        readsOffsets.append(matchOffsets) 
        
        readsCount += 1
        if (readsCount % 100) == 0:
            print '*'
    return readsMatched, readsCount, readsOffsets  
  
#This function aligns the reads to the genome using the k-mer indexing method.
def alignKmer(reads, genome):
    readsMatched = 0
    readsCount = 0
    readsOffsets = []
    index = kIdx.KmerIndex(genome, 10)
    for read in reads: 
        reverseRead = reverseComplement(read)
        matchOffsets = kIdx.kmerIndexExact(read, genome, index) #check if read matches in forward direction of genome
        matchOffsets.extend(kIdx.kmerIndexExact(reverseRead, genome, index)) #add results of any matches in reverse complement of genome
        
        readsCount += 1
        if (readsCount % 100) == 0:
            print '*'
        if len(list(matchOffsets)) > 0: #match - read aligned in at least one place
            readsMatched += 1
        readsOffsets.append(matchOffsets) 
    return readsMatched, readsCount, readsOffsets  
   
#This function aligns the reads to the genome using the FM indexing method.
def alignFM(reads, genome):#, output):
    readsMatched = 0
    readsCount = 0
    readsOffsets = []
    fm = matchingString.fmIndex(genome)
    for read in reads: 
        reverseRead = reverseComplement(read)
        matchOffsets = fm.occurrences(read) #check if read matches in forward direction of genome
        matchOffsets.extend(fm.occurrences(reverseRead)) #add results of any matches in reverse complement of genome
        
        readsCount += 1
        if (readsCount % 1000000) == 0:
            print '*'
        if len(list(matchOffsets)) > 0: #match - read aligned in at least one place
            readsMatched += 1
        readsOffsets.append(matchOffsets) 
    return readsMatched, readsCount, readsOffsets  
 