#!/usr/bin/env python
#This file includes functions for aligning reads with the reference genome by using the different matching algorithms. 
    
import matching

#The genome is double stranded and so the reads can come from one strand or the other.    
#To match both the read and the reverse complement of the read to the genome: 
def reverseComplement(read):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'} #each base is associated with its complementary base
    reverseRead = ''
    for base in read:
        reverseRead = complement[base] + reverseRead #complement added to beginning in order to reverse the read from end to start
    return reverseRead

#This aligns using the Hamming distance approximate matching method:
def alignHamming(reads, genome):
    readsMatched = 0
    readsCount = 0
    readsOffsets = []
    nextReads = next(reads)
    for read in nextReads: 
        nextReads = nextReads[:50] #prefix of read as all 100 bases have a smaller chance of matching
        nextReverseReads = reverseComplement(nextReads)
        #encodedReads = huffmanCompression.encode(nextReads, output)
        #encodedReverseReads = huffmanCompression.encode(nextReverseReads, output)
        matchOffsets = matching.naiveApproxHamming(nextReads, genome) #check if read matches in forward direction of genome
        matchOffsets.extend(matching.naiveApproxHamming(nextReverseReads, genome)) #add results of any matches in reverse complement of genome
        readsCount += 1
        if (readsCount % 50) == 0:
            print "*"
        if len(list(matchOffsets)) > 0: #match - read aligned in at least one place
            readsMatched += 1
        readsOffsets.append(matchOffsets) 
        nextReads = next(reads)
    return readsMatched, readsCount, readsOffsets      

#This aligns using the Edit distance approximate matching method:
def alignEdit(reads, genome):
    readsMatched = 0
    readsCount = 0
    readsOffsets = []
    nextReads = next(reads)
    for read in nextReads: 
        nextReads = nextReads[:50] #prefix of read as all 100 bases have a smaller chance of matching
        nextReverseReads = reverseComplement(nextReads)
        #encodedReads = huffmanCompression.encode(nextReads, output)
        #encodedReverseReads = huffmanCompression.encode(nextReverseReads, output)
        matchOffsets = matching.approxEdit(nextReads, genome) #check if read matches in forward direction of genome
        matchOffsets.extend(matching.approxEdit(nextReverseReads, genome)) #add results of any matches in reverse complement of genome
        readsCount += 1
        if (readsCount % 50) == 0:
            print "*"
        if len(list(matchOffsets)) > 0: #match - read aligned in at least one place
            readsMatched += 1
        readsOffsets.append(matchOffsets) 
        nextReads = next(reads)
    return readsMatched, readsCount, readsOffsets  
  
#This aligns using the k-mer indexing method:
def alignKmer(reads, genome):
    readsMatched = 0
    readsCount = 0
    readsOffsets = []
    nextReads = next(reads)
    for read in nextReads: 
        nextReads = nextReads[:50] #prefix of read as all 100 bases have a smaller chance of matching
        nextReverseReads = reverseComplement(nextReads)
        #encodedReads = huffmanCompression.encode(nextReads, output)
        #encodedReverseReads = huffmanCompression.encode(nextReverseReads, output)
        index = matching.kmerIndex(genome, 10)
        matchOffsets = matching.queryKmerIndex(nextReads, genome, index) #check if read matches in forward direction of genome
        matchOffsets.extend(matching.queryKmerIndex(nextReverseReads, genome, index)) #add results of any matches in reverse complement of genome
        readsCount += 1
        if (readsCount % 50) == 0:
            print "*"
        if len(list(matchOffsets)) > 0: #match - read aligned in at least one place
            readsMatched += 1
        readsOffsets.append(matchOffsets) 
        nextReads = next(reads)
    return readsMatched, readsCount, readsOffsets  
   
#This aligns using the FM indexing method:
def alignFM(reads, genome):#, output):
    readsMatched = 0
    readsCount = 0
    readsOffsets = []
    nextReads = next(reads)
    for read in nextReads: 
        nextReads = nextReads[:50] #prefix of read as all 100 bases have a smaller chance of matching
        nextReverseReads = reverseComplement(nextReads)
        #encodedReads = huffmanCompression.encode(nextReads, output)
        #encodedReverseReads = huffmanCompression.encode(nextReverseReads, output)
        fm = matching.fmIndex(genome)
        matchOffsets = fm.occurrences(nextReads) #check if read matches in forward direction of genome
        matchOffsets.extend(fm.occurrences(nextReverseReads)) #add results of any matches in reverse complement of genome
        readsCount += 1
        if (readsCount % 50) == 0:
            print "*"
        if len(list(matchOffsets)) > 0: #match - read aligned in at least one place
            readsMatched += 1
        readsOffsets.append(matchOffsets) 
        nextReads = next(reads)
    return readsMatched, readsCount, readsOffsets  
 