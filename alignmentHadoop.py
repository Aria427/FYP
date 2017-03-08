#!/usr/bin/env python
#This file includes functions for aligning reads with the reference genome via Hadoop.
    
import matchingDistances as dist
import matchingBoyerMoore as bm
import matchingKmerIndex as kIdx
import matchingFmIndex as fmIdx
import matchingSmithWaterman as sw
import matchingBurrowsWheeler as bwa

#This function finds the reverse complement of a sequencing read.
#The genome is double stranded, so the reads can come from one strand or the other.    
def reverseComplement(read):
    complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'} #each base is associated with its complementary base
    reverseRead = ''
    for base in read:
        reverseRead = complement[base] + reverseRead #complement added to beginning in order to reverse the read from end to start
    return reverseRead

#Quality values are ASCII-encoded -> a character encodes an integer.
#Low quality score Q => low confidence in value.
#High quality score Q => high confidence in value.
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
    matchOffset = dist.naiveApproxHamming(read, genome, 2) #check if read matches in forward direction of genome
    matchOffset.extend(dist.naiveApproxHamming(reverseRead, genome, 2)) #add results of any matches in reverse complement of genome

    if len(list(matchOffset)) > 0: #match - read aligned in at least one place
        qualityQ = []
        for q in quality:
            qualityQ.append(phred33ToQ(q))
        readQualityDictionary[read] = qualityQ
    return matchOffset, readQualityDictionary
    
#This function aligns the reads to the genome using the Edit distance approximate matching method.
def alignEdit(read, quality, genome):
    readQualityDictionary = {} #key:read, value:list of quality integers

    reverseRead = reverseComplement(read)
    matchOffset = dist.naiveApproxEdit(read, genome) #check if read matches in forward direction of genome
    matchOffset.extend(dist.naiveApproxEdit(reverseRead, genome)) #add results of any matches in reverse complement of genome
        
    if len(list(matchOffset)) > 0: #match - read aligned in at least one place
        qualityQ = []
        for q in quality:
            qualityQ.append(phred33ToQ(q))
        readQualityDictionary[read] = qualityQ

    return matchOffset, readQualityDictionary  
  
#This function aligns the reads to the genome using the Boyer-Moore approximate algorithm.
def alignBoyerMoore(read, quality, genome):
    readQualityDictionary = {} #key:read, value:list of quality integers
    
    #maximum number of mismatches = 2
    reverseRead = reverseComplement(read)
    matchOffset = bm.boyerMooreApproximate(read, genome, 2) #check if read matches in forward direction of genome
    matchOffset.extend(bm.boyerMooreApproximate(reverseRead, genome, 2)) #add results of any matches in reverse complement of genome
        
    if len(list(matchOffset)) > 0: #match - read aligned in at least one place
        qualityQ = []
        for q in quality:
            qualityQ.append(phred33ToQ(q))
        readQualityDictionary[read] = qualityQ
    
    return matchOffset, readQualityDictionary     
    
#This function aligns the reads to the genome using the k-mer approximate indexing method.
def alignKmer(read, quality, genome):
    readQualityDictionary = {} #key:read, value:list of quality integers

    index = kIdx.KmerIndex(genome, 10) #k-mer of length 10
    
    #maximum number of mismatches = 2
    reverseRead = reverseComplement(read)
    matchOffset = kIdx.kmerIndexApproximate(read, genome, 2, index) #check if read matches in forward direction of genome
    matchOffset.extend(kIdx.kmerIndexApproximate(reverseRead, genome, 2, index)) #add results of any matches in reverse complement of genome
        
    if len(list(matchOffset)) > 0: #match - read aligned in at least one place
        qualityQ = []
        for q in quality:
            qualityQ.append(phred33ToQ(q))
        readQualityDictionary[read] = qualityQ
     
    return matchOffset, readQualityDictionary  
   
#This function aligns the reads to the genome using the FM indexing method.
def alignFM(read, quality, genome):
    readQualityDictionary = {} #key:read, value:list of quality integers
    
    fm = fmIdx.fmIndex(genome)
    reverseRead = reverseComplement(read)
    matchOffset = fm.occurrences(read) #check if read matches in forward direction of genome
    matchOffset.extend(fm.occurrences(reverseRead)) #add results of any matches in reverse complement of genome
        
    if len(list(matchOffset)) > 0: #match - read aligned in at least one place
        qualityQ = []
        for q in quality:
            qualityQ.append(phred33ToQ(q))
        readQualityDictionary[read] = qualityQ
        
    return matchOffset, readQualityDictionary  
    
#This function aligns the reads to the genome using the Smith Waterman local alignment algorithm.
def alignSmithWaterman(read, quality, genome):
    readQualityDictionary = {} #key:read, value:list of quality integers
    
    #maximum number of mismatches = 2
    reverseRead = reverseComplement(read)
    matchOffset = sw.smithWatermanApproximate(genome, read, 2) #check if read matches in forward direction of genome
    matchOffset.extend(sw.smithWatermanApproximate(genome, reverseRead, 2)) #add results of any matches in reverse complement of genome
        
    if len(list(matchOffset)) > 0: #match - read aligned in at least one place
        qualityQ = []
        for q in quality:
            qualityQ.append(phred33ToQ(q))
        readQualityDictionary[read] = qualityQ

    return matchOffset, readQualityDictionary   
    
#This function aligns the reads to the genome using the Burrows Wheeler approximate algorithm.
def alignBurrowsWheeler(read, quality, genome):
    readQualityDictionary = {} #key:read, value:list of quality integers
    
    #sa = bwa.suffixArray(genome)
    bw = bwa.bwt(genome) 
    bwr = bwa.bwt(genome[::-1]) 
    
    #maximum number of mismatches = 2
    reverseRead = reverseComplement(read)
    matchOffset = bwa.burrowsWheelerApproximate(bw, bwr, read, 2) #check if read matches in forward direction of genome
    matchOffset.extend(bwa.burrowsWheelerApproximate(bw, bwr, reverseRead, 2)) #add results of any matches in reverse complement of genome
        
    if len(list(matchOffset)) > 0: #match - read aligned in at least one place
        qualityQ = []
        for q in quality:
            qualityQ.append(phred33ToQ(q))
        readQualityDictionary[read] = qualityQ
        
    return matchOffset, readQualityDictionary      
