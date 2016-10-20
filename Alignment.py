"""
This file includes functions for naive exact & approximate matching along with the alignment implementation.
"""

#The following is a naive algorithm for exact matching where all occurrences are recorded:
def naiveExact(pattern, text):
    matchOffsets = []
    #loop through every position from where P could start without running past the end of T 
    for i in range(len(text) - len(pattern) + 1): #loop over all possible alignments of P in T from left to right	 
        match = True
        for j in range(len(pattern)): #loop over characters in P from left to right
            #i'th alignment and j'th character            
            if text[i+j] != pattern[j]: #compare characters of P with T
                match = False #mismatch; reject alignment
                break
        if match:
            matchOffsets.append(i)  #all characters matched; record
    return matchOffsets

#Hamming distance = minimum no of substitutions required to change one string into another.
#The following is also a naive algorithm but for approximate matching using the Hamming distance:
def naiveApproxHamming(pattern, text, maxHammingDist=1):
    matchOffsets = []
    matches = []
    for i in range(len(text) - len(pattern) + 1):
        mismatches = 0
        for j in range(len(pattern)):
            if text[i+j] != pattern[j]: #mismatch
                mismatches += 1     
                if mismatches > maxHammingDist:
                    break           #exceeded maximum distance
        if mismatches <= maxHammingDist: #approximate match
            matchOffsets.append(i)
            matches.append(pattern)
    return matchOffsets, matches

def approxMatchOffsets(pattern, text):
    matchOffsets, _ = naiveApproxHamming(pattern, text)    
    return matchOffsets
    
def approxMatches(pattern, text):
    _, matches = naiveApproxHamming(pattern, text)
    return matches

#The genome is double stranded and so the reads can come from one strand or the other.    
#To match both the read and the reverse complement of the read to the genome: 
def reverseComplement(read):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'} #each base is associated with its complementary base
    reverseRead = ''
    for base in read:
        reverseRead = complement[base] + reverseRead #complement added to beginning in order to reverse the read from end to start
    return reverseRead

#To align the reads against the genome to see how many match:
def align(reads, genome):
    readsMatched = 0
    readsCount = 0
    readsOffsets = []
    readsMatches = []
    nextReads = next(reads)
    for read in nextReads: 
        nextReads = nextReads[:50] #prefix of read as all 100 bases have a smaller chance of matching
        matchOffsets = approxMatchOffsets(nextReads, genome) #check if read matches in forward direction of genome
        matchOffsets.extend(approxMatchOffsets(reverseComplement(nextReads), genome)) #add results of any matches in reverse complement of genome
        matches = approxMatches(read, genome)
        matches.extend(approxMatches(reverseComplement(read), genome))
        readsCount += 1
        if len(matchOffsets) > 0: #match - read aligned in at least one place
            readsMatched += 1
        readsOffsets.append(matchOffsets) 
        readsMatches.append(matches)
        nextReads = next(reads)
    return readsMatched, readsCount, readsOffsets, readsMatches
    