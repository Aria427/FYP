"""
This file includes functions for naive exact & approximate matching along with the alignment implementation.
"""
import bisect
import sys

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
    possibleMatch = ''
    for i in range(len(text) - len(pattern) + 1):
        mismatches = 0
        for j in range(len(pattern)):
            if text[i+j] != pattern[j]: #mismatch
                mismatches += 1  
                if mismatches > maxHammingDist:
                    possibleMatch = ''
                    break           #exceeded maximum distance
            possibleMatch = ''.join(pattern[j])
        if mismatches <= maxHammingDist: #approximate match
            matchOffsets.append(i)
            matches.append(possibleMatch)
            possibleMatch = ''
    return matchOffsets#, matches

#Edit distance = minimum no of edits (substitutions, insertions, deletions) required to change one string into another.
def editDistance(pattern, text):
    D = [] #distance matrix
    for i in range(len(pattern) + 1): #initialize D as x+1 by y+1 array of 0s 
        D.append([0] * (len(text) + 1))
    for i in range(len(pattern) + 1): #first row of ascending integers
        D[i][0] = i
    for i in range(len(text) + 1): #first column of ascending integers
        D[0][i] = i
    for i in range(1, len(pattern) + 1): #fill in rest of matrix row by row
        for j in range(1, len(text) + 1):
            distHorizontal = D[i][j-1] + 1
            distVertical = D[i-1][j] + 1
            if pattern[i-1] == text[j-1]:
                distDiagonal = D[i-1][j-1]
            else:
                distDiagonal = D[i-1][j-1] + 1
            D[i][j] = min(distHorizontal, distVertical, distDiagonal)
    
    return D[-1][-1] #edit distance = bottom-right

#The following is also a naive algorithm but for approximate matching using the Edit distance:
def naiveApproxEdit(pattern, text, maxEditDist=5500):
    matchOffsets = []
    D = []
    D = editDistance(pattern, text)
    if D > maxEditDist: #exceeded maximum distance
        return []
    if D <= maxEditDist: #approximate match
        matchOffsets.append(D)
    return matchOffsets      

class Index(object):
    def __init__(self, text, k): #initialise index from all k-mer length substrings
        self.k = k  # k-mer length
        self.index = []
        for i in range(len(text) - k + 1): #for each k-mer
            self.index.append((text[i:i+k], i)) #add (k-mer, offset) tuple
        self.index.sort()  #sort in ascending order according to k-mer
    
    def query(self, pattern): #returns no of index hits for first k-mer of P
        kmer = pattern[:self.k]  #query with first k-mer
        i = bisect.bisect_left(self.index, (kmer, -1)) #binary search
        hits = []
        while i < len(self.index): #collect matching index entries
            if self.index[i][0] != kmer:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits

def queryIndex(pattern, text, index):
    k = index.k
    matchOffsets = []
    for i in index.query(pattern):
        if pattern[k:] == text[i+k:i+len(pattern)]:  #verify rest of P matches
            matchOffsets.append(i)
    return matchOffsets
  
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
        index = Index(genome, 10)
        matchOffsets = queryIndex(nextReads, genome, index)
        matchOffsets.extend(queryIndex(reverseComplement(nextReads), genome, index))
        #matchOffsets = naiveApproxHamming(nextReads, genome) #check if read matches in forward direction of genome
        #matchOffsets.extend(naiveApproxHamming(reverseComplement(nextReads), genome)) #add results of any matches in reverse complement of genome
        #matches = approxMatches(read, genome)
        #matches.extend(approxMatches(reverseComplement(read), genome))
        readsCount += 1
        if (readsCount % 50) == 0:
            print "*"
        if len(list(matchOffsets)) > 0: #match - read aligned in at least one place
            readsMatched += 1
        readsOffsets.append(matchOffsets) 
        #readsMatches.append(matches)
        nextReads = next(reads)
    return readsMatched, readsCount, readsOffsets, readsMatches    
 