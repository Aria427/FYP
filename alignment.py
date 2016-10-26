"""
This file includes functions for naive exact & approximate matching along with the alignment implementation.
"""
import bisect
import sys
import numpy

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
#The following is also an online naive algorithm but for approximate matching using the Hamming distance:
def naiveApproxHamming(pattern, text, maxHammingDist=1):
    matchOffsets = []
    for i in range(len(text) - len(pattern) + 1):
        mismatches = 0
        for j in range(len(pattern)):
            if text[i+j] != pattern[j]: #mismatch
                mismatches += 1  
                if mismatches > maxHammingDist:   
                    break           #exceeded maximum distance
        if mismatches <= maxHammingDist: #approximate match
            matchOffsets.append(i)
    return matchOffsets

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

#The following is also an online naive algorithm for approximate matching but using the Edit distance:
def naiveApproxEdit(pattern, text, maxEditDist=5500):
    matchOffsets = []
    D = []
    D = editDistance(pattern, text)
    if D > maxEditDist: #exceeded maximum distance
        return []
    if D <= maxEditDist: #approximate match
        matchOffsets.append(D)
    return matchOffsets      

#Edit distance backtrace
def editDistanceTrace(pattern, text, D):
    i = len(pattern)
    j = len(text) 
    while i > 0:
        distDiagonal, distVertical, distHorizontal = sys.maxint, sys.maxint, sys.maxint
        delta = None
        if i > 0 and j > 0:
            delta = 0 if pattern[i-1] == text[j-1] else 1
            distDiagonal = D[i-1, j-1] + delta
        if i > 0:
            distVertical = D[i-1, j] + 1
        if j > 0:
            distHorizontal = D[i, j-1] + 1
        if distDiagonal <= distVertical and distDiagonal <= distHorizontal: #diagonal wins
            i -= 1; j -= 1
        elif distVertical <= distHorizontal: #vertical win - insertion in P wrt T
            i -= 1
        else: #horizontal wins
            j -= 1
    return j #offset of leftmost character of T in match

#The following is an approximate matching algorithm using the backtrace of the edit distance:
def approxEdit(pattern, text): #if multiple alignments tie for best, report leftmost
    matchOffsets = []
    distanceJ = None
    distance = len(pattern) + len(text) 
    D = numpy.zeros((len(pattern)+1, len(text)+1), dtype=int)
    D[1:, 0] = range(1, len(pattern)+1) #first row = 0s, first column = usual
    for i in xrange(1, len(pattern)+1):
        for j in xrange(1, len(text)+1):
            delta = 1 if pattern[i-1] != text[j-1] else 0
            distHorizontal = D[i][j-1] + 1
            distVertical = D[i-1][j] + 1
            distDiagonal = D[i-1][j-1] + delta
            D[i, j] = min(distHorizontal, distVertical, distDiagonal) 
    for j in xrange(0, len(text)+1):
        if D[len(pattern), j] < distance:
            distanceJ = j
            distance = D[len(pattern), j] #minimum edit distance in last row
    matchOffset = editDistanceTrace(pattern, text[:distanceJ], D) #backtrace stops as it gets to first row
    matchOffsets.append(matchOffset) 
    return matchOffsets  
    
#The following is an index object which is applied for k-mer indexing:    
class Index(object):
    def __init__(self, text, k): #initialise index from all k-mer length substrings - preprocesses string text
        self.k = k  #k-mer length
        self.index = []
        for i in range(len(text) - k + 1): #for each k-mer 
            self.index.append((text[i:i+k], i)) #add (k-mer, offset) tuple
        self.index.sort()  #sort in ascending order according to k-mer
    
    def query(self, pattern): #returns no of index hits for first k-mer of P
        kmer = pattern[:self.k]  #query with first k-mer
        i = bisect.bisect_left(self.index, (kmer, -1)) #binary search for 1st position in list where k-mer occurs
        indexHits = [] #all indices in T where the first k bases of P match
        while i < len(self.index): #collect matching index entries
            if self.index[i][0] != kmer:
                break
            indexHits.append(self.index[i][1]) #append 2nd value of tuple
            i += 1
        return indexHits

#The following is the offline algorithm using k-mer indexing by the above class:
def queryIndex(pattern, text, index):
    k = index.k #length of k
    matchOffsets = []
    for i in index.query(pattern): #returns a list of possible places where P could start (where 1st k bases of P match 1st k bases of T)
        if pattern[k:] == text[i+k:i+len(pattern)]:  #verify rest of P matches
            matchOffsets.append(i)
    return matchOffsets


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
        #matchOffsets = approxEdit(nextReads, genome) #check if read matches in forward direction of genome
        #matchOffsets.extend(approxEdit(reverseComplement(nextReads), genome)) #add results of any matches in reverse complement of genome
        readsCount += 1
        if (readsCount % 50) == 0:
            print "*"
        if len(list(matchOffsets)) > 0: #match - read aligned in at least one place
            readsMatched += 1
        readsOffsets.append(matchOffsets) 
        nextReads = next(reads)
    return readsMatched, readsCount, readsOffsets, readsMatches    
 