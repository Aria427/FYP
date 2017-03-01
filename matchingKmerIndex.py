#!/usr/bin/env python
#This file includes functions for the k-mer indexing algorithm.

import bisect #allows binary search on a list

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

#This function implements the k-mer indexing exact offline algorithm.
def kmerIndexExact(pattern, text, indexObject):
    k = indexObject.k #length of k
    matchOffsets = []

    for i in indexObject.query(pattern): #returns a list of possible places where P could start (where 1st k bases of P match 1st k bases of T)
        if pattern[k:] == text[i+k:i+len(pattern)]:  #verify rest of P matches
            matchOffsets.append(i)
            
    return matchOffsets
 
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
                
    print 'index hits = {', indexHits, '}'
    return list(allMatches)    
    