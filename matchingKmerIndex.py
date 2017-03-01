#!/usr/bin/env python
#This file includes functions for the k-mer indexing algorithm.

import bisect #allows binary search on a list

#This class is an index object which is applied for k-mer indexing. 
class kmerIndex(object):
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

#This function is the offline algorithm using k-mer indexing by the above class.
def queryKmerIndex(pattern, text, index):
    k = index.k #length of k
    matchOffsets = []

    for i in index.query(pattern): #returns a list of possible places where P could start (where 1st k bases of P match 1st k bases of T)
        if pattern[k:] == text[i+k:i+len(pattern)]:  #verify rest of P matches
            matchOffsets.append(i)
            
    return matchOffsets
    