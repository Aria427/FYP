#!/usr/bin/env python
#This file includes functions for exact & approximate matching of reads with the reference genome using integers.

import fileParsing

import bisect
import sys
import numpy as np
#import gmpy

def bitLength(integer): #length of 1s and 0s
    length = 0
    while integer:
        length += 1
        integer >>= 1
    return length

def bitCount(integer): #length of 1s only
    count = 0
    while integer:
        count += (integer & 1)
        integer >>= 1
    return count 

def ffs(x): #find first set = index of LSB (set => 1) 0 bit
    return (x&-x).bit_length()-1

#def hammingWeight(x, y):
#    return gmpy.popcount(x ^ y)
 
#Hamming distance = minimum no of substitutions required to change one string into another.
#The following is also an online naive algorithm but for approximate matching using the Hamming distance:
def naiveApproxHamming(pattern, text, maxHammingDist=10):
    matchOffsets = []
    mismatches = 0 #hamming distance - no of differing bits
    for i in xrange(len(pattern)):
        bitwise = pattern[i] ^ text #bitwise exclusive or
        #mismatches = hammingWeight(pattern[i], text)
        #Hamming weight (no of non-zero bits) found using Wegner algorithm
        while np.all(bitwise != 0): #bit is set => mismatch
            mismatches += 1
            #index = ffs(bitwise)
            bitwise &= bitwise - 1 #clear lowest order non-zero bit
            if mismatches > maxHammingDist: #exceeded maximum distance
                break
        if mismatches <= maxHammingDist:
            matchOffsets.append(i)
    return matchOffsets

#This might not be possible as k-mer indexing was constructed for strings.    
#The following is an index object which is applied for k-mer indexing:    
class kmerIndex(object):
    def __init__(self, text, k): #initialise index from all k-mer length substrings - preprocesses string text
        self.k = k  #k-mer length
        self.index = []
        for i in xrange(bitCount(text) - k + 1): #for each k-mer 
            self.index.append((text & i+k, i)) #add (k-mer, offset) tuple
        self.index.sort()  #sort in ascending order according to k-mer
    
    def query(self, pattern): #returns no of index hits for first k-mer of P
        kmer = pattern >> self.k  #query with first k-mer
        i = bisect.bisect_left(self.index, (kmer, -1)) #binary search for 1st position in list where k-mer occurs
        indexHits = [] #all indices in T where the first k bases of P match
        while i < len(self.index): #collect matching index entries
            if self.index[i][0] != kmer:
                break
            indexHits.append(self.index[i][1]) #append 2nd value of tuple
            i += 1
        return indexHits

#The following is the offline algorithm using k-mer indexing by the above class:
def queryKmerIndex(pattern, text, index):
    k = index.k #length of k
    matchOffsets = []
    for i in index.query(pattern): #returns a list of possible places where P could start (where 1st k bases of P match 1st k bases of T)
        if ((pattern & k) == (text & i+k)):# == 0): #verify rest of P matches
            matchOffsets.append(i)
    return matchOffsets

#The following generates the suffix array for some text:
def suffixArray(text):
    #sorted() => simple but inefficient
    satups = sorted([(text[i:], i) for i in xrange(0, len(text))]) #sorted list of all rotations
    return map(lambda x: x[1], satups) #extract offsets from last column

#The following generates the Burrows-Wheeler Transform of some text using the suffix array:
def bwtSa(text, sa=None):
    bwt = []
    dollarRow = None
    if sa is None:
        sa = suffixArray(text)
    for si in sa: #take integers just to the left of sorted suffixes
        if si == 0:
            dollarRow = len(bwt)
            bwt.append('$') 
        else:
            bwt.append(text[si-1])
    return (''.join(bwt), dollarRow) #integer-ized version of list bw   

#This might not be possible as FM indexing, suffix arrays and BWT were constructed for strings.
#The following is an object which evaluates the rank checkpoints and the rank queries for FM indexing:
#evaluation = O(1) time; checkpoints = O(m) space where m = length of T
class fmCheckpoints(object): 
    def __init__(self, bwt, spacing=4): #create checkpoints periodically while scanning BWT
        self.checkpts = {}
        self.spacing = spacing #spacing between checkpoints
        total = {}
        for char in bwt: #for each unique character in T
            if char not in total:
                total[char] = 0 #entry in tally dictionary
                self.checkpts[char] = [] #entry in checkpoint map
        for i in xrange(0, len(bwt)):
            total[bwt[i]] += 1 #up to and including
            if (i % spacing) == 0:
                for char in total.iterkeys():
                    self.checkpts[char].append(total[char]) #construct checkpoints
    
    def rank(self, bwt, char, row): #characters in BWT up to and including row - ascending B-rank
        if row < 0 or char not in self.checkpts:
            return 0
        i, delta = row, 0
        while (i % self.spacing) != 0: #walk up to left to calculate rank
            if bwt[i] == char:
                delta += 1
            i -= 1
        return self.checkpts[char][i // self.spacing] + delta #parallel list of ranks
    
#The following is an index object which is applied for FM indexing:
#index = O(m) size where m = length of T; checkpoints & SA samples = O(1) spacing
#queries = O(n) where n = query length; search of k occurrences = O(n+k)
class fmIndex():        
    def __init__(self, text, spacing=4, ssaIval=4):
        text = format(text, 'b')
        text = fileParsing.binaryToBase(text)
        if (text[-1]) != '$':
            text += '$' #add $ if not present
        sa = suffixArray(text)
        self.sa = sa
        self.bwt, self.dollarRow = bwtSa(text, sa) #BWT string and 2 offset
        self.slen = len(self.bwt)
        self.checkpts = fmCheckpoints(self.bwt, spacing) #rank checkpoints
        self.first = {}
        occurrences = dict()
        occCount = 0
        for char in self.bwt:
            occurrences[char] = occurrences.get(char, 0) + 1 #occurrences of every character
        for char, count in sorted(occurrences.iteritems()):
            self.first[char] = occCount #compact view of first column
            occCount += count
 
    def count(self, char): #occurrences of characters < char
        if char not in self.first: #char does not occur in T (rare)
            for c in sorted(self.first.iterkeys()):
                if char < c: 
                    return self.first[c]
            return self.first[c]
        else:
            return self.first[char]        
            
    def interval(self, prefix): #range of BWT rows
        l, r = 0, self.slen - 1 #closed (inclusive) interval
        for i in xrange(len(prefix)-1, -1, -1): #from right to left
            l = self.checkpts.rank(self.bwt, prefix[i], l-1) + self.count(prefix[i])
            r = self.checkpts.rank(self.bwt, prefix[i], r) + self.count(prefix[i]) - 1
            if r < l:
                break
        return l, r+1
    
    def resolve(self, row): #offset of BWT row wrt T
        steps = 0
        def stepLeft(row): #respective of character in BWT row, move left
            char = self.bwt[row]
            return self.checkpts.rank(self.bwt, char, row-1) + self.count(char)
        print '3'
        while row not in self.sa:
            row = stepLeft(row)
            steps += 1
        return self.sa[row] + steps
    
    def occurrences(self, pattern): #offsets of all occurrences of P
        pattern = format(pattern, 'b')
        pattern = fileParsing.binaryToBase(pattern)
        l, r = self.interval(pattern)
        return [ self.resolve(x) for x in xrange(l, r) ]       
