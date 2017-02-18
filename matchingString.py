#!/usr/bin/env python
#This file includes functions for exact & approximate matching of reads with the reference genome using strings.

import bisect
import sys
import numpy

#This function is a naive algorithm for exact matching where all occurrences are recorded.
def naiveExact(pattern, text):
    matchOffsets = []
    #loop through every position from where P could start without running past the end of T 
    for i in xrange(len(text) - len(pattern) + 1): #loop over all possible alignments of P in T from left to right	 
        match = True
        for j in xrange(len(pattern)): #loop over characters in P from left to right
            #i'th alignment and j'th character            
            if text[i+j] != pattern[j]: #compare characters of P with T
                match = False #mismatch; reject alignment
                break
        if match:
            matchOffsets.append(i)  #all characters matched; record
    return matchOffsets

#Hamming distance = minimum no of substitutions required to change one string into another.
#This function is also an online naive algorithm but for approximate matching using the Hamming distance.
def naiveApproxHamming(pattern, text, maxHammingDist=1):
    matchOffsets = []
    for i in xrange(len(text) - len(pattern) + 1):
        mismatches = 0
        for j in xrange(len(pattern)):
            if text[i+j] != pattern[j]: #mismatch
                mismatches += 1  
                if mismatches > maxHammingDist:   
                    break           #exceeded maximum distance
        if mismatches <= maxHammingDist: #approximate match
            matchOffsets.append(i)
    return matchOffsets  
   
#This function calculates the edit distance of a read against the genome.
#Edit distance = minimum no of edits (substitutions, insertions, deletions) required to change one string into another.
def editDistance(pattern, text): #pattern=read, text=genome
    D = [] #distance matrix
    for i in xrange(len(pattern) + 1): #initialize D as x+1 by y+1 array of 0s 
        D.append([0] * (len(text) + 1))
    for i in xrange(len(pattern) + 1): #first row of ascending integers
        D[i][0] = i
    for i in xrange(len(text) + 1): #first column of ascending integers
        D[0][i] = i
    for i in xrange(1, len(pattern) + 1): #fill in rest of matrix row by row
        for j in xrange(1, len(text) + 1):
            distHorizontal = D[i][j-1] + 1
            distVertical = D[i-1][j] + 1
            if pattern[i-1] == text[j-1]:
                distDiagonal = D[i-1][j-1]
            else:
                distDiagonal = D[i-1][j-1] + 1
            D[i][j] = min(distHorizontal, distVertical, distDiagonal)
    
    return D[-1][-1] #edit distance = bottom-right

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

#This function is an approximate matching algorithm using the backtrace of the edit distance.
def approxEdit(pattern, text): #if multiple alignments tie for best, report leftmost
    matchOffsets = []
    distanceJ = None
    distance = len(pattern) + len(text) 
    D = numpy.zeros((len(pattern)+1, len(text)+1), dtype=int)
    D[1:, 0] = xrange(1, len(pattern)+1) #first row = 0s, first column = usual
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
    
#This class is an index object which is applied for k-mer indexing:. 
class kmerIndex(object):
    #This function initialises the object.
    def __init__(self, text, k): #initialise index from all k-mer length substrings - preprocesses string text
        self.k = k  #k-mer length
        self.index = []
        for i in xrange(len(text) - k + 1): #for each k-mer 
            self.index.append((text[i:i+k], i)) #add (k-mer, offset) tuple
        self.index.sort()  #sort in ascending order according to k-mer
    
    #This function finds the number of inex hits for first k-mer of pattern P
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

#This function generates the suffix array for some text.
def suffixArray(text):
    #sorted() => simple but inefficient
    satups = sorted([(text[i:], i) for i in xrange(0, len(text))]) #sorted list of all rotations
    return map(lambda x: x[1], satups) #extract offsets from last column

#This function generates the Burrows-Wheeler Transform of some text using the suffix array.
def bwtSa(text, sa=None):
    bwt = []
    dollarRow = None
    if sa is None:
        sa = suffixArray(text)
    for si in sa: #take characters just to the left of sorted suffixes
        if si == 0:
            dollarRow = len(bwt)
            bwt.append('$')
        else:
            bwt.append(text[si-1])
    return (''.join(bwt), dollarRow) #string-ized version of list bw   

#This class is an object which evaluates the rank checkpoints and the rank queries for FM indexing:
#evaluation = O(1) time; checkpoints = O(m) space where m = length of T
class fmCheckpoints(object): 
    #This function initialises the object.
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
    
    #This function ranks the characters in BWT up to and including row - ascending B-rank
    def rank(self, bwt, char, row):
        if row < 0 or char not in self.checkpts:
            return 0
        i, delta = row, 0
        while (i % self.spacing) != 0: #walk up to left to calculate rank
            if bwt[i] == char:
                delta += 1
            i -= 1
        return self.checkpts[char][i // self.spacing] + delta #parallel list of ranks
    
#This class is an index object which is applied for FM indexing:
#index = O(m) size where m = length of T; checkpoints & SA samples = O(1) spacing
#queries = O(n) where n = query length; search of k occurrences = O(n+k)
class fmIndex(): 
    #This function initialises the object.   
    def __init__(self, text, spacing=4, ssaIval=4):
        if text[-1] != '$':
            text += '$' #add $ if not present
        sa = suffixArray(text)
        self.sa = sa
        self.bwt, self.dollarRow = bwtSa(text, sa) #BWT string and $ offset
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
    
    #This function calculates the number of occurences of characters < char
    def count(self, char): 
        if char not in self.first: #char does not occur in T (rare)
            for c in sorted(self.first.iterkeys()):
                if char < c: 
                    return self.first[c]
            return self.first[c]
        else:
            return self.first[char]
    
    #This function generates the range of BWT rows.
    def interval(self, prefix):
        l, r = 0, self.slen - 1 #closed (inclusive) interval
        for i in xrange(len(prefix)-1, -1, -1): #from right to left
            l = self.checkpts.rank(self.bwt, prefix[i], l-1) + self.count(prefix[i])
            r = self.checkpts.rank(self.bwt, prefix[i], r) + self.count(prefix[i]) - 1
            if r < l:
                break
        return l, r+1
    
    #This function generates the offset of the BWT row wrt text T
    def resolve(self, row): 
        steps = 0
        def stepLeft(row): #respective of character in BWT row, move left
            char = self.bwt[row]
            return self.checkpts.rank(self.bwt, char, row-1) + self.count(char)
        while row not in self.sa:
            row = stepLeft(row)
            steps += 1
        return self.sa[row] + steps
    
    #This function generates the offsets of all occurrences of pattern P
    def occurrences(self, pattern): 
        l, r = self.interval(pattern)
        return [ self.resolve(x) for x in xrange(l, r) ]       
