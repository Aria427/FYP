#!/usr/bin/env python
#This file includes functions for the FM indexing algorithm.

from itertools import chain, islice

#This function inefficiently generates the suffix array for some text.
def suffixarray(text):
    #sorted() => simple but inefficient
    satups = sorted([(text[i:], i) for i in xrange(0, len(text))]) #sorted list of all rotations
    return map(lambda x: x[1], satups) #extract offsets from last column

#This function generates the suffix array for some text.
#It implements the algorithm of Vladu and Negruseri:
#   http://web.stanford.edu/class/cs97si/suffix-array.pdf
def suffixArray(text):
    L = sorted((b, i) for i, b in enumerate(text)) #list of pairs = (base, index)
    n = len(text)
    count = 1
    
    while count < n:
        U = [0] * n
        for (r, i), (s, j) in zip(L, islice(L, 1, None)):
            U[j] = U[i] + (r != s)
        #Invariant: U[i] 
        #   index of text[i:i+count] in sorted list of the text's unique substrings of length count

        L = sorted(chain((((U[i],  U[i+count]), i) for i in range(n - count)),
                         (((U[i], -1), i) for i in range(n - count, n))))
        #Invariant: L[i][1]
        #   starting index in text of substring i of length count in sorted order
        
        count *= 2
        
    return [i for _, i in L]

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
    def __init__(self, bwt, spacing=5): #create checkpoints periodically while scanning BWT
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
    def __init__(self, text, spacing=5, ssaIval=5):
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
