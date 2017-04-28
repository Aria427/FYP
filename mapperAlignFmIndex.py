#!/usr/bin/env python
#This file contains the map step for the aligment MapReduce implementation.

import sys
import gzip
from itertools import chain, islice
import pickle

#This function reads the genome in chunks.
def readGenomeChunks(s3File, bytesNum=100):
    with gzip.open(s3File, 'r') as f:
        f.readline() #ignore header line with genome information
        for chunk in iter(lambda: f.read(bytesNum), ''):
            data = chunk.rstrip().upper().replace('N', '').replace('\n', '').replace(' ', '')
            yield data

#This function reads the whole genome from the file.
def readGenome(s3File):
    genome = ''
    with gzip.open(s3File, 'r') as f: 
        for line in f:
            if line and line[0] != '>': #ignore header line with genome information
                l = line.rstrip().upper().replace('N', '').replace(' ', '')
                genome += l
        return genome
                  
#This function reads the sequencing reads after optimisation as input to the mapper.      
def readOptimisedReads(file):
    sequence, quality = '', ''
    while True: #runs until EOF
    	line = file.readline() 
        if not line: #reached EOF
             break
            
        line = line.split()
        sequence = line[0]
        quality = line[1]

        yield sequence, quality

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

#This function aligns the reads to the genome using the FM indexing method.
def alignFM(read, quality, index):
    readQualityDictionary = {} #key:read, value:list of quality integers
    
    matchOffset = index.occurrences(read) #check if read matches in forward/backward direction of genome
        
    if len(list(matchOffset)) > 0: #match - read aligned in at least one place
        readQualityDictionary[read] = quality
        
    return matchOffset, readQualityDictionary 
    
def main():
    genomeFile = '/home/aria427/test/data/HumanGenome_Part100Update.gz' #-cacheArchive s3://fyp-input/HumanGenome.fa.gz#human  
    preprocGenome = '/home/aria427/test/preprocessed/FmIndexGenome'
    
    genome = readGenome(genomeFile)
    readSeq = readOptimisedReads(sys.stdin) #Human reads=28,094,847
    
    f = open(preprocGenome, 'r')
    index = pickle.load(f)
    
    for read, quality in readSeq:
        #Human genome=64,185,939 lines -> 3,273,481,150 bytes
        offset, rqDict = alignFM(read, quality, index)
         
        #write results to STDOUT (standard output)
        for o in offset: #to remove empty list and '[' ']' characters
            #tab-delimited, key:offset of match with reads, value:<default count of 1, genome subsequence matched, read matched, corresponding quality> 
            print '%s\t%s\t%s\t%s\t%s' % (o, 1, genome[o:o+len(read)], read, quality) 
            #The output here will be the input for the reduce step 
 
    f.close()    
            
if __name__ == '__main__':
    main()  
    