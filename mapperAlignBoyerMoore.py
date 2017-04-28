#!/usr/bin/env python
#This file contains the map step for the aligment MapReduce implementation.

import sys
import gzip

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

#This function uses the Z algorithm to preprocess a pattern P -> returns its Z array.
#Zi(P) = length of longest substring of P that starts at i > 1 and matches a prefix of P.
def zArray(pattern):
    assert len(pattern) > 1
    z = [len(pattern)] + [0] * (len(pattern)-1)
    
    #initial comparison of pattern[1:] with prefix
    for i in range(1, len(pattern)):
        if pattern[i] == pattern[i-1]:
            z[1] += 1
        else:
            break
    
    r, l = 0, 0
    if z[1] > 0:
        r, l = z[1], 1
    
    for k in range(2, len(pattern)):
        assert z[k] == 0

        if k > r: #Case 1
            for i in range(k, len(pattern)):
                if pattern[i] == pattern[i-k]:
                    z[k] += 1
                else:
                    break
            r, l = k + z[k] - 1, k

        else: #Case 2
            #calculate length of beta
            nbeta = r - k + 1
            zkp = z[k - l]

            if nbeta > zkp: #Case 2a: Zkp wins
                z[k] = zkp

            else: #Case 2b: Compare characters just past r
                nmatch = 0
                for i in range(r+1, len(pattern)):
                    if pattern[i] == pattern[i - k]:
                        nmatch += 1
                    else:
                        break
                l, r = k, r + nmatch
                z[k] = r - k + 1

    return z

#This function generates the N array from the Z array.
def nArray(pattern):
    return zArray(pattern[::-1])[::-1]

#This function generates the L' array using the pattern and the N array.
# L'[i] = largest index j < n such that N[j] = |P[i:]|
def bigLprimeArray(pattern, n):
    lprime = [0] * len(pattern)
    
    for j in range(len(pattern)-1):
        i = len(pattern) - n[j]
        if i < len(pattern):
            lprime[i] = j + 1

    return lprime

#This function generates the L array using the pattern and the L' array.
#L[i] = largest index j < n such that N[j] >= |P[i:]|
def bigLarray(pattern, lprime):
    l = [0] * len(pattern)
    l[1] = lprime[1]

    for i in range(2, len(pattern)):
        l[i] = max(l[i-1], lprime[i])
        
    return l

#This function generates the l' array using the N array.
def smallLprimeArray(n):
    smallLprime = [0] * len(n)
    
    for i in range(len(n)):
        if n[i] == i+1:  #prefix matching a suffix
            smallLprime[len(n)-i-1] = i+1
    
    for i in range(len(n)-2, -1, -1):  #'smear' them out to the left
        if smallLprime[i] == 0:
            smallLprime[i] = smallLprime[i+1]

    return smallLprime

#This function returns the tables needed to apply the good suffix rule.
def goodSuffixTable(pattern):
    n = nArray(pattern)
    lprime = bigLprimeArray(pattern, n)
    
    return lprime, bigLarray(pattern, lprime), smallLprimeArray(n)

#This function returns the amount to shift as determined by the good suffix rule;
#given a mismatch at offset i and the L/L' and l' arrays.
def goodSuffixMismatch(i, bigLprime, smallLprime):
    length = len(bigLprime)
    assert i < length
    if i == length - 1:
        return 0
    i += 1  #i points to leftmost matching position of P
    
    if bigLprime[i] > 0:
        return length - bigLprime[i]
    return length - smallLprime[i]

#This function also returns the amount to shift as determined by the good suffix rule;
#but given a full match of P to T.
def goodSuffixMatch(smallLprime):
    return len(smallLprime) - smallLprime[1]

#This function creates and returns a dense bad character table indexed by offset and then by character;
#given a pattern string and a list with ordered alphabet characters.
def denseBadCharTable(pattern, amap):
    tab = []
    nxt = [0] * len(amap)
    
    for i in range(0, len(pattern)):
        c = pattern[i]
        assert c in amap
        tab.append(nxt[:])
        nxt[amap[c]] = i+1

    return tab

#This class encapsulates a pattern and its associated Boyer-Moore preprocessing. 
class BoyerMoore(object):
    #This function initialises the Boyer Moore object to be used.     
    def __init__(self, pattern, alphabet='ACGT'):
        self.pattern = pattern
        self.alphabet = alphabet
        
        self.amap = {}
        for i in range(len(self.alphabet)): 
            self.amap[self.alphabet[i]] = i #create map from alphabet characters to integers
        
        self.badChar = denseBadCharTable(pattern, self.amap) #create bad character rule table
        _, self.bigL, self.smallLprime = goodSuffixTable(pattern) #create good suffix rule table
    
    #This function returns the amount to shift as determined by the bad character rule at offset i.
    #Upon mismatch, skip alignments until mismatch becomes match or P moves past mismatched character.
    def badCharacterRule(self, i, c):
        assert c in self.amap
        ci = self.amap[c]
        assert i > (self.badChar[i][ci]-1)
        return i - (self.badChar[i][ci]-1)
    
    #This function returns the amount to shift as determined by the 'weak' good suffix rule given a mismatch at offset i.
    #Let t = substring matched by inner loop; skip until there are no mismtaches between P and t or P moves past t.
    def goodSuffixRule(self, i):
        length = len(self.bigL)
        assert i < length
        if i == length - 1:
            return 0
        i += 1  #i points to leftmost matching position of P
        
        if self.bigL[i] > 0:
            return length - self.bigL[i]
        return length - self.smallLprime[i]
    
    #This function returns the amount to shift for the case where P matches T exactly.
    #A special case of the good suffix rule.
    def matchSkip(self):
        return len(self.smallLprime) - self.smallLprime[1]

#This function implements the Boyer-Moore exact matching algorithm.
def boyerMooreExact(pattern, bmObject, text): #the Boyer-Moore object does the preprocessing for P
    i = 0 #keeps track of where we are in the text
    occurrences = []
    alignments = 0 #number of aligned characters
    characterComparisons = 0 #number of character comparisons

    #loop through all the positions in T where P could start w/o running past the end of T
    while i < len(text) - len(pattern) + 1: 
        shift = 1 
        mismatched = False
        alignments += 1
        
        for j in range(len(pattern)-1, -1, -1): #loop through P from end to start
            characterComparisons += 1
            if pattern[j] != text[i+j]: #mismatch
                bcShift = bmObject.badCharacterRule(j, text[i+j])
                gsShift = bmObject.goodSuffixRule(j)
                shift = max(shift, bcShift, gsShift) #largest shift will save most time and comparisons
                mismatched = True
                break
        
        if not mismatched: #matched T exactly
            occurrences.append(i)
            gsShift = bmObject.matchSkip()
            shift = max(shift, gsShift)
            
        i += shift #update position by shift
        
    return occurrences#, alignments, characterComparisons
    
#This function implements approximate matching on Boyer-Moore using the pigeonhole principle. 
def boyerMooreApproximate(pattern, text, n): #n = max number of mismatches
    segmentLength = int(round(len(pattern) / (n+1))) #n+1 segments where at least 1 must match perfectly against T
    allMatches = set() #all the indices where matches where found
    
    for i in range(n+1): #for each segment in P
        #bounds of P for segment being searched
        start = i*segmentLength
        end = min((i+1)*segmentLength, len(pattern)) #min() to not run past end of P
        
        bmObject = BoyerMoore(pattern[start:end], alphabet='ACGT') #does preprocessing of Boyer-Moore
        matches = boyerMooreExact(pattern[start:end], bmObject, text)
        
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
                
    return list(allMatches)               

#This function aligns the reads to the genome using the Boyer-Moore approximate algorithm.
def alignBoyerMoore(read, quality, genome):
    readQualityDictionary = {} #key:read, value:list of quality integers
    
    #maximum number of mismatches = 5
    matchOffset = boyerMooreApproximate(read, genome, 5) #check if read matches in forward/backward direction of genome
        
    if len(list(matchOffset)) > 0: #match - read aligned in at least one place
        readQualityDictionary[read] = quality
    
    return matchOffset, readQualityDictionary 
    
def main():
    genomeFile = '/home/aria427/test/data/HumanGenome_Part100Update.gz' #-cacheArchive s3://fyp-input/HumanGenome.fa.gz#human  
    genome = readGenome(genomeFile)
    readSeq = readOptimisedReads(sys.stdin) #Human reads=28,094,847
    
    for read, quality in readSeq:
	#Human genome=64,185,939 lines -> 3,273,481,150 bytes
	offset, rqDict = alignBoyerMoore(read, quality, genome)
         
        #write results to STDOUT (standard output)
        for o in offset: #to remove empty list and '[' ']' characters
            #tab-delimited, key:offset of match with reads, value:<default count of 1, genome subsequence matched, read matched, corresponding quality> 
            print '%s\t%s\t%s\t%s\t%s' % (o, 1, genome[o:o+len(read)], read, quality) 
            #The output here will be the input for the reduce step 
  
 
if __name__ == '__main__':
    main()
    