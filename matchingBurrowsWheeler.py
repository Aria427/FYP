#!/usr/bin/env python
#This file includes functions for the Burrows-Wheeler alignment algorithm.

from itertools import chain, islice
from operator import itemgetter

ALPHABET = ['A', 'C', 'G', 'T']
O = {} #dictionary with keys = $, A, C, G, T and values = arrays of counts
C = {} #number of lexographically greater symbols in the reference genome
D = [] # estimated lower bounds of mismatches in the short read

#rewards/penalties => score
GAP = 1
MISMATCH = 1
MATCH = 0

prunes = 0

#This function inefficiently generates the suffix array for some text.
def suffixarray(text):
    sa = sorted([(text[i:], i) for i in xrange(0, len(text)+1)]) #sorted list of all rotations
    return map(lambda x: x[1], sa) #extract offsets from last column

#This function generates the suffix array for some text.
#It implements the algorithm of Vladu and Negruseri:
#   http://web.stanford.edu/class/cs97si/suffix-array.pdf
def suffixArray(text):
    L = sorted((b, i) for i, b in enumerate(text)) #list of pairs = (base, index)
    n = len(text)+1
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
def bwt(text):
    bwt = []

    for i in suffixArray(text):
        if i == 0:
            bwt.append('$')
        else:
            bwt.append(text[i-1])
    
    return ''.join(bwt)

#This function generates a dictionary of characters to range of occurrences in the first column.
def firstColumn(totals):
    col = {}
    temp = 0
    
    for i, j in sorted(totals.iteritems()):
        col[i] = (temp, temp+j)
        temp += j
        
    return col

#This function generates a list of the number of occurences of a character for each reference subsequence.
def rank(bwt):
    ranks, totals = {}, {}

    for char in ALPHABET:
        if (char not in totals) and (char != '$'):
            ranks[char] = []
            totals[char] = 0

    for char in bwt:
        if char != '$':
            totals[char] += 1
        for t in totals.iterkeys():
            ranks[t].append(totals[t])

    return ranks, totals
     
#This function decoded the Burrows-Wheeler Transform.
def ibwt(bwt):
    ranks, totals = rank(bwt)
    fc = firstColumn(totals)
    row = 0
    t = '$'

    while bwt[row] != '$':
        char = bwt[row]
        t = char+t
        row = fc[char][0] + ranks[row]

    return t

#This function implements the Burrows Wheeler exact matching algorithm.
def burrowsWheelerExact(bwt, pattern):
    ranks, totals = rank(bwt)
    fc = firstColumn(totals)

    if pattern[-1] not in fc:
        return 0 #character not in bwt

    l, r = fc[pattern[-1]]
    i = len(pattern)-2

    while i >= 0 and r > 1:
        char = pattern[i]
        l = fc[char][0] + ranks[char][l-1] #R(aW)
        r = fc[char][0] + ranks[char][r-1] #Rbar(aW)
        i -= 1
        #print('l: '+str(l)+' r: '+str(r))
        #print(''.join(bw[l:r]))
        #print(''.join(sorted(bw)[l:r]))

    return r-l #number of exact matches of text to bwt 

#This function calculates the number of lexographically greater symbols C in the reference genome.
def computeC(totals):
    C = {'A': 0, 'C': 0, 'G': 0, 'T': 0}

    for k in ALPHABET:
        for ref in ALPHABET:
            if ref < k: C[k] += totals[ref]

    return C

#This function calculates the estimated lower bounds of mismatches in the short read.
def computeD(pattern, C, Oprime, bw):
    k = 1
    l = len(bw)-2
    z = 0
    D = [0] * len(pattern)

    for i in range(0, len(pattern)):
        k = C[pattern[i]] + Oprime[pattern[i]][k-1] + 1
        l = C[pattern[i]] + Oprime[pattern[i]][l]
        if k > l:
            k = 1
            l = len(bw)-1
            z = z+1
        D[i] = z

    return D

#This function ensures that if D[i] = -1, its value will be considered as 0.
def getD(i):
    if i < 0:
        return 0
    else:
        return D[i]

#This function operates in a similar manner to getD().
def getO(char, index):
    if index < 0:
        return 0
    else:
        return O[char][index]

#This function implements the Burrows Wheeler approximate matching algorithm using recursion.
def burrowsWheelerApproximateRecursion(pattern, i, maxMismatches, k, l):    
    global prunes
    saIndex = set() #set of suffix array indices at which a match starts 
    
    #pruning based on estimated mistakes
    if maxMismatches < getD(i):
        prunes += 1
        return set()

    #end of query condition
    temp = set()
    if i < 0:
        for j in range(k,l+1):
            temp.add((j, maxMismatches))
        return temp
    
    #Insertion
    saIndex = saIndex.union(burrowsWheelerApproximateRecursion(
                            pattern, i-1, maxMismatches-GAP, k, l))

    for char in ALPHABET:
        tempk = C[char] + getO(char, k-1) + 1
        templ = C[char] + getO(char, l)

        if tempk <= templ:
            #Deletion
            saIndex = saIndex.union(burrowsWheelerApproximateRecursion(
                                    pattern, i, maxMismatches-GAP, tempk, templ))
            
            if char == pattern[i]: #match
                saIndex = saIndex.union(burrowsWheelerApproximateRecursion(
                                        pattern, i-1, maxMismatches+MATCH, tempk, templ))
                
            else: #mismatch
               saIndex = saIndex.union(burrowsWheelerApproximateRecursion(
                                        pattern, i-1, maxMismatches-MISMATCH, tempk, templ))

    return saIndex

#This function implements the Burrows Wheeler approximate matching algorithm.
#bw = bwt of genome, bwr = bwt of reverse genome, pattern = short read to match, maxMismatches = max number of differences 
def burrowsWheelerApproximate(bw, bwr, pattern, maxMismatches):    
    global prunes
    global O
    global C
    global D
    
    O, totals = rank(bw)
    Oprime, _ = rank(bwr) #reverse ranks
    C = computeC(totals) 
    D = computeD(pattern, C, Oprime, bw)

    saIndexSet = burrowsWheelerApproximateRecursion(pattern, len(pattern)-1, maxMismatches, 0, len(bw)-1)
    indexDict = {}

    for (i, j) in saIndexSet:
        if i in indexDict: #index already exists
            if indexDict[i] < j:
                indexDict[i] = j #pick higher maxMismatches value
                prunes += 1
        else:
            indexDict[i] = j

    return sorted(indexDict.items(), key=itemgetter(1), reverse=True) #sort by maxMismtaches in descending order 

#This function finds the best match position.
#sa = suffix array
def bestMatchPosition(bw, bwr, pattern, maxMismatches, sa):
    saIndexList = burrowsWheelerApproximate(bw, bwr, pattern, maxMismatches)
    
    if len(saIndexList) != 0:
        bestIndex, score = saIndexList[0]
        return sa[bestIndex]+1, score
    else:
        return -1, -1    
    
#This function displays the results.
def displayOutput(saIndexList, sa, s, read):
    saValues = [(sa[i], j) for (i, j) in saIndexList]

    print '\n-------------------------------------'
    print 'Reference: ' + s
    print 'Read: \t   ' + read + '\n'

    print str(len(saValues)) + ' match/es found!\n'
    print 'Score\tPosition\tSuffix\n'
    for v, x in saValues:
        print str(x) + '\t' + str(v) + '\t\t' + s[v:v+35]

    print '-------------------------------------'

"""    
ref = 'ACGTACGTACGTAAACCCGGGTTTACGT' #reference
read = 'ACGTAACCGGTTACGTAAGGTT' #read
    
sa = suffixArray(ref)
bw = bwt(ref) 
bwr = bwt(ref[::-1]) 

print burrowsWheelerExact(bw, read)     
    
threshold = 8 #error score upper bound
displayOutput(burrowsWheelerApproximate(bw, bwr, read, threshold), sa, ref, read)

threshold = 9 #error score upper bound
displayOutput(burrowsWheelerApproximate(bw, bwr, read, threshold), sa, ref, read)

threshold = 10 #error score upper bound
displayOutput(burrowsWheelerApproximate(bw, bwr, read, threshold), sa, ref, read)
"""