#!/usr/bin/env python
#This file includes functions for the Burrows-Wheeler alignment algorithm.

alphabet = ['A', 'C', 'G', 'T']

#This function inefficiently generates the suffix array for some text.
def suffixArray(text):
    sa = sorted([(text[i:], i) for i in xrange(0, len(text)+1)]) #sorted list of all rotations
    return map(lambda x: x[1], sa) #extract offsets from last column

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

    for char in alphabet:
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
def burrowsWheelerExact(bwt, text):
    ranks, totals = rank(bwt)
    fc = firstColumn(totals)

    if text[-1] not in fc:
        return 0 #character not in bwt

    l, r = fc[text[-1]]
    i = len(text)-2

    while i >= 0 and r > 1:
        char = text[i]
        l = fc[char][0] + ranks[char][l-1] #R(aW)
        r = fc[char][0] + ranks[char][r-1] #Rbar(aW)
        i -= 1
        print('l: '+str(l)+' r: '+str(r))
        #print(''.join(bw[l:r]))
        #print(''.join(sorted(bw)[l:r]))

    return r-l #number of exact matches of text to bwt

seq1 = 'ATAGACGACATACAGACAGCATACAGACAGCATACAGA' #reference
seq2 = 'TTTAGCATGCGCATATCAGCAATACAGACAGATACG' #read
    
sa = suffixArray(seq1)
bw = bwt(seq1)
bwr = bwt(seq1[::-1])

burrowsWheelerExact(bw, seq2)    
    