#!/usr/bin/env python
#This file includes functions for the Hamming and Edit distance algorithms.

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

    #loop through every position from where P could start without running past the end of T 
    for i in xrange(len(text) - len(pattern) + 1): #loop over all possible alignments of P in T from left to right	 
        mismatches = 0
        for j in xrange(len(pattern)): #loop over characters in P from left to right
            #i'th alignment and j'th character   
            if text[i+j] != pattern[j]: #mismatch
                mismatches += 1  
                if mismatches > maxHammingDist:   
                    break #exceeded maximum distance
                    
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
#Takes too long => inefficient.
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

#This function is a naive algorithm for approximate matching using the Edit distance.
#Takes too long => inefficient.
def naiveApproxEdit(pattern, text):
    matchOffsets = []
    editDist = editDistance(pattern, text)

    #loop through every position from where P could start without running past the end of T 
    for i in xrange(len(text) - len(pattern) + 1): #loop over all possible alignments of P in T from left to right	 
        mismatches = 0
        for j in xrange(len(pattern)): #loop over characters in P from left to right
            #i'th alignment and j'th character   
            if text[i+j] != pattern[j]: #mismatch
                mismatches += 1  
                if mismatches > editDist:   
                    break #exceeded maximum distance
                    
        if mismatches <= editDist: #approximate match
            matchOffsets.append(i)
            
    return matchOffsets
    