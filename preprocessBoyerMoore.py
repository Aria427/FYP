#!/usr/bin/env python
#This file contains the preprocessing of the reads as done by the Boyer Moore aligner.

from path import path
import mapperAlignBoyerMoore as bm
import pickle

#This function reads the sequencing reads after optimisation as input to the mapper.      
def readOptimisedReads(readsFile):
    sequence, quality = '', ''
    with open(readsFile, 'r') as file:
        while True: #runs until EOF
            line = file.readline() 
            if not line: #reached EOF
                 break
                
            line = line.split()
            sequence = line[0]
            quality = line[1]
    
            yield sequence, quality

readsFile = path('Data/HumanReads_100').abspath()
readSeq = readOptimisedReads(readsFile)

preprocFile = path('Preprocessed/BoyerMooreReads').abspath()
with open(preprocFile, 'w') as f:
    n = 5 #maximum number of mismatches
    for read, quality in readSeq:
        segmentLength = int(round(len(read) / (n+1)))
        allMatches = set() 
        
        for i in range(n+1): 
            start = i*segmentLength
            end = min((i+1)*segmentLength, len(read)) 
            
            bmObject = bm.BoyerMoore(read[start:end], alphabet='ACGT') #does preprocessing of Boyer-Moore
            pickle.dump(bmObject, f) #write object to file
           
#Not sure if this is needed yet!            
