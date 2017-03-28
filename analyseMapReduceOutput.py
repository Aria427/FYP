#!/usr/bin/env python
#This file includes functions for analysing the MapReduce alignment output.

import fileinput
from glob import glob
from path import path
from itertools import groupby
from operator import itemgetter
from collections import defaultdict, Counter

#This function reads in the output from the reducer using a generator.
def readReducerOutput(file, separator='\t'):
    with open(file, 'r') as f:
        for line in f:
            yield line.rstrip().split(separator, 3)  
            
#This function reads in the output from the reducers using a generator.
def readReducerOutputs(files, separator='\t'):
    fnames = glob(files)
    for line in fileinput.input(fnames):
        yield line.rstrip().split(separator, 3)
 
resultFiles = path('Output Test Files\PhiXHamming-output\part-000**').abspath() #* for multiple files
#resultFile = path('Output Test Files\PhiXHamming-output\part-00000').abspath()
data = readReducerOutputs(resultFiles)

offsetCountDict, readQualityDict, readQualityList = {}, defaultdict(list), []
for currentOffset, group in groupby(data, itemgetter(0)):
    try:
        for offset, count, rq in group:	
            offsetCountDict[int(offset)] = int(count)
            readQuality = rq.rstrip().replace('{', '').replace('}', '').split(':', 1)
            readQualityDict[readQuality[0]].append(readQuality[1])
            readQuality = rq.rstrip().replace('{', '').replace('}', '').split(',')
            readQualityList.append(readQuality)
    except ValueError:
        pass 

#1. At which offset do reads match the most?
#print offsetCountDict   
maxCount = max(v for k, v in offsetCountDict.iteritems())  
countOffsetDict = dict((v, k) for k, v in offsetCountDict.iteritems())
offsetMostMatches = countOffsetDict[maxCount]    
print 'The sequencing reads matched the most at offset %d.' % offsetMostMatches

#2. How many reads matched respective of quality?
flattenList = sum(readQualityList, [])
totalMatches = len(set(flattenList))
print 'The total number of reads matched = %d.' % totalMatches

#3. How many reads mismatched respective of quality?
totalReads = 1000 #PhiX=1000, Human=28,094,847
totalMismatches = totalReads - totalMatches
print 'The total number of reads mismatched = %d.' % totalMismatches

#4. Which read matched the most respective of quality?
rqFrequency = Counter(flattenList).most_common(2)
print 'The sequencing read which matched the most is %s.' % rqFrequency

