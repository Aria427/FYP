#!/usr/bin/env python
#This file contains the reduce step for the alignment MapReduce implementation.       

import alignmentString
import sys

#This function reads in the output from the combiner using a generator.
def readCombinerOutput(file, separator='\t'):
    for line in file:
        yield line.rstrip().split(separator, 1)

#This function iterates through the genome lines.
def genomeLine(data):
    for isGR, grValue in data:
        if isGR == 'G':
            for g in grValue:
                genSeq = g[0] 
                yield genSeq   
        else: 
            pass

#This function iterates through the read lines.               
def readLine(data):
    for isGR, grValue in data:
        if isGR == 'R':
            for r in grValue:
                readSeq = r[1] 
                yield readSeq  
        else:
            pass
        
def main():
    #input comes from STDIN (standard input)
    data = readCombinerOutput(sys.stdin) 
  
    genome = genomeLine(data)
    reads = readLine(data)
    
    totalMatches, totalCount, totalOffsets = 0, 0, []
    for g in genome:
         matchesCount, count, offsets = alignmentString.alignFM(reads, g)
            
         totalMatches += matchesCount
         totalCount = count #number of reads stays the same as every chunk goes through each read again
         totalOffsets.append(offsets)
         print '%d/%d reads matched the genome.' % (totalMatches, totalCount) #write result to STDOUT

if __name__ == '__main__':
    main()

