#!/usr/bin/env python
#This file contains the reduce step for the alignment MapReduce implementation.       

#import alignmentString
import sys

#This function reads in the output from the mapper using a generator.
def readMapperOutput(file, separator='\t'):
    for line in file:
        yield line.rstrip()#.split(separator, 1)

#This function iterates through the genome lines.
def genomeLine(data):
    for isGR, grValue in data:
        if isGR == 'G':
            #for g in grValue:
            #    genSeq = g[0] 
            genSeq = grValue[0]
            yield genSeq   
        else: 
            pass

#This function iterates through the read lines.               
def readLine(data):
    for isGR, grValue in data:
        if isGR == 'R':
            #for r in grValue:
            #    readSeq = r[0] 
            readSeq = grValue[0]
            yield readSeq  
        else:
            pass
        
def main():
    #input comes from STDIN (standard input)
    data = readMapperOutput(sys.stdin) 
  
    for d in data:
        print d      
    
if __name__ == '__main__':
    main()

"""
def main():
    #input comes from STDIN (standard input)
    data = readMapperOutput(sys.stdin) 
  
    genome = genomeLine(data)  

    totalMatches, totalCount, totalOffsets = 0, 0, []
    overlap = ''
    for g in genome:
         g = overlap + g #overlap is appended to the start of the next chunk
         reads = readLine(data) 
         matchesCount, count, offsets = alignmentString.alignFM(reads, g)
            
         totalMatches += matchesCount
         totalCount = count #number of reads stays the same as every chunk goes through each read again
         totalOffsets.append(offsets)
         overlap = g[-100:] #100 for PhiX, 60 for Human
    print '%d/%d reads matched the genome.' % (totalMatches, totalCount) #write result to STDOUT
"""