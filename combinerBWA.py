#!/usr/bin/env python
#This file contains the combine step for the aligment MapReduce implementation.

import sys

#This function reads in the output from the mapper using a generator.
def readMapperOutput(file, separator='\t'):
    for line in file:
        yield line.rstrip().split(separator, 1)
                    
def main():
    #input comes from STDIN (standard input)
    data = readMapperOutput(sys.stdin) 
  
    genomeValue = []
    readValue = []
    for isGR, tupleValue in data:
        if isGR == 'G':
            genomeValue.append(tupleValue)
        elif isGR == 'R':
            readValue.append(tupleValue)
        else:
            pass
        
    #write result to STDOUT    
    print "%s\t%s" % ('G', genomeValue) #tab-delimited key:value
    print "%s\t%s" % ('R', readValue)   #tab-delimited key:value
    #key = character identifying genome
    #value = list of tuples with sequence and its start position

if __name__ == '__main__':
    main()  
