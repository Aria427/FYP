#!/usr/bin/env python
#This file contains the partition step for the aligment MapReduce implementation.

import sys

#This function reads in the output from the mapper using a generator.
def readMapperOutput(file, separator='\t'):
    for line in file:
        yield line.rstrip().split(separator, 1)
        
def main():
    #input comes from STDIN (standard input)
    data = readMapperOutput(sys.stdin) 
  
    tuples = []
    for o1 in data:
        tuples.append((o1[0], o1[1])) #(offset, 1)
    
    offsetOnes = {} #dictionary for partitions - 0: [1,1,1]
    for t in tuples:
        try:
            offsetOnes[t[0]].append(t[1]) #offset: list of 1s
        except KeyError:
            offsetOnes[t[0]] = [t[1]]

    #write results to STDOUT (standard output)
    for offset, ones in offsetOnes.items():
        print '%s\t%s' % (offset, ones) #tab-delimited, key:offset of match with reads, value:list of 1 counts 
    #The output here will be the input for the reduce step.

if __name__ == '__main__':
    main()
