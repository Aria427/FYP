#!/usr/bin/env python
#This file contains the reduce step for the alignment MapReduce implementation.       

import sys
from itertools import groupby
from operator import itemgetter

#This function reads in the output from the mapper using a generator.
def readMapperOutput(file, separator='\t'):
    for line in file:
        yield line.rstrip().split(separator, 1)
        
def main():
    #input comes from STDIN (standard input)
    data = readMapperOutput(sys.stdin) 
  
    #for offsets, count in groupby(data, itemgetter(0)):
        #write results to STDOUT (standard output)
    #    print '%s\t%s' % (offsets, count)  #tab-delimited, key:list of offsets of match with read, value:number of matches
    
    for currentOffset, group in groupby(data, itemgetter(0)):
        try:
            totalCount = sum(int(count) for offset, count in group)
            print '%s\t%s' % (currentOffset, totalCount) #write result to STDOUT
        except ValueError:
            pass #count = NAN

if __name__ == '__main__':
    main()
