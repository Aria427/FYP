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
        
    for currentOffset, group in groupby(data, itemgetter(0)):
        try:
            totalCount = sum(int(count) for offset, count in group)
            #write result to STDOUT
            print "%s\t%s" % (currentOffset, totalCount) #tab-delimited, key:offset of match with reads, value:count of match
        except ValueError:
            pass #count = NAN
    
if __name__ == '__main__':
    main()
