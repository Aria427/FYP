#!/usr/bin/env python
#This file contains the reduce step for the MapReduce implementation of the genome word count.       

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
  
    for currentWord, group in groupby(data, itemgetter(0)):
        try:
            totalCount = sum(int(count) for word, count in group)
            print "%s\t%s" % (currentWord, totalCount) #write result to STDOUT
        except ValueError:
            pass #count = NAN

if __name__ == '__main__':
    main()
