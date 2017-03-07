#!/usr/bin/env python
#This file contains the reduce step for the alignment MapReduce implementation.       

import sys
#from itertools import groupby
#from operator import itemgetter

#This function reads in the output from the partitioner using a generator.
def readPartitionerOutput(file, separator='\t'):
    for line in file:
        yield line.rstrip().split(separator, 1)
        
def main():
    #input comes from STDIN (standard input)
    data = readPartitionerOutput(sys.stdin) 
  
    dictionary = {}
    for o1s in data:
        dictionary[o1s].append(o1s[1]) #offset: list of 1s
        
    offsetCount = {}
    for offset, ones in dictionary.items():
        total = sum(ones)
        offsetCount[offset] = total

    #write results to STDOUT (standard output)
    for offset, count in offsetCount.items():
        print '%s\t%s' % (offset, count) #tab-delimited, key:offset of match with reads, value:count of match
    #The output here will be the input for the reduce step.
        
    """
    for currentOffset, group in groupby(data, itemgetter(0)):
        try:
            totalCount = sum(int(count) for offset, count in group)
            #write result to STDOUT
            print "%s\t%s" % (currentOffset, totalCount) #tab-delimited, key:offset of match with reads, value:count of match
        except ValueError:
            pass #count = NAN
    """

if __name__ == '__main__':
    main()
