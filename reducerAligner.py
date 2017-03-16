#!/usr/bin/env python
#This file contains the reduce step for the alignment MapReduce implementation.       

import sys
from itertools import groupby
from operator import itemgetter

#This function reads in the output from the partitioner using a generator.
def readPartitionerOutput(file, separator='\t'):
    for line in file:
        yield line.rstrip().split(separator, 4)
        
def main():
    #input comes from STDIN (standard input)
    data = readPartitionerOutput(sys.stdin) 

    totalCount, readQuality = 0, {}
    for currentOffset, group in groupby(data, itemgetter(0)):
        try:
            for offset, count, read, quality in group:	
                totalCount += int(count)
                readQuality[read] = quality
            #write result to STDOUT
            #tab-delimited, key:offset of match with reads, value:<count of match, reads:qualities matched>
            print "%s\t%s\t%s" % (currentOffset, totalCount, readQuality) 
        except ValueError:
            pass #count = NAN

if __name__ == '__main__':
    main()
