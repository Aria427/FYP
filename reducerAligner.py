#!/usr/bin/env python
#This file contains the reduce step for the alignment MapReduce implementation.       

import sys
from itertools import groupby
from operator import itemgetter

#This function reads in the output from the mapper using a generator.
def readMapperOutput(file, separator='\t'):
    for line in file:
        yield line.rstrip().split(separator, 4)
        
def main():
    #input comes from STDIN (standard input)
    data = readMapperOutput(sys.stdin) 

    totalCount, genomeSubSeq, readQuality = 0, '', {}
    for currentOffset, group in groupby(data, itemgetter(0)):
        try:
            for offset, count, genomeSeq, read, quality in group:	
                totalCount += int(count)
                genomeSubSeq = genomeSeq #not concatenated as only 1 subsequence for an offset
                readQuality[read] = quality
            #write result to STDOUT
            #tab-delimited, key:offset of match with reads, value:<count of match, genome subsequence matched, reads:qualities matched>
            #offset corresponds to start position of genome subsequence
            print "%s\t%s\t%s\t%s" % (currentOffset, totalCount, genomeSubSeq, readQuality) 
        except ValueError:
            pass #count = NAN

if __name__ == '__main__':
    main()
