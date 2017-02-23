#!/usr/bin/env python
#This file contains the reduce step for the MapReduce implementation.       

import sys

#This function reads in the output from the mapper using a generator.
def readMapperOutput(file, separator='\t'):
    for line in file:
        yield line.rstrip().split(separator)

def main():
    #input comes from STDIN (standard input)
    data = readMapperOutput(sys.stdin) 
  
    for isGR, tuple in data:
        print "%s\t%s" % (isGR, tuple) #write result to STDOUT

if __name__ == '__main__':
    main()

