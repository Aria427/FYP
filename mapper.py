#!/usr/bin/env python
#This file contains a basic map step for the MapReduce implementation of word count.       

import sys

def readInput(file):
    for line in file:
        yield line.split()

def main():
    #input file comes from STDIN (standard input)
    data = readInput(sys.stdin)
    for words in data:
        #write results to STDOUT (standard output)
        for word in words:
            print '%s\t%d' % (word, 1) #tab-delimited; trivial word count = 1
    #The output here will be the input for the reduce step.

if __name__ == '__main__':
    main()        