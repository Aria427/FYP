#!/usr/bin/env python

import sys

#This function reads in the output from the mapper using a generator.
def readMapperOutput(file, separator='\t'):
    for line in file:
        yield line.rstrip().split(separator, 1)

def main():
    #input comes from STDIN (standard input)
    data = readMapperOutput(sys.stdin, separator='\t') 
  
    for word, count in data:
        try:
            if count == '1':
                pass #ignore words which occur only once
            else:
                print "%s\t%s" % (word, count) #write result to STDOUT
        except ValueError:
            pass #when count not a number

if __name__ == '__main__':
    main()
