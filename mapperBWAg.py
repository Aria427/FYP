#!/usr/bin/env python
#This file contains the genome map step for the aligment MapReduce implementation.

import sys

#This function reads the reference genome as input to the mapper.
def readInputGenome(file):
    for line in file:
        if line and line[0] != '>':
            yield line.rstrip().upper().replace('N', '').replace(' ','') 
                  
def main():
    #input files comes from STDIN (standard input)
    genomeSeq = readInputGenome(sys.stdin) #reads genome file as well as simple text  
    overlap = '' 
    genomeStartIndex = 0
        
    for g in genomeSeq:
     	g = overlap + g #append to start of next line to handle patterns found between lines
    
       	#write results to STDOUT (standard output)
  	print '%s\t%s' % ('G', (g, genomeStartIndex)) #tab-delimited key:value
        #key = character identifying genome
        #value = tuple with sequence and its start position
        #The output here will be the input for the reduce step.
            
        overlap = g[-100:] #overlap = length of read
        genomeStartIndex += 50 #index of each subsequence is incremented by length of genome line 
    
if __name__ == '__main__':
    main()       

