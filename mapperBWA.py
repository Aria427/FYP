#!/usr/bin/env python
#This file contains the map step for the aligment MapReduce implementation.

import sys

def readInputGenome(file):
    for line in file:
        if line and line[0] != '>':
            yield line.rstrip().upper().replace('N', '').replace(' ','') 
            
def readInputReads(file):
    flag, sequence, quality = '', '', ''
    while True: #runs until EOF
        line = file.readline() 
        if not line: #reached EOF
            break
        
        if line.startswith('#'): #read details
            line = file.readline()
            pass
        
        elif line.startswith('>'): #>flags reads scores
            line = file.readline()
            pass
        
        else:
            line = line.split()
            flag = line[0]
            sequence = line[1]
            quality = line[2]
            
            if sequence == 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN':
                line = file.readline()
                pass
            
            else:
                sequence = sequence.rstrip().upper().replace('N', '').replace(' ', '')
                yield sequence  
                    
def main():
    #input file comes from STDIN (standard input)
    genomeSeq = readInputGenome(sys.stdin)
    readSeq = readInputReads(sys.stdin)
    
    overlap = '' 
    genomeStartIndex = 0
    for g in genomeSeq:
        g = overlap + g #append to start of next line to handle patterns found between lines

        #write results to STDOUT (standard output)
        print '%s\t%s' % ('G', (g, genomeStartIndex)) #tab-delimited key:value
        #key = character identifying genome
        #value = tuple with sequence and its start position
        
        overlap = g[-60:] #overlap = length of read
        genomeStartIndex += 50 #index of each subsequence is incremented by length of genome line 
        #The output here will be the input for the reduce step.
      
    readStartIndex = 0
    for r in readSeq:
        #write results to STDOUT (standard output)
        print '%s\t%s' % ('R', (r, readStartIndex)) #tab-delimited key:value
        #key = character identifying read
        #value = tuple with sequence and its start position
        
        readStartIndex += 60 #index of each read is incremented by length of read

if __name__ == '__main__':
    main()       

