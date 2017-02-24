#!/usr/bin/env python
#This file contains the map step for the aligment MapReduce implementation.

import sys

#This function reads the reference genome as input to the mapper.
def readInputGenome(file):
    for line in file:
        if line and line[0] != '>':
            yield line.rstrip().upper().replace('N', '').replace(' ','') 

#This function reads the sequencing reads as input to the mapper.        
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

def readInputPhiXReads(file):  
    readID, sequence, quality = '', '', ''
    while True: #runs until EOF
        line = file.readline() 
        if not line: #reached EOF
            break

        if line.startswith('@'): #first line of read/record
            #reset to default values
            readID = line.rstrip()
            sequence = ''
            quality = ''   

        elif not readID: #if no previous line starts with @
            readID = line.rstrip() #get first ID
            continue

        elif not sequence:
            sequenceLines = [] 
            while not line.startswith('+'): #not placeholder line (third line)
                #rstrip() - removes leading/trailing whitespace
                #replace() - removes whitespace from within string
                line = line.rstrip().upper().replace('N', '').replace(' ', '')
                sequenceLines.append(line) #no whitespace in string sequence
                line = file.readline()
            sequence = ''.join(sequenceLines) #merge lines to form sequence
            yield sequence
        
        elif not quality:
            quality = []
            while True: #collect base qualities
                quality += line.rstrip().replace(' ', '') 
                if len(quality) >= len(sequence): #bases and qualities line up
                    break
                else:
                    line = file.readline()
                  
def main():
    #input file comes from STDIN (standard input)
    if sys.stdin.name.endswith('fastq'):#'.bz2'):
        readSeq = readInputPhiXReads(sys.stdin)
        readStartIndex = 0
        
        for r in readSeq:
            #write results to STDOUT (standard output)
            print '%s\t%s' % ('R', (r, readStartIndex)) #tab-delimited key:value
            #key = character identifying read
            #value = tuple with sequence and its start position
            
            readStartIndex += 100 #index of each read is incremented by length of read
    
    elif sys.stdin.name.endswith('.bz2'):
        readSeq = readInputReads(sys.stdin)
        readStartIndex = 0
        
        for r in readSeq:
            #write results to STDOUT (standard output)
            print '%s\t%s' % ('R', (r, readStartIndex)) #tab-delimited key:value
            #key = character identifying read
            #value = tuple with sequence and its start position
            
            readStartIndex += 60 #index of each read is incremented by length of read
                   
    else:
        genomeSeq = readInputGenome(sys.stdin) #reads genome file as well as simple text  
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

if __name__ == '__main__':
    main()       

