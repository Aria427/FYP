#!/usr/bin/env python
#This file contains the map step for the aligment MapReduce implementation.

import alignmentString

import sys
from path import path

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
    #file = open(file, 'r')
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
    #file.close()
                  
def main():
    #hard-coded reference genome
    genomeFile = '/home/hduser/PhiXGenome.fa'
    with open(genomeFile, 'r') as f:
        genomeSeq = (line.rstrip().upper().replace('N', '').replace(' ','')
                    for line in f if line and line[0] != '>') #ignore header line

    	totalMatches, totalCount, totalOffsets = 0, 0, []
    	overlap = '' 
       
    	for g in genomeSeq:
         g = overlap + g #overlap is appended to the start of the next chunk 
  
         readSeq = readInputPhiXReads(sys.stdin)
         matchesCount, count, offsets = alignmentString.alignFM(readSeq, g)
	
         totalMatches += matchesCount
         totalCount = count #number of reads stays the same as every chunk goes through each read again
         totalOffsets.append(offsets)
         overlap = g[-100:] #100 for PhiX, 60 for Human  
        
         print '%d/%d reads matched the genome.' % (totalMatches, totalCount) 

if __name__ == '__main__':
    main()       

#Run MapReduce job on Hadoop using:
#   bin/hadoop jar share/hadoop/tools/lib/hadoop-streaming-2.7.3.jar -file /home/hduser/mapperAligner.py -mapper /home/hduser/mapperAligner.py -file /home/hduser/reducerAligner.py -reducer /home/hduser/reducerAligner.py -file /home/hduser/alignmentString.py -file /home/hduser/matchingString.py -input /user/hduser/PhiXSequencingReads1000.fastq -output /user/hduser/bwaPhiX-output

"""                
def main():
    #hard-coded reference genome
    genomeFile = path('Data\HumanGenome.fa.gz').abspath()
    genomeSeq = readInputGenome(genomeFile)
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
    
    #input file - sequencing reads - comes from STDIN (standard input)
    if sys.stdin.readline().startswith('#'): #if the first line of the input starts with '#' => human sequencing reads
        readSeq = readInputReads(sys.stdin)
        readStartIndex = 0
        
        for r in readSeq:
            #write results to STDOUT (standard output)
            print '%s\t%s' % ('R', (r, readStartIndex)) #tab-delimited key:value
            #key = character identifying read
            #value = tuple with sequence and its start position
            #The output here will be the input for the reduce step.
            
            readStartIndex += 60 #index of each read is incremented by length of read
      
    else: #if the first line of the input starts with '@' => PhiX sequencing reads
        readSeq = readInputPhiXReads(sys.stdin)
        readStartIndex = 0
        
        for r in readSeq:
            #write results to STDOUT (standard output)
            print '%s\t%s' % ('R', (r, readStartIndex)) #tab-delimited key:value
            #key = character identifying read
            #value = tuple with sequence and its start position
            #The output here will be the input for the reduce step.
            
            readStartIndex += 100 #index of each read is incremented by length of read 
"""   
