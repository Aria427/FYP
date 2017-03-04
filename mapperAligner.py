#!/usr/bin/env python
#This file contains the map step for the aligment MapReduce implementation.

import alignmentMatch
import sys

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
            
            #Each read has length = 60
            if sequence == 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN':
                line = file.readline()
                pass
            
            else:
                qualityNoNs = ''
                N = [pos for pos, char in enumerate(sequence) if char == 'N'] #positions of N in read
                sequenceNoNs = sequence.rstrip().upper().replace('N', '').replace(' ', '') 
                
                for i in range(len(quality)):
                    if i not in N: #remove indices corresponding to Ns in read
                        qualityNoNs = qualityNoNs + quality[i]
                
                yield sequenceNoNs, qualityNoNs

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
                N = [pos for pos, char in enumerate(sequence) if char == 'N'] #positions of N in read
                line = line.rstrip().upper().replace(' ', '')
                sequenceLines.append(line) #no whitespace in string sequence
                line = file.readline()
            sequence = ''.join(sequenceLines) #merge lines to form original sequence
            sequenceNoNs = sequence.replace('N', '') #remove Ns
            temp = sequenceNoNs
        
        elif not quality:
            qualityLines = []
            qualityNoNs = ''
            
            while True: #collect base qualities
                line = line.rstrip().replace(' ', '')
                qualityLines.append(line) 
                quality = ''.join(qualityLines) #merge lines to form quality
                if len(quality) >= len(sequence): #bases and qualities line up
                    break
                else:
                    line = file.readline()
            
            for i in range(len(quality)): 
                if i not in N: #remove indices corresponding to Ns in read
                    qualityNoNs = qualityNoNs + quality[i]
                
            yield temp, qualityNoNs
                  
def main():
    #hard-coded reference genome
    genomeFile = '/home/hduser/PhiXGenome.fa'
    with open(genomeFile, 'r') as f:
        genomeSeq = (line.rstrip().upper().replace('N', '').replace(' ','')
                    for line in f if line and line[0] != '>') #ignore header line

        totalMatches, totalCount, totalOffsets, completeRQDict = 0, 0, [], {}
        overlap = '' 
        genomeIndex = 0
       
        for g in genomeSeq:
            g = overlap + g #overlap is appended to the start of the next chunk 
  
            readSeq = readInputPhiXReads(sys.stdin)
            matchesCount, count, offsets, rqDict = alignmentMatch.alignFM(readSeq, g)
    	
            totalMatches += matchesCount
            totalCount = count #number of reads stays the same as every chunk goes through each read again
            totalOffsets.append(offsets)
            completeRQDict.update(rqDict)
            
            offsets = [o for offset in offsets for o in offset] #flatten list
            for i in range(len(offsets)):
                offsets[i] += genomeIndex #as each genome line has its own offsets
            
            #write results to STDOUT (standard output)
            print '%s\t%s' % (g, offsets) #tab-delimited, key:genome line, value:list of offsets of match with read
            #The output here will be the input for the reduce step.
            
            overlap = g[-100:] #100 for PhiX & 60 for Human  
            genomeIndex += 70 #index of each subsequence is incremented by length of genome line, 70 for PhiX & 50 for Human  
    
if __name__ == '__main__':
    main()       

#Run MapReduce job on Hadoop using:
#   bin/hadoop jar share/hadoop/tools/lib/hadoop-streaming-2.7.3.jar -file /home/hduser/mapperAligner.py -mapper /home/hduser/mapperAligner.py -file /home/hduser/reducerAligner.py -reducer /home/hduser/reducerAligner.py -file /home/hduser/alignmentMatch.py -file /home/hduser/matchingDistances.py -file /home/hduser/matchingBoyerMoore.py -file /home/hduser/matchingKmerIndex.py -file /home/hduser/matchingFmIndex.py -file /home/hduser/matchingSmithWaterman.py -file /home/hduser/matchingBurrowsWheeler.py -input /user/hduser/PhiXSequencingReads1000.fastq -output /user/hduser/bwaPhiX-output
