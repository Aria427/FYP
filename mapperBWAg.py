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

#Run MapReduce job on Hadoop using:
#   bin/hadoop jar share/hadoop/tools/lib/hadoop-streaming-2.7.3.jar -D mapreduce.job.maps=2 -file /home/hduser/mapperBWAg.py -mapper /home/hduser/mapperBWAg.py -input /user/hduser/PhiXGenome.fa -file /home/hduser/mapperBWAr.py -mapper /home/hduser/mapperBWAr.py -input /user/hduser/PhiXSequencingReads1000.fastq  -file /home/hduser/reducerBWA.py -reducer /home/hduser/reducerBWA.py -file /home/hduser/alignmentString.py -file /home/hduser/matchingString.py  -output /user/hduser/bwaPhiX-output
