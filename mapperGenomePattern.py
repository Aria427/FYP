#!/usr/bin/env python
#This file contains the map step for the MapReduce implementation of the genome word count.       

import sys

#This function reads the reference genome as input to the mapper.
def readInput(file):
    for line in file:
        if line and line[0] != '>':
            yield line.upper().replace('\n', '').replace('N', '')

#This function returns a list of k-mers, similar to the sliding window approach.
def kmerList(sequence, k):
    kmers = []
    for i in xrange(0, len(sequence) + 1 - k):
        kmers.append( sequence[i:i+k] )
    return kmers            
            
def main():
    #input file comes from STDIN (standard input)
    subseqs = readInput(sys.stdin)
    
    last3 = '' #store the last 3 (4-1) bases of each line
    for s in subseqs:
        s = last3 + s #append to start of next line to handle patterns found between lines
        kmers = kmerList(s, 4)  #generate 4-mers of line
        last3 = s[-3:]
    
        #write results to STDOUT (standard output)
        for k in kmers:
            print '%s\t%s' % (k, 1) #tab-delimited
        #The output here will be the input for the reduce step.

if __name__ == '__main__':
    main()
    
#Run MapReduce job on Hadoop using:
    #bin/hadoop jar share/hadoop/tools/lib/hadoop-streaming-2.7.3.jar -file /home/hduser/mapperGenomePattern.py -mapper /home/hduser/mapperGenomePattern.py -file /home/hduser/reducerGenomePattern.py -reducer /home/hduser/reducerGenomePattern.py -input /user/hduser/HumanGenome.fa.gz  -output /user/hduser/4mercount-output
