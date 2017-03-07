#!/usr/bin/env python
#This file contains the partition step for the aligment MapReduce implementation.

import sys

#This function reads in the output from the mapper using a generator.
def readMapperOutput(file, separator='\t'):
    for line in file:
        yield line.rstrip().split(separator, 1)
        
def main():
    #input comes from STDIN (standard input)
    data = readMapperOutput(sys.stdin) 
  
    tuples = []
    for o1 in data:
        tuples.append((o1[0], o1[1])) #(offset, 1)
    
    offsetOnes = {} #dictionary for partitions - 0: [1,1,1]
    for t in tuples:
        try:
            offsetOnes[t[0]].append(t[1]) #offset: list of 1s
        except KeyError:
            offsetOnes[t[0]] = [t[1]]

    #write results to STDOUT (standard output)
    for offset, ones in offsetOnes.items():
        print '%s\t%s' % (offset, ones) #tab-delimited, key:offset of match with reads, value:list of 1 counts 
    #The output here will be the input for the reduce step.

if __name__ == '__main__':
    main()

#Run MapReduce job on Hadoop using:
#   bin/hadoop jar share/hadoop/tools/lib/hadoop-streaming-2.7.3.jar -file /home/hduser/mapperAligner.py -mapper /home/hduser/mapperAligner.py -file /home/hduser/partitionerAligner.py -partitioner /home/hduser/partitionerAligner.py -file /home/hduser/reducerAligner.py -reducer /home/hduser/reducerAligner.py -file /home/hduser/alignmentMatch.py -file /home/hduser/matchingDistances.py -file /home/hduser/matchingBoyerMoore.py -file /home/hduser/matchingKmerIndex.py -file /home/hduser/matchingFmIndex.py -file /home/hduser/matchingSmithWaterman.py -file /home/hduser/matchingBurrowsWheeler.py -input /user/hduser/PhiXSequencingReads1000.fastq -output /user/hduser/bwaPhiX-output
#=> unkown class: partitioner    