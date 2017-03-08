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
  
    #possible offsets for PhiX: 0-77
    #possible offsets for Human: 0-64,185,939
    reducers = 5
    for o1 in data: #o1[0]=offset, o1[1]=1
        if o1[0] <= 15:
            return 0
         
        elif o1[0] > 15 and o1[0] <= 30:
            return 1 % reducers
         
        elif o1[0] > 30 and o1[0] <= 45:
            return 2 % reducers
            
        elif o1[0] > 45 and o1[0] <= 60:
            return 3 % reducers
            
        else:
            return 4 % reducers
         
    #write results to STDOUT (standard output)
    #print '%s\t%s' % (offset, ones) #tab-delimited, key:offset of match with reads, value:list of 1 counts 
    #The output here will be the input for the reduce step.

if __name__ == '__main__':
    main()

#Run MapReduce job on Hadoop using:
#   bin/hadoop jar share/hadoop/tools/lib/hadoop-streaming-2.7.3.jar -file /home/hduser/mapperAligner.py -mapper /home/hduser/mapperAligner.py -file /home/hduser/partitionerAligner.py -partitioner /home/hduser/partitionerAligner.py -file /home/hduser/reducerAligner.py -reducer /home/hduser/reducerAligner.py -file /home/hduser/alignmentHadoop.py -file /home/hduser/matchingDistances.py -file /home/hduser/matchingBoyerMoore.py -file /home/hduser/matchingKmerIndex.py -file /home/hduser/matchingFmIndex.py -file /home/hduser/matchingSmithWaterman.py -file /home/hduser/matchingBurrowsWheeler.py -input /user/hduser/PhiXSequencingReads1000.fastq -output /user/hduser/bwaPhiX-output
#=> unkown class: partitioner    