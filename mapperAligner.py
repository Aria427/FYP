#!/usr/bin/env python
#This file contains the map step for the aligment MapReduce implementation.

import alignmentHadoop
import sys
import gzip

#This function reads the genome into a string.
def readGenome(file): 
    with open(file, 'r') as f: 
        for line in f:
            if line and line[0] != '>': #ignore header line with genome information
                l = line.rstrip().upper().replace('N', '').replace(' ', '') 
                yield l

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
        
        elif not sequence or not quality:
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
    #genomeFile = '/home/hduser/HumanGenome.fa.gz'            
    
    readSeq = readInputPhiXReads(sys.stdin)
    
    for read, quality in readSeq:
        genome = readGenome(genomeFile)
        overlap = ''
        
        for g in genome:
            g = overlap + g
            offset, rqDict = alignmentHadoop.alignHamming(read, quality, g)
            
            #write results to STDOUT (standard output)
            for o in offset: #to remove empty list and '[' ']' characters
                print '%s\t%s' % (o, 1) #tab-delimited, key:offset of match with reads, value:default count of 1 
            #The output here will be the input for the reduce step  
            
            overlap = g[-100:] #100 for PhiX, 60 for Human
    
    
if __name__ == '__main__':
    main()            

#Run MapReduce job on Hadoop using:
# PhiX
#   bin/hadoop jar share/hadoop/tools/lib/hadoop-streaming-2.7.3.jar -D mapred.reduce.tasks=4 -D mapred.output.key.comparator.class=org.apache.hadoop.mapred.lib.KeyFieldBasedComparator -D  mapred.text.key.comparator.options=-n -file /home/hduser/mapperAligner.py -mapper /home/hduser/mapperAligner.py -file /home/hduser/reducerAligner.py -reducer /home/hduser/reducerAligner.py -file /home/hduser/alignmentMatch.py -file /home/hduser/matchingDistances.py -file /home/hduser/matchingBoyerMoore.py -file /home/hduser/matchingKmerIndex.py -file /home/hduser/matchingFmIndex.py -file /home/hduser/matchingSmithWaterman.py -file /home/hduser/matchingBurrowsWheeler.py -input /user/hduser/PhiXSequencingReads1000.fastq -output /user/hduser/PhiXHamming-output -partitioner org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner -jobconf stream.map.output.field.separator=\t -jobconf stream.num.map.output.key.fields=1 -jobconf map.output.key.field.separator=\t -jobconf num.key.fields.for.partition=1
"""
    bin/hadoop jar share/hadoop/tools/lib/hadoop-streaming-2.7.3.jar \
    -D mapred.reduce.tasks=4 \
    -D mapred.output.key.comparator.class=org.apache.hadoop.mapred.lib.KeyFieldBasedComparator 
    -D  mapred.text.key.comparator.options=-n \ #sort numerically in shuffle & sort phase
    -file /home/hduser/mapperAligner.py -mapper /home/hduser/mapperAligner.py \
    -file /home/hduser/reducerAligner.py -reducer /home/hduser/reducerAligner.py \
    -file /home/hduser/alignmentHadoop.py \
    -file /home/hduser/matchingDistances.py -file /home/hduser/matchingBoyerMoore.py \
    -file /home/hduser/matchingKmerIndex.py -file /home/hduser/matchingFmIndex.py \
    -file /home/hduser/matchingSmithWaterman.py -file /home/hduser/matchingBurrowsWheeler.py \
    -input /user/hduser/PhiXSequencingReads1000.fastq \
    -output /user/hduser/PhiXHamming-output \
    -partitioner org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner \
    -jobconf stream.map.output.field.separator=\t \ #seperator between key and value
    -jobconf stream.num.map.output.key.fields=1 \ #no of \t before end of key
    -jobconf map.output.key.field.separator=\t \
    -jobconf num.key.fields.for.partition=1 #no of fields used for partitioning
"""
# Human
#   bin/hadoop jar share/hadoop/tools/lib/hadoop-streaming-2.7.3.jar -D mapred.reduce.tasks=5 -D mapred.output.key.comparator.class=org.apache.hadoop.mapred.lib.KeyFieldBasedComparator -D  mapred.text.key.comparator.options=-n -file /home/hduser/mapperAligner.py -mapper /home/hduser/mapperAligner.py -file /home/hduser/reducerAligner.py -reducer /home/hduser/reducerAligner.py -file /home/hduser/alignmentHadoop.py -file /home/hduser/matchingDistances.py -file /home/hduser/matchingBoyerMoore.py -file /home/hduser/matchingKmerIndex.py -file /home/hduser/matchingFmIndex.py -file /home/hduser/matchingSmithWaterman.py -file /home/hduser/matchingBurrowsWheeler.py -input /user/hduser/HumanSequencingReads.tsv.bz2 -output /user/hduser/HumanHamming-output -partitioner org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner -jobconf stream.map.output.field.separator=\t -jobconf stream.num.map.output.key.fields=1 -jobconf map.output.key.field.separator=\t -jobconf num.key.fields.for.partition=1
