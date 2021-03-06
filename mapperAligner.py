#!/usr/bin/env python
#This file contains the map step for the aligment MapReduce implementation.

import alignmentHadoop
import sys
import smart_open #Add bootstrap action to s3://fyp-input/installPythonModules.sh

#This function reads the genome in chunks (groups of bytes) from the file stored in S3.
def readGenomeChunks(s3File, bytesNum=100):
    with gzip.open(s3File, 'r') as f:
        f.readline()
        for chunk in iter(lambda: f.read(bytesNum), ''):
            data = chunk.rstrip().upper().replace('N', '').replace('\n', '').replace(' ', '')
            yield data          

#This function reads the whole genome from the file stored in S3.
def readGenome(s3File):
    genome = ''
    with gzip.open(s3File, 'r') as f: 
        for line in f:
            if line and line[0] != '>': #ignore header line with genome information
                l = line.rstrip().upper().replace('N', '').replace(' ', '')
                genome += l
        return genome

#This function reads the genome in chunks (groups of bytes) from the file stored in S3 using smart_open.
def readGenomeEMR(s3FileUrl, bytesNum=100):
    with smart_open.smart_open(s3FileUrl, 'r') as dataFile:
        firstLine = dataFile.readline()
        if firstLine.startswith('>'):
            firstLine = ''
            pass #ignore header information
        data = firstLine + dataFile.read(bytesNum).rstrip().upper().replace('N', '').replace(' ', '')
        yield data
        
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
    #genomeFile = '/home/hduser/PhiXGenome.fa'
    #genomeFile = '/home/hduser/HumanGenome.fa.gz' 
    
    #For Amazon EMR:
    genomeFile = 'https://s3.eu-central-1.amazonaws.com/fyp-input/HumanGenome.fa' 
       
    readSeq = readInputReads(sys.stdin) #Human reads=28,094,847
    
    #readsMatched, readsMismatched = 0, 0
    for read, quality in readSeq: 
        #Human genome=64,185,939 lines -> 3,273,481,150 bytes
        genome = readGenomeEMR(genomeFile, 500000000) #500,000,000 bytes = 0.5GB
        overlap = ''
        
        for g in genome:
            g = overlap + g
            offset, rqDict = alignmentHadoop.alignHamming(read, quality, g)
            
            #write results to STDOUT (standard output)
            for o in offset: #to remove empty list and '[' ']' characters
                #if o == '':
                #    readsMismatched += 1
                #else:
                #tab-delimited, key:offset of match with reads, value:<default count of 1, genome subsequence matched, read matched, corresponding quality> 
                print '%s\t%s\t%s\t%s\t%s' % (o, 1, g[o:o+len(read)], read, quality)  
                    #readsMatched += 1
                #The output here will be the input for the reduce step  
            
            overlap = g[-100:] #100 for PhiX, 60 for Human   
 
if __name__ == '__main__':
    main()            

#Run MapReduce job on Hadoop using:
# PhiX
#   bin/hadoop jar share/hadoop/tools/lib/hadoop-streaming-2.7.3.jar -D stream.map.output.field.separator=\t -D stream.num.map.output.key.fields=1 -D mapreduce.map.output.key.field.separator=\t -D mapreduce.partition.keypartitioner.options=-k1,1 -D  mapreduce.job.reduces=2 -D  mapreduce.job.output.key.comparator.class=org.apache.hadoop.mapred.lib.KeyFieldBasedComparator -D  mapreduce.partition.keycomparator.options=-n -file /home/hduser/mapperAligner.py -mapper /home/hduser/mapperAligner.py -file /home/hduser/reducerAligner.py -reducer /home/hduser/reducerAligner.py -file /home/hduser/alignmentMatch.py -file /home/hduser/matchingDistances.py -file /home/hduser/matchingBoyerMoore.py -file /home/hduser/matchingKmerIndex.py -file /home/hduser/matchingFmIndex.py -file /home/hduser/matchingSmithWaterman.py -file /home/hduser/matchingBurrowsWheeler.py -input /user/hduser/PhiXSequencingReads1000.fastq -output /user/hduser/PhiXHamming-output -partitioner org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner 

# Human
#   bin/hadoop jar share/hadoop/tools/lib/hadoop-streaming-2.7.3.jar -D stream.map.output.field.separator=\t -D stream.num.map.output.key.fields=1 -D mapreduce.map.output.key.field.separator=\t -D mapreduce.partition.keypartitioner.options=-k1,1 -D mapreduce.job.reduces=10 -D mapreduce.job.output.key.comparator.class=org.apache.hadoop.mapred.lib.KeyFieldBasedComparator -D mapreduce.partition.keycomparator.options=-n -file /home/hduser/mapperAligner.py -mapper /home/hduser/mapperAligner.py -file /home/hduser/reducerAligner.py -reducer /home/hduser/reducerAligner.py -file /home/hduser/alignmentHadoop.py -file /home/hduser/matchingDistances.py -file /home/hduser/matchingBoyerMoore.py -file /home/hduser/matchingKmerIndex.py -file /home/hduser/matchingFmIndex.py -file /home/hduser/matchingSmithWaterman.py -file /home/hduser/matchingBurrowsWheeler.py -input /user/hduser/HumanSequencingReads.tsv -output /user/hduser/HumanHamming-output -partitioner org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner
