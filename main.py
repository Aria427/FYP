#!/usr/bin/env python
#This file contains the main functionality of the program.

import fileCompressionAndParsing
import lzwCompression
import genomeCompressionComparison

import genomePattern
import intGenomePattern

import alignment
import alignmentMatch
import alignmentInt
#import analyseAlignment

#import visualisation

from path import path
import numpy as np
import gzip
import sys
import pdb
import time
import struct

#pdb.set_trace()
#genomeFile = path('Data\HumanGenome.fa.gz').abspath()
#binaryGenomeFile = path('Output Data\HumanGenomeLong.bin').abspath()
genomeFile = path('Data\PhixGenome.fa').abspath()
binaryGenomeFile = path('Output Data\PhixGenomeIntZip.bin').abspath() 
#fileCompressionAndParsing.parseGenomeInt(genomeFile, binaryGenomeFile)
#lzwCompression.compress(genomeFile, 'Output Data\HumanGenomeLZW.txt')
#genomeCompressionComparison.compressionComparison(path('Output Analysis Results\GenomeCompressionAnalysis.png').abspath()) 

#with open(binaryGenomeFile , 'rb') as f:
    #decodedGenome = np.fromfile(f, dtype=np.int)
    #print decodedGenome #length = 161856144
#d = reduce(lambda x,y: x+str(y), decodedGenome, '')
#decodedGenome = int(d)
    
#count = genomePattern.countIntWords(genomeFile)    
#with open('Output Test Files\intWordsCount.txt', 'w') as out:
#    out.write(str(count))

#readsFile = path('Data\HumanSequencingReads.tsv.bz2').abspath()
#reads = fileCompressionAndParsing.parseReadsString(readsFile)
#reads = fileCompressionAndParsing.parseReadsInt(readsFile)  
readsFile = path('Data\PhiXSequencingReads1000.fastq').abspath()
reads = fileCompressionAndParsing.parseReadsPhiXString(readsFile)
#reads = fileCompressionAndParsing.parseReadsPhiXInt(readsFile)
#print len(next(reads)) #length of PhiX read = 100; length of human read = 58
#There are 28,094,847 human reads in total.

with open(genomeFile, 'r') as f:
    genomeSeq = (line.rstrip().upper().replace('N', '').replace(' ','')
                    for line in f if line and line[0] != '>') #ignore header line

    totalMatches, totalCount, totalOffsets, completeRQDict = 0, 0, [], {}
    overlap = '' 
    genomeIndex = 0
       
    for g in genomeSeq:
        g = overlap + g #overlap is appended to the start of the next chunk 
  
        readSeq = fileCompressionAndParsing.parseReadsPhiXString(readsFile)
        matchesCount, count, offsets, rqDict = alignmentMatch.alignFM(readSeq, g)
    
        offsets = [o for offset in offsets for o in offset] #flatten list
        for i in range(len(offsets)):
            offsets[i] += genomeIndex #as each genome line has its own offsets
                
        totalMatches += matchesCount
        totalCount = count #number of reads stays the same as every chunk goes through each read again
        totalOffsets.append(offsets)
        completeRQDict.update(rqDict)
            
        overlap = g[-100:] #100 for PhiX & 60 for Human  
        #genomeIndex += 70 #index of each subsequence is incremented by length of genome line, 70 for PhiX & 50 for Human  

    print '%d/%d reads matched the genome.' % (totalMatches, totalCount) #The result is not 100% but this is to be expected due to sequencing errors. 
    totalOffsets = [o for offset in totalOffsets for o in offset] #flatten list
    print len(sorted(totalOffsets, key=int))
    print len(completeRQDict)

genome = fileCompressionAndParsing.parseGenomeString(genomeFile)
reads = fileCompressionAndParsing.parseReadsPhiXString(readsFile)
#for r, q in reads:
#    print r, ':', q

matchesCount, totalCount, offsets, rqDict = alignmentMatch.alignHamming(reads, genome)
print '%d/%d reads matched the genome.' % (matchesCount, totalCount) #The result is not 100% but this is to be expected due to sequencing errors. 
offsets = [o for offset in offsets for o in offset]
print offsets #len(offsets) = matchesCount
print len(rqDict) #=> number unique reads which matches

#=> aligning line by line causes loss of data


#alignment.alignUncompressed(readsFile, genomeFile)
#alignment.alignCompressed(readsFile, binaryGenomeFile)

#with open(binaryGenomeFile , 'rb') as f:
    #decodedGenome = np.fromfile(f, dtype=np.int)
#matchesCount, totalCount, offsets = alignmentInt.alignHamming(reads, decodedGenome)
#print '%d/%d reads matched the genome.' % (matchesCount, totalCount)
#print offsets

"""
textFile = path('Output Test Files\DataVisualisationTest.txt').abspath()
textFile = visualisation.visualisationText(genome, offsets, textFile)

jpgFile = path('Output Test Files\DataVisualisationTest.jpg').abspath()
jpgFile = visualisation.visualisationJPG(genome, offsets, jpgFile) 

pngFile = path('Output Test Files\DataVisualisationTest.png').abspath()
pngFile = visualisation.visualisationGD(genome, offsets, pngFile)
"""
#visualisation.visualisationTkinter(decodedGenome, offsets)
