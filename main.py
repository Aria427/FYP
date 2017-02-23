#!/usr/bin/env python
#This file contains the main functionality of the program.

import fileCompressionAndParsing
import lzwCompression
import genomeCompressionComparison

import intGenomePattern
import genomePattern

import alignmentString
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
genomeFile = path('Data\HumanGenome.fa.gz').abspath()
#binaryGenomeFile = path('Output Data\HumanGenomeLong.bin').abspath()
#genomeFile = path('Data\PhixGenome.fa').abspath()
#binaryGenomeFile = path('Output Data\PhixGenomeIntZip.bin').abspath() 
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

#reads = fileCompressionAndParsing.parseReadsInt(path('Data\HumanSequencingReads.tsv.bz2').abspath()) 
reads = fileCompressionAndParsing.parseReadsString(path('Data\HumanSequencingReads.tsv.bz2').abspath()) 
#reads = fileCompressionAndParsing.parseReadsPhiXInt(path('Data\PhiXSequencingReads1000.fastq').abspath())
#reads = fileCompressionAndParsing.parseReadsPhiXString(path('Data\PhiXSequencingReads1000.fastq').abspath())
#print len(next(reads)) #length of PhiX read = 100; length of human read = 58
#There are 28,094,847 human reads in total.

#genome = fileCompressionAndParsing.parseGenomeString(genomeFile)
#matchesCount, totalCount, offsets = alignmentString.alignHamming(reads, genome)
#print '%d/%d reads matched the genome.' % (matchesCount, totalCount)
#print offsets

def readInChunks(genomeFile, chunkSize=65536):
   while True:
       data = f.read(chunkSize).rstrip().upper().replace('N', '').replace(' ', '')
       if not data:
           break
       yield data

def chunked(file, chunk_size):
    return iter(lambda: file.read(chunk_size).rstrip().upper().replace('N', '').replace(' ', ''), '')     
   
totalMatches, totalCount, totalOffsets = 0, 0, []
with gzip.open(genomeFile, 'r') as f:
    subseqs = chunked(f, 100000) #file size = 3,273,481,150 bytes => chunk = 100,000 bytes
    lastRead = '' #overlap
    chunkCount = 0
    for s in subseqs:
        print 'Chunk read'
        chunkCount += 1
        reads = fileCompressionAndParsing.parseReadsString(path('Data\HumanSequencingReads.tsv.bz2').abspath()) 
        #reads = fileCompressionAndParsing.parseReadsPhiXString(path('Data\PhiXSequencingReads1000.fastq').abspath())
        s = lastRead + s
        matchesCount, count, offsets = alignmentString.alignFM(reads, s)
        totalMatches += matchesCount
        totalCount = count
        totalOffsets.append(offsets)
        lastRead = s[-58:]
        print 'Chunk %d/32735 aligned' % chunkCount        
        print '%d/%d reads matched the genome.' % (totalMatches, totalCount)
#print totalOffsets
     

#matchesCount, totalCount, offsets = alignmentInt.alignHamming(reads, decodedGenome)
#print '%d/%d reads matched the genome.' % (matchesCount, totalCount) #The result is not 100% but this is to be expected due to sequencing errors. 
#print offsets

#analyseAlignment.plotTimeVsMatches(reads, decodedGenome, path('Output Test Files\AlignmentAnalysis.png').abspath())

"""
textFile = path('Output Test Files\DataVisualisationTest.txt').abspath()
textFile = visualisation.visualisationText(genome, offsets, textFile)

jpgFile = path('Output Test Files\DataVisualisationTest.jpg').abspath()
jpgFile = visualisation.visualisationJPG(genome, offsets, jpgFile) 

pngFile = path('Output Test Files\DataVisualisationTest.png').abspath()
pngFile = visualisation.visualisationGD(genome, offsets, pngFile)
"""
#visualisation.visualisationTkinter(decodedGenome, offsets)
