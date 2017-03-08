#!/usr/bin/env python
#This file contains the main functionality of the program.

import fileCompressionAndParsing
#import lzwCompression
#import genomeCompressionComparison

#import genomePattern
#import intGenomePattern

#import alignment
#import alignmentMatch
#import alignmentInt
import alignmentHadoop
#import analyseAlignment

#import visualisation

from path import path
#import gzip
#import sys
#import pdb
#import time

#pdb.set_trace()
#genomeFile = path('Data\HumanGenome.fa.gz').abspath()
#binaryGenomeFile = path('Output Data\HumanGenomeLong.bin').abspath()
genomeFile = path('Data\PhixGenome.fa').abspath()
#binaryGenomeFile = path('Output Data\PhixGenomeIntZip.bin').abspath() 
#fileCompressionAndParsing.parseGenomeInt(genomeFile, binaryGenomeFile)
#lzwCompression.compress(genomeFile, 'Output Data\HumanGenomeLZW.txt')
#genomeCompressionComparison.compressionComparison(path('Output Analysis Results\GenomeCompressionAnalysis.png').abspath()) 
    
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

#alignment.alignUncompressed(readsFile, genomeFile)
#alignment.alignCompressed(readsFile, binaryGenomeFile)


genome = fileCompressionAndParsing.parseGenomeString(genomeFile)

offsets, completeRQDict = [], {}
for read, quality in reads:
    offset, rqDict = alignmentHadoop.alignHamming(read, quality, genome)
    offsets.append(offset)
    completeRQDict.update(rqDict)
offsets = [o for oset in offsets for o in oset] #flatten list
print offsets
print completeRQDict


#textFile = path('Output Test Files\DataVisualisationTest.txt').abspath()
#textFile = visualisation.visualisationText(genome, offsets, textFile)
#jpgFile = path('Output Test Files\DataVisualisationTest.jpg').abspath()
#jpgFile = visualisation.visualisationJPG(genome, offsets, jpgFile) 
#pngFile = path('Output Test Files\DataVisualisationTest.png').abspath()
#pngFile = visualisation.visualisationGD(genome, offsets, pngFile)
#visualisation.visualisationTkinter(decodedGenome, offsets)