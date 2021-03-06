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

import analyseMapReduceOutput
import visualisation

from path import path
import gzip
import sys
#import pdb
#import time

#pdb.set_trace()
#genomeFile = path('Data\HumanGenome.fa.gz').abspath() #64,185,939, line=50
#binaryGenomeFile = path('Output Data\HumanGenomeLong.bin').abspath()
genomeFile = path('Data\PhixGenome.fa').abspath() #77, line=70
binaryGenomeFile = path('Output Data\PhixGenomeIntZip.bin').abspath() 
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

"""
genome = fileCompressionAndParsing.parseGenomeString(genomeFile)
reads = fileCompressionAndParsing.parseReadsPhiXString(readsFile)

offsets, completeRQDict = [], {}
for read, quality in reads:
    offset, rqDict = alignmentHadoop.alignHamming(read, quality, genome)
    offsets.append(offset)
    completeRQDict.update(rqDict)
offsets = [o for oset in offsets for o in oset] #flatten list
print len(offsets)
#print completeRQDict                 

def readGenomeLines(file, lines=100):    
    with open(file, 'r') as f:
        firstLine = f.readline()
        if firstLine.startswith('>'):
            firstLine = ''
            pass #ignore header information
        while True:
            data = firstLine + ''.join(f.next().rstrip().upper().replace('N', '').replace(' ', '') for x in xrange(lines))
            if not data:
                break
            yield data
            firstLine = '' #first line in file is only needed for first iteration if not header
    
reads = fileCompressionAndParsing.parseReadsPhiXString(readsFile)
offsets, completeRQDict = [], {}
for read, quality in reads:
    genome = readGenomeLines(genomeFile)
    overlap = ''
    for g in genome:
        g = overlap + g
        offset, rqDict = alignmentHadoop.alignHamming(read, quality, g)
        offsets.append(offset)
        completeRQDict.update(rqDict)
        overlap = g[-100:] #100 for PhiX, 60 for Human
offsets = [o for oset in offsets for o in oset] #flatten list
print len(sorted(offsets, key=int))
#print completeRQDict    
"""
 
resultFiles = path('Output Test Files\PhiXHamming-output\part-000**').abspath() #* for multiple files    
analyseMapReduceOutput.analyseMR(resultFiles)
visualisation.visualise(binaryGenomeFile, resultFiles)
