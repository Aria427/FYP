#!/usr/bin/env python
#This file contains the main functionality of the program.

import fileParsing
import alignmentInt
#import analyseAlignment
#import visualisation
import intGenomePatternTesting
import genomePatternTesting

from path import path
import numpy as np

import sys
import pdb
import time

#pdb.set_trace()
genomeFile = path('Data\TenMB.fa').abspath()
#binaryGenomeFile = path('Output Data\HumanGenomeZip.bin').abspath()
#fileParsing.parseGenomeInt(genomeFile, binaryGenomeFile)
#binaryGenomeFile = path('Output Data\PhixGenome.bin').abspath() 
#fileParsing.parseGenomeInt(path('Data\PhixGenome.fa').abspath(), binaryGenomeFile)

#with open(binaryGenomeFile , 'rb') as f:
    #decodedGenome = np.fromfile(f, dtype=np.int)
    #print decodedGenome #length = 161856144
#d = reduce(lambda x,y: x+str(y), decodedGenome, '')
#decodedGenome = int(d)

#integers, integersCount = intGenomePatternTesting.countIntegers(decodedGenome)
#pairsCount = intGenomePatternTesting.countIntegerPairs(decodedGenome) 
#histogram, bins = intGenomePatternTesting.createHist(decodedGenome)
#intGenomePatternTesting.createIntegerHistogram(integersCount)
#intGenomePatternTesting.createPairHistogram(pairsCount)

genomePatternTesting.countIntWordsChunks(genomeFile)
genomePatternTesting.countIntWordsBytes(genomeFile)
genomePatternTesting.countIntWords(genomeFile)


reads = fileParsing.parseReadsInt(path('Data\HumanSequencingReads.tsv.bz2').abspath()) 
#reads = fileParsing.parseReadsPhiXInt(path('Data\PhiXSequencingReads1000.fastq').abspath())


#analyseAlignment.plotTimeVsMatches(reads, decodedGenome, path('Output Test Files\AlignmentAnalysis.png').abspath())

#start = time.time()
#matchesCount, totalCount, offsets = alignmentInt.alignHamming(reads, decodedGenome)
#print '%d/%d reads matched the genome.' % (matchesCount, totalCount) #The result is not 100% but this is to be expected due to sequencing errors. 
#print offsets
#end = time.time()
#print end-start

"""
textFile = path('Output Test Files\DataVisualisationTest.txt').abspath()
textFile = visualisation.visualisationText(genome, offsets, textFile)

jpgFile = path('Output Test Files\DataVisualisationTest.jpg').abspath()
jpgFile = visualisation.visualisationJPG(genome, offsets, jpgFile) 

pngFile = path('Output Test Files\DataVisualisationTest.png').abspath()
pngFile = visualisation.visualisationGD(genome, offsets, pngFile)
"""
#visualisation.visualisationTkinter(decodedGenome, offsets)
