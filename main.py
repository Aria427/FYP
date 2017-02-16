#!/usr/bin/env python
#This file contains the main functionality of the program.

import fileParsing
import genomeCompressionComparison

import intGenomePattern
import genomePattern

import lz77
import lz78
import lzw

import alignmentInt
#import analyseAlignment
#import visualisation

from path import path
import numpy as np
import gzip

import sys
import pdb
import time

#pdb.set_trace()
genomeFile = path('Data\HumanGenome.fa.gz').abspath()
#binaryGenomeFile = path('Output Data\HumanGenomeLong.bin').abspath()
#fileParsing.parseGenomeLong(genomeFile, binaryGenomeFile)
#binaryGenomeFile = path('Output Data\PhixGenomeLongZip.bin').abspath() 
#fileParsing.parseGenomeLong(path('Data\PhixGenome.fa').abspath(), binaryGenomeFile)

#with open(binaryGenomeFile , 'rb') as f:
    #decodedGenome = np.fromfile(f, dtype=np.int)
    #print decodedGenome #length = 161856144
#d = reduce(lambda x,y: x+str(y), decodedGenome, '')
#decodedGenome = int(d)

#count = genomePattern.countIntWords(genomeFile)    
#with open('Output Test Files\intWordsCount.txt', 'w') as out:
#    out.write(str(count))
   
genomeCompressionComparison.compressionComparison(path('Output Test Files\GenomeCompressionAnalysis.png').abspath()) 
   

reads = fileParsing.parseReadsInt(path('Data\HumanSequencingReads.tsv.bz2').abspath()) 
#reads = fileParsing.parseReadsPhiXInt(path('Data\PhiXSequencingReads1000.fastq').abspath())


#analyseAlignment.plotTimeVsMatches(reads, decodedGenome, path('Output Test Files\AlignmentAnalysis.png').abspath())


#matchesCount, totalCount, offsets = alignmentInt.alignHamming(reads, decodedGenome)
#print '%d/%d reads matched the genome.' % (matchesCount, totalCount) #The result is not 100% but this is to be expected due to sequencing errors. 
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
