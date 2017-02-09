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
genomeFile = path('Data\HumanGenome.fa.gz').abspath()
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

start = time.time()
genomePatternTesting.countIntWordsChunks(genomeFile) 
#Counter({'TTTT': 298456206, 'AAAA': 295878805, 'AAAT': 233204521, 'TTTG': 221603116, 
#'CCCA': 220519406, 'AAAG': 213147885, 'CCCT': 212879264, 'TTTA': 196975835, 
#'GGGA': 182228315, 'TTTC': 181239017, 'GGGG': 158825594, 'CCCC': 157575620, 
#'GGGT': 153734106, 'AAAC': 153370668, 'GGGC': 129676917, 'CCCG': 30886900, 
#'CCCH': 433, 'HHHR': 430, 'IIIV': 325, 'KKKI': 325, 'VVVA': 291, 'AAAL': 259, 
#'LLLT': 259, 'RRRK': 212, 'RRRU': 125, 'UUUK': 112, 'GGGL': 88, 'LLLV': 88, 
#'RRRG': 79, 'VVVG': 45, 'AAAD': 42, 'VVVR': 42, 'OOOM': 42, 'DDDO': 42, 
#'RRRA': 42, 'VVVT': 26, 'MMMG': 25, 'VVVC': 21, 'UUUG': 9, 'MMMC': 7, 'JJJH': 6, 
#'HHHV': 6, 'MMMA': 6, 'RRRJ': 6, 'MMMT': 5, 'XXXK': 3, 'RRRX': 3, 'BBBV': 2, 
#'KKKB': 2, 'YYYK': 1, 'RRRY': 1, 'RRRM': 1}
end = time.time()
print end-start


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
