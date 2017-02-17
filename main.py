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

out = open('Output Test Files\PhixLz77IntWholeFile.txt', 'w')
with open(path('Data\PhixGenome.fa').abspath()) as f:
    data = ''
    for line in f:
        if line and line[0] != '>':
            data = data + line.upper().replace('\n', '').replace('N', '')
    encoding = lz77.encode(data, 4)
    out.write(str(encoding)) 
out.close()

out = open('Output Test Files\PhixLz77IntPartsFile.txt', 'w')
with open(path('Data\PhixGenome.fa').abspath()) as f:
    #read genome into line-by-line generator     
    subseqs = (line.upper().replace('\n', '').replace('N', '') 
                for line in f if line and line[0] != '>') #ignore header line
            
    encoding = ''
    lineCount = 0
    for s in subseqs:
        try:
            #encoding every two lines is more accurate than every one line
            #window will lose some characters if every one line is considered
            s = s + next(subseqs) 
        except StopIteration:
            pass
        if lineCount == 0: #normal encoding
            encoding = encoding + str(lz77.encode(s, 4)).strip('[]') + ', '
        else: 
            #ignore first codeword as it produces non-sensical output;
            #this is due to slight reading difficulties as every 2 lines are taken
            encoding = encoding + str(lz77.encode(s, 4)[1:]).strip('[]') + ', '
        lineCount += 1
    out.write(encoding)
out.close()


#genomeCompressionComparison.compressionComparison(path('Output Test Files\GenomeCompressionAnalysis.png').abspath()) 
   

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
