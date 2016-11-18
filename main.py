#!/usr/bin/env python
#This file contains the main functionality of the program.

import fileParsing
import matching
import alignment
import alignmentInt
import analyseAlignment
import visualisation

from path import path
from bitarray import bitarray
import struct
import numpy 
import sys
import pdb

#pdb.set_trace()

#binaryGenomeFile = path('Output Data\HumanGenomeZip.bin').abspath()
#fileParsing.parseGenomeInt(path('Data\HumanGenome.fa.gz').abspath(), binaryGenomeFile)
binaryGenomeFile = path('Output Data\PhixGenome.bin').abspath() 
#fileParsing.parseGenomeInt(path('Data\PhixGenome.fa').abspath(), binaryGenomeFile)

with open(binaryGenomeFile , 'rb') as f:
    decodedGenome = 0
    #decodedGenome = ''
    for chunk in iter(lambda: f.read(4), ''):
        fileBytes = struct.unpack('i', chunk)[0]     
        #print 'Unpacked data lambda: %s' % fileBytes
        decodedGenome += fileBytes #MUST CONCATENATE NOT ADD
        #decodedGenome += '{0:b}'.format(fileBytes) 
        #decodedGenome += format(fileBytes, 'b')
    print 'Integer decoding of data: %d' % decodedGenome            
    
#reads = fileParsing.parseReads(path('Data\HumanSequencingReads.tsv.bz2').abspath())        
reads = fileParsing.parseReads(path('Data\PhiXSequencingReads1000.fastq').abspath())
print next(reads)

#analyseAlignment.plotTimeVsMatches(reads, decodedGenome, path('Output Test Files\AlignmentAnalysis.png').abspath())

#matchesCount, totalCount = alignmentInt.alignHamming(reads, decodedGenome)
#print '%d/%d reads matched the genome.' % (matchesCount, totalCount) #The result is not 100% but this is to be expected due to sequencing errors. 

"""
textFile = path('Output Test Files\DataVisualisationTest.txt').abspath()
textFile = visualisation.visualisationText(genome, offsets, textFile)

jpgFile = path('Output Test Files\DataVisualisationTest.jpg').abspath()
jpgFile = visualisation.visualisationJPG(genome, offsets, jpgFile) 

pngFile = path('Output Test Files\DataVisualisationTest.png').abspath()
pngFile = visualisation.visualisationGD(genome, offsets, pngFile)
"""
#visualisation.visualisationTkinter(genome, offsets)
