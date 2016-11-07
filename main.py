#!/usr/bin/env python
#This file contains the main functionality of the program.

import fileParsing
import alignment
import analyseAlignment
import visualisation

from path import path
from bitarray import bitarray
import pdb

#pdb.set_trace()

#binaryGenomeFile = path('Output Data\HumanGenome.bin').abspath()
#fileParsing.parseGenome(path('Data\HumanGenome.fa.gz').abspath(), binaryGenomeFile)
binaryGenomeFile = path('Output Data\PhixGenome.bin').abspath() 
fileParsing.parseGenome(path('Data\PhixGenome.fa').abspath(), binaryGenomeFile)

ba = bitarray()
with open(binaryGenomeFile, 'rb') as fh:
    ba.fromfile(fh)

decodedGenome = ''.join(ba.decode(fileParsing.bases))  

#reads = fileParsing.parseReads(path('Data\HumanSequencingReads.tsv.bz2').abspath())        
reads = fileParsing.parseReads(path('Data\PhiXSequencingReads1000.fastq').abspath())

#analyseAlignment.plotTimeVsMatches(reads, decodedGenome, path('Output Test Files\AlignmentAnalysis.png').abspath())

matchesCount, totalCount, offsets = alignment.alignHamming(reads, decodedGenome)
print "%d/%d reads matched the genome." % (matchesCount, totalCount) #The result is not 100% but this is to be expected due to sequencing errors. 

"""
textFile = path('Output Test Files\DataVisualisationTest.txt').abspath()
textFile = visualisation.visualisationText(genome, offsets, textFile)

jpgFile = path('Output Test Files\DataVisualisationTest.jpg').abspath()
jpgFile = visualisation.visualisationJPG(genome, offsets, jpgFile) 

pngFile = path('Output Test Files\DataVisualisationTest.png').abspath()
pngFile = visualisation.visualisationGD(genome, offsets, pngFile)
"""
#visualisation.visualisationTkinter(genome, offsets)
