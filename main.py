#!/usr/bin/env python
#This file contains the main functionality of the program.

import fileParsing
import alignment
import visualisation

from path import path
import collections
import pdb

#pdb.set_trace()

#binaryGenomeFile = path('Data\HumanGenome.bin').abspath()
#only need to parseGenome() once; use binary file directly afterwards
#genome = fileParsing.parseGenome(path('Data\HumanGenome.fa.gz').abspath(), binaryGenomeFile)

binaryGenomeFile = path('Data\PhixGenome.bin').abspath() 
#genome = fileParsing.parseGenome(path('Data\PhixGenome.fa').abspath(), binaryGenomeFile)

with open(binaryGenomeFile, 'rb') as f:
    byte = f.read(1)
    genome = byte
    while byte != '':
        byte = f.read(1)
        genome += byte
       
#reads = fileParsing.parseReads(path('Data\HumanSequencingReads.tsv.bz2').abspath())        
reads = fileParsing.parseReads(path('Data\PhiXSequencingReads1000.fastq').abspath())

matchesCount, totalCount, offsets, matches = alignment.align(reads, genome)
print "%d/%d reads matched the genome." % (matchesCount, totalCount) #The result is not 100% but this is to be expected due to sequencing errors. 

"""
textFile = path('Output Test Files\DataVisualisationTest.txt').abspath()
textFile = visualisation.visualisationText(genome, offsets, textFile)

jpgFile = path('Output Test Files\DataVisualisationTest.jpg').abspath()
jpgFile = visualisation.visualisationJPG(genome, offsets, jpgFile) 

pngFile = path('Output Test Files\DataVisualisationTest.png').abspath()
pngFile = visualisation.visualisationGD(genome, offsets, pngFile)
"""
#visualisation.visualisationTkinter(genome, offsets, matches)
