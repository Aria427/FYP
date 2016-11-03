#!/usr/bin/env python
#This file contains the main functionality of the program.

import fileParsing
import huffmanCompression
import alignment
import timeComplexity
import visualisation

from path import path
from time import time
import timeit
import pdb

#pdb.set_trace()

#textGenomeFile = path('Data\HumanGenome.txt').abspath()
#genome = fileParsing.parseGenome(path('Data\HumanGenome.fa.gz').abspath(), textGenomeFile)
textGenomeFile = path('Data\PhixGenome.txt').abspath() 
genome = fileParsing.parseGenome(path('Data\PhixGenome.fa').abspath(), textGenomeFile)

tree, codes = huffmanCompression.codeGeneration(genome)
print codes

binaryGenomeFile = path('Data\PhixGenome.bin').abspath() 
encodedGenome = huffmanCompression.encode(genome, binaryGenomeFile)

#with open(binaryGenomeFile, 'rb') as f:
    #byte = f.read(1)
    #genome = byte
    #while byte != '':
        #byte = f.read(1)
        #genome += byte
   
decodedGenome = huffmanCompression.decode(tree, encodedGenome)

#reads = fileParsing.parseReads(path('Data\HumanSequencingReads.tsv.bz2').abspath())        
reads = fileParsing.parseReads(path('Data\PhiXSequencingReads1000.fastq').abspath())
binaryReadsFile = path('Data\PhiXSequencingReads1000.bin').abspath()

t0 = time()
matchesCount, totalCount, offsets = alignment.alignFM(reads, decodedGenome) #, binaryReadsFile)
print "%d/%d reads matched the genome." % (matchesCount, totalCount) #The result is not 100% but this is to be expected due to sequencing errors. 
t1 = time()

print (t1-t0)

timeComplexity.plotTC(alignment.alignFM, reads, genome, 1)

"""
textFile = path('Output Test Files\DataVisualisationTest.txt').abspath()
textFile = visualisation.visualisationText(genome, offsets, textFile)

jpgFile = path('Output Test Files\DataVisualisationTest.jpg').abspath()
jpgFile = visualisation.visualisationJPG(genome, offsets, jpgFile) 

pngFile = path('Output Test Files\DataVisualisationTest.png').abspath()
pngFile = visualisation.visualisationGD(genome, offsets, pngFile)
"""
#visualisation.visualisationTkinter(genome, offsets, matches)
