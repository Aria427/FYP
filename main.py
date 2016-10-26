"""
This file contains the main functionality of the program.
"""

import fileParsing
import alignment
import visualisation

from path import path
import collections
import pdb

#pdb.set_trace()
genome = fileParsing.parseGenome(path('Data\HumanGenome.fa.gz').abspath(), path('Data\HumanGenome.txt').abspath())
#reads = fileParsing.parseReads(path('Data\HumanSequencingReads.tsv.bz2').abspath())
#genome = fileParsing.parseGenome(path('Data\PhixGenome.fa').abspath(), path('Data\PhixGenome.txt').abspath())
#reads = fileParsing.parseReads(path('Data\PhiXSequencingReads1000.fastq').abspath())

#matchesCount, totalCount, offsets, matches = alignment.align(reads, genome)
#print "%d/%d reads matched the genome." % (matchesCount, totalCount) #The result is not 100% but this is to be expected due to sequencing errors. 

"""
textFile = path('Output Test Files\DataVisualisationTest.txt').abspath()
textFile = visualisation.visualisationText(genome, offsets, textFile)

jpgFile = path('Output Test Files\DataVisualisationTest.jpg').abspath()
jpgFile = visualisation.visualisationJPG(genome, offsets, jpgFile) 

pngFile = path('Output Test Files\DataVisualisationTest.png').abspath()
pngFile = visualisation.visualisationGD(genome, offsets, pngFile)
"""
#visualisation.visualisationTkinter(genome, offsets, matches)
