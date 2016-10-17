"""
This file contains the main functionality of the program.
"""

import FileParsing
import Alignment
import Visualisation

from path import path
import collections
import pdb

#pdb.set_trace()
#genome = FileParsing.readGenomeInefficient(path('Data\HumanGenome.fa.gz').abspath())
#reads = FileParsing.readSequence1(path('Data\HumanSequencingReads.tsv.bz2').abspath())
genome = FileParsing.readGenomeInefficient(path('Data\PhixGenome.fa').abspath())
reads = FileParsing.readSequence1(path('Data\PhiXSequencingReads1000.fastq').abspath())
 
"""
print "Length of the genome: %d" % len(genome)
print "Frequency of each base in the genome:", collections.Counter(genome)

print "Length of the reads: %d" % len(reads)
readsFreq = collections.Counter()
for read in reads:
    readsFreq.update(read)
print "Frequency of each base in the reads:", readsFreq
"""
matchesCount, totalCount, offsets = Alignment.align(reads, genome)
print "%d/%d reads matched the genome." % (matchesCount, totalCount)
#The result is not 100% but this is to be expected due to sequencing errors. 

#Problem: naiveApproxHamming() is returning empty lists

"""
textFile = path('Output Test Files\DataVisualisationTest.txt').abspath()
textFile = Visualisation.visualisationText(genome, offsets, textFile)

jpgFile = path('Output Test Files\DataVisualisationTest.jpg').abspath()
jpgFile = Visualisation.visualisationJPG(genome, offsets, jpgFile) 

pngFile = path('Output Test Files\DataVisualisationTest.png').abspath()
pngFile = Visualisation.visualisationGD(genome, offsets, pngFile)
"""

#Visualisation.visualisationTkinter(genome, offsets)
