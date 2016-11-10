#!/usr/bin/env python
#This file contains the main functionality of the program.

import fileParsing
import alignment
import analyseAlignment
import visualisation

from path import path
from bitarray import bitarray
import struct
import numpy 
import pdb

#pdb.set_trace()

#binaryGenomeFile = path('Output Data\HumanGenome.bin').abspath()
#fileParsing.parseGenome(path('Data\HumanGenome.fa.gz').abspath(), binaryGenomeFile)
binaryGenomeFile = path('Output Data\PhixGenome.bin').abspath() 
fileParsing.parseGenome(path('Data\PhixGenome.fa').abspath(), binaryGenomeFile)

#each line in genome file has length = 70

with open('Output Test Files\longlong.bin', 'wb') as f:
    l = fileParsing.baseToBinary('GAGTTTTATCGCTTCCATGACGCAGAAGTTA')
    print len(l)
    x = int(l, 2)
    print x.bit_length()
    f.write(struct.pack('q', x)) #long long - 8 bytes

with open('Output Test Files\int.bin', 'wb') as f:
    l = fileParsing.baseToBinary('GAGTTTTATCGCTTC')
    print len(l)
    x = int(l, 2)
    print x.bit_length()
    f.write(struct.pack('i', x)) #int - 4 bytes

with open('Output Test Files\int.bin' , 'rb') as f:
    fileBytes = f.read(4) #read first 4 bytes 
    data = fileBytes
    while fileBytes != '':
        fileBytes = f.read(4) #read the rest 4 bytes by 4 bytes
        data += fileBytes    
    for i in range(0, len(data), 4):
        fileBytes = struct.unpack('i', data[i:i+4])[0] #[0] to remove unnecessary syntax  
        print fileBytes

with open('Output Test Files\int.bin', 'rb') as f:
    for chunk in iter(lambda: f.read(4), ''):
        fileBytes = struct.unpack('i', chunk)[0]        
        print fileBytes
        
with open('Output Test Files\int.bin' , 'rb') as f:
    fileBytes = struct.unpack('i', f.read(4))[0]
    print fileBytes
    
with open('Output Test Files\int.bin', 'rb') as f:
    fileBytes = numpy.fromfile(f, dtype=numpy.int)
    print fileBytes
       
#ba = bitarray()
#with open(binaryGenomeFile, 'rb') as fh:
    #ba.fromfile(fh)

#decodedGenome = ''.join(ba.decode(fileParsing.bases))  

#reads = fileParsing.parseReads(path('Data\HumanSequencingReads.tsv.bz2').abspath())        
#reads = fileParsing.parseReads(path('Data\PhiXSequencingReads1000.fastq').abspath())

#analyseAlignment.plotTimeVsMatches(reads, decodedGenome, path('Output Test Files\AlignmentAnalysis.png').abspath())

#matchesCount, totalCount, offsets = alignment.alignHamming(reads, decodedGenome)
#print "%d/%d reads matched the genome." % (matchesCount, totalCount) #The result is not 100% but this is to be expected due to sequencing errors. 

"""
textFile = path('Output Test Files\DataVisualisationTest.txt').abspath()
textFile = visualisation.visualisationText(genome, offsets, textFile)

jpgFile = path('Output Test Files\DataVisualisationTest.jpg').abspath()
jpgFile = visualisation.visualisationJPG(genome, offsets, jpgFile) 

pngFile = path('Output Test Files\DataVisualisationTest.png').abspath()
pngFile = visualisation.visualisationGD(genome, offsets, pngFile)
"""
#visualisation.visualisationTkinter(genome, offsets)
