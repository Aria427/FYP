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
fileParsing.parseGenomeInt(path('Data\PhixGenome.fa').abspath(), binaryGenomeFile)


with open('Output Test Files\longlong.bin', 'wb') as f:
    sequence = 'GAGTTTTATCGCTTCCATGACGCAGAAGTTA'
    print 'Length of long long sequence: %d' % len(sequence)
    binary = fileParsing.baseToBinary(sequence)
    print 'Length of binary converted long long sequence: %d' % len(binary)
    x = int(binary, 2)
    print 'Bit length of binary converted long long sequence: %d' % x.bit_length()
    f.write(struct.pack('q', x)) #long long - 8 bytes

with open('Output Test Files\int.bin', 'wb') as f:
    sequence = 'GAGTTTTATCGCTTC'
    print 'Length of int sequence: %d' % len(sequence)
    binary = fileParsing.baseToBinary(sequence)
    print 'Length of binary converted int sequence: %d' % len(binary)
    x = int(binary, 2)
    print 'Bit length of binary converted int sequence: %d' % x.bit_length()
    f.write(struct.pack('i', x)) #int - 4 bytes

print '\n'

with open('Output Test Files\int.bin' , 'rb') as f:
    fileBytes = f.read(4) #read first 4 bytes 
    data = fileBytes
    while fileBytes != '':
        fileBytes = f.read(4) #read the rest 4 bytes by 4 bytes
        data += fileBytes 
    print 'Encoded data: %s' % data
    print 'Representation of encoded data: %s' % repr(data)
    print 'Hexadecimal format of encoded data: %s' % data.encode('hex')
    print 'Integer format of encoded data: %s' % ' '.join([str(ord(a)) for a in data])
    for i in range(0, len(data), 4):
        fileBytes = struct.unpack('i', data[i:i+4])[0] #[0] to remove unnecessary syntax  
        print 'Unpacked data range: %s' % fileBytes

with open('Output Test Files\int.bin' , 'rb') as f:
    fileBytes = struct.unpack('i', f.read(4))[0]
    print 'Unpacked data direct: %s' % fileBytes       
        
with open('Output Test Files\int.bin', 'rb') as f:
    for chunk in iter(lambda: f.read(4), ''):
        fileBytes = struct.unpack('i', chunk)[0]        
        print 'Unpacked data lambda: %s' % fileBytes
    
with open('Output Test Files\int.bin', 'rb') as f:
    fileBytes = numpy.fromfile(f, dtype=numpy.int)
    print 'Unpacked data numpy: %s' % fileBytes
 
      
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
