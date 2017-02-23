#!/usr/bin/env python
#This file includes functions for aligning reads with overlapping regions of the reference genome. 

import fileCompressionAndParsing
import alignmentString
import gzip
import numpy as np
from path import path

#Human genome = 3,273,481,150 bytes
#Human sequencing reads = 28,094,847 reads

#This function reads the file in chunks.
def readInChunks(file, chunkSize=65536):
   while True:
       data = file.read(chunkSize).rstrip().upper().replace('N', '').replace(' ', '')
       if not data:
           break
       yield data

#This function also reads the file in chunks.       
def readInChunk(file, chunkSize):
    return iter(lambda: file.read(chunkSize).rstrip().upper().replace('N', '').replace(' ', ''), '')     

#The genome is read in chunks while storing overlapping regions between each chunk.
#The size of the overlap is the length of a read; i.e. 60.    
    
#This function aligns the sequencing reads against the uncompressed reference genome.
def alignUncompressed(reads, genome):
    totalMatches, totalCount, totalOffsets = 0, 0, []
    with open(genome, 'r') as f:
        #subseqs = readInChunk(f, 100000) #chunk = 100,000 bytes
        subseqs = (line.rstrip().upper().replace('N', '').replace(' ','')
                    for line in f if line and line[0] != '>') #ignore header line
        overlap = '' 
        chunkCount = 0
        for s in subseqs:
            print 'Chunk %d read' % chunkCount
            s = overlap + s #overlap is appended to the start of the next chunk
            
            #r = fileCompressionAndParsing.parseReadsString(reads) 
            r = fileCompressionAndParsing.parseReadsPhiXString(reads)
            matchesCount, count, offsets = alignmentString.alignFM(r, s)
            
            totalMatches += matchesCount
            totalCount = count #number of reads stays the same as every chunk goes through each read again
            totalOffsets.append(offsets)
            overlap = s[-100:] #100 for PhiX, 60 for Human
            print 'Chunk %d aligned' % chunkCount        
            chunkCount += 1
            print '%d/%d reads matched the genome.' % (totalMatches, totalCount)
            #print totalOffsets

#This function aligns the sequencing reads against the compressed reference genome.  
def alignCompressed(reads, genome):
    tempFile = path('Output Test Files/tmp.txt').abspath()
    totalMatches, totalCount, totalOffsets = 0, 0, []
    with open(genome, 'r') as f:
        fileCompressionAndParsing.deparseGenomeInt(genome, tempFile)
    
    with open(tempFile, 'r') as f:
        subseqs = (line.rstrip().upper().replace('N', '').replace(' ','')
                    for line in f if line and line[0] != '>') #ignore header line   
        overlap = '' 
        chunkCount = 0
        #integerCount = 0
        #for i in integerGenome: 
            #integerCount += 1
            #s = fileCompressionAndParsing.decompressInt(i) #integer is decompressed
            #for i in integerGenome[integerCount:]:
            #    while len(s) < 100:
            #        s += fileCompressionAndParsing.decompressInt(i)
        for s in subseqs:  
            print 'Chunk %d read' % chunkCount
            s = overlap + s #overlap is appended to the start of the next chunk
            
            #r = fileCompressionAndParsing.parseReadsString(reads) 
            r = fileCompressionAndParsing.parseReadsPhiXString(reads)
            matchesCount, count, offsets = alignmentString.alignFM(r, s)
            
            totalMatches += matchesCount
            totalCount = count #number of reads stays the same as every chunk goes through each read again
            totalOffsets.append(offsets)
            overlap = s[-100:]
            print 'Chunk %d aligned' % chunkCount 
            chunkCount += 1
            print '%d/%d reads matched the genome.' % (totalMatches, totalCount)
            #print totalOffsets

