#!/usr/bin/env python
#This file includes functions for efficient parsing of a genome and sequencing reads.

import gzip
import bz2
 
import struct
from bitarray import bitarray

binaryBases = {'A' : '00',
               'C' : '01',
               'G' : '10',
               'T' : '11' }

bases = {'A' : bitarray('00'),
         'C' : bitarray('01'),
         'G' : bitarray('10'),
         'T' : bitarray('11') }               
               
def baseToBinary(line):
    for base, binary in binaryBases.items():
        line = line.replace(base, binary)
    return line             
 
#To read the genome into an integer:
def parseGenomeInt(input, output):
    binary = open(output, 'wb')
    with open(input, 'r') as file: 
    #with gzip.open(input, 'r') as file:
        for line in file:
            if line and line[0] != '>': #ignore header line with genome information
                l = line[0:15].rstrip().upper().replace('N', '') #15 = allowed amount for int 
                l = baseToBinary(l)
                x = int(l, 2)
                binary.write(struct.pack('i', x))
                l = line[15:30].rstrip().upper().replace('N', '') #15->29, 15 included 
                l = baseToBinary(l)
                x = int(l, 2)
                binary.write(struct.pack('i', x))
                l = line[30:45].rstrip().upper().replace('N', '') #30->44, 30 included
                l = baseToBinary(l)
                x = int(l, 2)
                binary.write(struct.pack('i', x))
                l = line[45:60].rstrip().upper().replace('N', '') #45->59, 40 included
                l = baseToBinary(l)
                x = int(l, 2)
                binary.write(struct.pack('i', x))
                #each line has length = 70
                l = line[60:71].rstrip().upper().replace('N', '') #60->70, 60 included
                l = baseToBinary(l)
                x = int(l, 2)
                binary.write(struct.pack('i', x))
    binary.close()  
        
#To read the genome into a bitarray:
def parseGenomeBitArray(input, output): #file is compressed by ~200,000KB
    binary = open(output, 'wb')
    with open(input, 'r') as file: 
    #with gzip.open(input, 'r') as file:
        for line in file:
            if line and line[0] != '>': #ignore header line with genome information
                l = line.rstrip().upper().replace('N', '') 
                ba = bitarray()
                ba.encode(bases, l)
                ba.tofile(binary)
    binary.close()    
  
#To read the genome into a string:
def parseGenomeString(input, output): #file is not compressed
    genome = ''
    with open(input, 'r') as file: 
    #with gzip.open(input, 'r') as file:
        for line in file:
            if line and line[0] != '>': #ignore header line with genome information
                l = line.rstrip().upper().replace('N', '') 
                genome += l 
    return genome
    
#To efficiently read the sequencing reads:        
def parseReads(filename): #fastQ 
    readID, sequence, quality = '', '', ''
    ba = bitarray()
    file = open(filename, 'r')
    #file = bz2.BZ2File(filename, 'r')
    while True: #runs until EOF
        line = file.readline()
        if not line: #reached EOF
            break

        if line.startswith('@'): #first line of read/record
            #reset to default values
            readID = line.rstrip()
            sequence = ''
            quality = ''   

        elif not readID: #if no previous line starts with @
            readID = line.rstrip() #get first ID
            continue

        elif not sequence:
            sequenceLines = []
            while not line.startswith('+'): #not placeholder line (third line)
                #rstrip() - removes leading/trailing whitespace
                #replace() - removes whitespace from within string
                sequenceLines.append(line.rstrip().replace(' ', '')) #no whitespace in sequence
                line = file.readline()
            sequence = ''.join(sequenceLines) #merge lines to form sequence
            ba.encode(bases, sequence)
            yield sequence
            #yield ba #ValueError: symbol not in prefix code
        
        elif not quality:
            quality = []
            while True: #collect base qualities
                quality += line.rstrip().replace(' ', '') 
                if len(quality) >= len(sequence): #bases and qualities line up
                    break
                else:
                    line = file.readline()
    
    file.close() 
       

"""
The following functions were also tested but proved to be insufficient: 

from itertools import groupby   
    
def readGenome1(filename): #fastA 
    sequence = None
    sequenceLines = []
    file = open(filename, 'r')
    #file = gzip.open(filename, 'r')
    line = file.readline()
    if not line: #reached EOF
        return None

    while True: #read  sequence lines up to blank line
        line = file.readline().rstrip()
        if line == "": #reached end of record or end of file
            break
        sequenceLines.append(line)     
    sequence = ''.join(sequenceLines) #merge lines to form sequence
    
    file.close()
    return sequence

def readGenome2(filename): #faster than readGenome1()        
    #filehandle = open(filename, 'r')
    filehandle = gzip.open(filename, 'r')
    fhBuffer = io.BufferedReader(filehandle)
    #ignore boolean (x[0]) and hold header or sequence since they alternate
    iteration = (x[1] for x in groupby(fhBuffer, lambda line: line[0] == ">"))
    for header in iteration:
        header.next()[1:].strip() #drop '>'
        seq = ''.join(s.strip() for s in iteration.next()) #join all sequence lines
        yield seq    
    filehandle.close()    

def readSequence1(filename):
    sequenceLines = []
    with open(filename) as file:  
    #with bz2.BZ2File(filename) as file:
        while True: #loops every 4 lines (each read is a set of 4) indefinitely until EOF
            next(file) #skip tag line 
            seq = next(file) #string of DNA bases
            next(file) #skip + line
            next(file) #skip quality sequence line          
            if len(seq) == 0: #reached EOF
                break
            sequenceLines.append(seq)
            yield seq
        sequences = ''.join(sequenceLines)
    yield sequences
                     
def readSequence2(filename):
    sequences = []
    with open(filename) as file:  
        while True: #loops every 4 lines (each read is a set of 4) indefinitely until EOF
            file.readline() #read tag
            seq = file.readline().rstrip() #string of DNA bases
            file.readline() #+
            file.readline() #quality sequence - each quality score is Phredd33 encoded (converts to ASCII character)
            if len(seq) == 0: #reached EOF
                break
            sequences.append(seq)
    return sequences     

"""