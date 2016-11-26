#!/usr/bin/env python
#This file includes functions for efficient parsing of a genome and sequencing reads.

import gzip
import bz2
 
import random
import struct
from bitarray import bitarray

binaryBases = {'A' : '00',
               'C' : '01',
               'G' : '10',
               'T' : '11',
               }

bases = {'A' : bitarray('00'),
         'C' : bitarray('01'),
         'G' : bitarray('10'),
         'T' : bitarray('11') }               
               
def baseToBinary(line):
    for base, binary in binaryBases.items():
        line = line.replace(base, binary)
    return line             

def binaryToBase(line):
    for base, binary in binaryBases.items():
        line = line.replace(binary, base)
    return line      
  
def compressGenome(line, lineLength, binaryFile):
    for i in range(0, lineLength+1):
        l = line[i:i+15].rstrip().upper().replace('N', '').replace(' ', '')  #15 = allowed amount for int 
        l = baseToBinary(l)
        try:
            bytes = int(l, 2) #create 4 bytes from base 2 integer
            binaryFile.write(struct.pack('=i', bytes))
        except ValueError:
            pass
            #print "Invalid string found in bytes: %s" % format(bytes) 

def compressReads(line, lineLength):
    subsequences = [] #array of integers
    for i in range(0, lineLength+1):
        l = line[i:i+15].rstrip().upper().replace('N', '').replace(' ', '')  #15 = allowed amount for int 
        l = baseToBinary(l)
        try:
            bytes = int(l, 2) #create 4 bytes from base 2 integer
            subsequences.append(bytes)
        except ValueError:
            pass
            #print "Invalid string found in bytes: %s" % format(bytes)  
    #print subsequences
    #s = reduce(lambda x,y: x+str(y), subsequences, '')
    #sequence = int(s)
    #sequence = int(''.join(map(str, subsequences))) #join integer array into one single integer
    return subsequences#sequence
                       
#To read the genome into an integer:
def parseGenomeInt(input, output): #file is compressed by ~70%
    binary = open(output, 'wb')
    with open(input, 'r') as file: 
    #with gzip.open(input, 'r') as file:
        for line in file:
            if line and line[0] != '>': #ignore header line with genome information
                #Each line has length = 70 in PhiX and length = 50 in Human
                #Last line has length = 67 in PhiX and length = 41 in Human
                #print len(line.rstrip()) 
                compressGenome(line, 50, binary)
                #pass
        #last = line
        #print len(last)
    binary.close()  
        
#To read the genome into a bitarray:
def parseGenomeBitArray(input, output): #file is compressed by ~75%
    binary = open(output, 'wb')
    with open(input, 'r') as file: 
    #with gzip.open(input, 'r') as file:
        for line in file:
            if line and line[0] != '>': #ignore header line with genome information
                l = line.rstrip().upper().replace('N', '') #N => no confidence
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
    
#To parse the sequencing reads into an integer generator:        
def parseReadsInt(filename):  
    readID, sequence, quality = '', '', ''
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
            #sequenceLines = []
            while not line.startswith('+'): #not placeholder line (third line)
                #Each line has length = 123 in Human
                sequenceLine = compressReads(line, 123)
                #sequenceLines.append(sequenceLine) #no whitespace in integer sequence
                line = file.readline()
            #s = reduce(lambda x,y: x+str(y), sequenceLines, '')
            #sequence = int(s)
            #sequence = int(''.join(map(str, sequenceLines))) #merge lines to form sequence
            line = file.readline()
            yield sequenceLine#sequence
        
        elif not quality:
            quality = []
            while True: #collect base qualities
                quality += line.rstrip().replace(' ', '') 
                if len(quality) >= len(sequence): #bases and qualities line up
                    break
                else:
                    line = file.readline()
    file.close() 
       
#To parse the sequencing reads into a string generator:        
def parseReadsString(filename):  
    readID, sequence, quality = '', '', ''
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
                line = line.rstrip().upper().replace('N', '').replace(' ', '')
                sequenceLines.append(line) #no whitespace in string sequence
                line = file.readline()
            sequence = ''.join(sequenceLines) #merge lines to form sequence
            yield sequence
        
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