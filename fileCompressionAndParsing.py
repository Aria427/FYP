#!/usr/bin/env python
#This file includes functions for efficient compression and parsing of a genome and sequencing reads.

import gzip
import bz2
 
import struct
import numpy as np
from bitarray import bitarray

binaryBases = {'A' : '00',
               'C' : '01',
               'G' : '10',
               'T' : '11' }

bases = {'A' : bitarray('00'),
         'C' : bitarray('01'),
         'G' : bitarray('10'),
         'T' : bitarray('11') }               
 
#This function converts each nucleotide base to its binary representation.              
def baseToBinary(line):
    for base, binary in binaryBases.items():
        line = line.replace(base, binary)
    return line             

#This function converts each binary encoding to its equivalent nucleotide base.
def binaryToBase(line):
    result = ''
    bins = [''.join(t) for t in zip(*[iter(line)]*2)] #seperate each two bins in a list

    for b in bins: #for each two bits
        for base, binary in binaryBases.items():
            if binary == b:
                result += b.replace(b, base)

    return result     
    
#This function compresses a sequence using int's into a binary file.      
def compressInt(sequence, binaryFile):
    s = baseToBinary(sequence)
    try:
        byte = int(s, 2) #create 4 bytes from base 2 integer
        binaryFile.write(struct.pack('=i', byte)) #= => native byte-order
    except ValueError:
        pass
        #print "Invalid string found in bytes: %s" % format(byte) 

#This function compresses a sequence using long's into a binary file.            
def compressLong(sequence, binaryFile): 
    s = baseToBinary(sequence)
    try:
        byte = long(s, 2) #create 8 bytes from base 2 integer
        binaryFile.write(struct.pack('=q', byte)) #= => native byte-order
    except ValueError:
        pass      
                       
#This function reads the genome sequence into an integer using compressInt().
def parseGenomeInt(input, output): #file is compressed by ~70%
    binary = open(output, 'wb')
    with gzip.open(input, 'r') as f:
        #read genome into line-by-line generator
        subseqs = (line.rstrip().upper().replace('N', '').replace(' ','')
                    for line in f if line and line[0] != '>') #ignore header line
        #Each line has length = 70 in PhiX and length = 50 in Human
        #Last line has length = 67 in PhiX and length = 41 in Human
        for seq in subseqs:
            s = seq[0:15] #15 = allowed amount for int 
            compressInt(s, binary)
            s = seq[15:30] #15->29, 15 included 
            compressInt(s, binary)
            s = seq[30:45] #30->44, 30 included
            compressInt(s, binary)
            s = seq[45:51] #45->50, 45 included
            compressInt(s, binary)
            #s = seq[60:71] #60->70, 60 included
            #compressInt(s, binary)
    binary.close()
    
#This function reads the genome sequence into a long using compressLong().
def parseGenomeLong(input, output): #file is compressed by ~70%
    binary = open(output, 'wb')
    with gzip.open(input, 'r') as f:
        #read genome into line-by-line generator
        subseqs = (line.rstrip().upper().replace('N', '').replace(' ','')
                    for line in f if line and line[0] != '>') #ignore header line
        #Each line has length = 70 in PhiX and length = 50 in Human
        #Last line has length = 67 in PhiX and length = 41 in Human
        for seq in subseqs:
            s = seq[0:15] #15 = allowed amount for long
            compressLong(s, binary)
            s = seq[15:30] #15->29, 15 included 
            compressLong(s, binary)
            s = seq[30:45] #30->44, 30 included
            compressLong(s, binary)
            s = seq[45:51] #45->50, 45 included
            compressLong(s, binary)
            #s = seq[60:71] #60->70, 60 included
            #compressLong(s, binary)
    binary.close()
 
#This function decompresses an integer sequence into a text file (if a filename is given).  
def decompressInt(sequence, textFile=None): 
    integer = struct.unpack('=i', sequence) 
    integer = integer[0] #as unpack returns a tuple
    
    binary = '{0:08b}'.format(integer) #convert integer to binary format
    seq = binaryToBase(binary)
    
    if textFile != None:
        textFile.write(seq)
    else:
        return seq
    
#This function decompresses a long sequence into a text file (if a filename is given).  
def decompressLong(sequence, textFile=None):
    longS = struct.unpack('=q', sequence) 
    longS = longS[0] #as unpack returns a tuple
    
    binary = '{0:08b}'.format(longS) #convert long to binary format
    seq = binaryToBase(binary)
    
    if textFile != None:
        textFile.write(seq)
    else:
        return seq
 
#This function reads the integer compressed genome sequence back into its original format.
def deparseGenomeInt(input, output):
    text = open(output, 'w')
    with open(input , 'rb') as f:
        integerGenome = np.fromfile(f, dtype=np.int32) #read genome into list of ints
        integers = 0
        for i in integerGenome:
            decompressInt(i, text)
            integers += 1
            if (integers % 5) == 0: #4 for human, 5 for PhiX
                text.write('\n')
    text.close()
    
#This function reads the long compressed genome sequence back into its original format.
def deparseGenomeLong(input, output):
    text = open(output, 'w')
    with open(input , 'rb') as f:
        longGenome = np.fromfile(f, dtype=np.int64) #read genome into list of longs
        longs = 0
        for l in longGenome:
            decompressLong(l, text)
            longs += 1
            if (longs % 4) == 0:
                text.write('\n')
    text.close()    
    
#This function reads the genome sequence into a bitarray.
#Ambiguities arise when a subsequence does not completely fill the bit array;
#i.e. 0's are appended => additional A's are encoded
def parseGenomeBitArray(input, output): #file is compressed by ~75%
    binary = open(output, 'wb')
    with open(input, 'r') as file:
        for line in file:
            if line and line[0] != '>': #ignore header line with genome information
                l = line.rstrip().upper().replace('N', '') #N => no confidence
                ba = bitarray()
                ba.encode(bases, l)
                ba.tofile(binary)
    binary.close()    
  
#This function reads the genome into a string.
#It leads to a memory error for the human genome.
def parseGenomeString(input): #file is not compressed
    genome = ''
    with open(input, 'r') as file: 
    #with gzip.open(input, 'r') as file:
        for line in file:
            if line and line[0] != '>': #ignore header line with genome information
                l = line.rstrip().upper().replace('N', '').replace(' ', '') 
                genome += l 
    return genome

#This function compresses the sequencing reads in the same manner as the genome sequence is compressed.            
def compressReadsInt(line, lineLength):
    subsequences = [] #array of integers
    for i in range(0, lineLength+1):
        l = line[i:i+15].rstrip().upper().replace('N', '') #15 = allowed amount for int 
        l = baseToBinary(l)
        try:
            bytes = int(l, 2) #create 4 bytes from base 2 integer
            subsequences.append(bytes)
        except ValueError:
            pass
            #print "Invalid string found in bytes: %s" % format(bytes)  
    #print subsequences
    #s = reduce(lambda x,y: x+str(y), subsequences, '')
    #sequence = int(s) #join integer array into one single integer
    return subsequences#sequence    
    
#This function parses the Human sequencing reads into an integer generator.   
def parseReadsInt(filename):  
    flag, sequence, quality = '', '', ''
    file = bz2.BZ2File(filename, 'r') 
    while True: #runs until EOF
        line = file.readline() 
        if not line: #reached EOF
            break
        
        if line.startswith('#'): #read details
            line = file.readline()
            pass
        
        elif line.startswith('>'): #>flags reads scores
            line = file.readline()
            pass
        
        else:
            line = line.split()
            flag = line[0]
            sequence = line[1]
            quality = line[2]
            
            if sequence == 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN':
                line = file.readline()
                pass
            
            else:
                #Each read has length = 60
                sequence = compressReadsInt(sequence, 60)
                #seq = reduce(lambda x,y: x+str(y), sequence, '')
                #s = int(seq) #merge lines to form sequence
                yield sequence
        
    file.close()     
 
#This function parses the Human sequencing reads into a string generator.   
def parseReadsString(filename):  
    flag, sequence, quality = '', '', ''
    file = bz2.BZ2File(filename, 'r') 
    while True: #runs until EOF
        line = file.readline() 
        if not line: #reached EOF
            break
        
        if line.startswith('#'): #read details
            line = file.readline()
            pass
        
        elif line.startswith('>'): #>flags reads scores
            line = file.readline()
            pass
        
        else:
            line = line.split()
            flag = line[0]
            sequence = line[1]
            quality = line[2]
            
            #Each read has length = 60
            if sequence == 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN':
                line = file.readline()
                pass
            
            else:
                qualityNoNs = ''
                N = [pos for pos, char in enumerate(sequence) if char == 'N'] #positions of N in read
                sequenceNoNs = sequence.rstrip().upper().replace('N', '').replace(' ', '') 
                
                for i in range(len(quality)):
                    if i not in N: #remove indices corresponding to Ns in read
                        qualityNoNs = qualityNoNs + quality[i]
                
                yield sequenceNoNs, qualityNoNs
        
    file.close()       
    
#This function parses the PhiX sequencing reads into an integer generator.        
def parseReadsPhiXInt(filename):  
    readID, sequence, quality = '', '', ''
    file = open(filename, 'r')
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
                #Each line has length = 123 in Human
                sequenceLine = compressReadsInt(line, 123)
                sequenceLines.append(sequenceLine) #no whitespace in integer sequence
                line = file.readline()
            #s = reduce(lambda x,y: x+str(y), sequenceLine, '')
            #sequence = int(s) #merge lines to form sequence
            line = file.readline()
            yield sequenceLine
        
        elif not quality:
            quality = []
            while True: #collect base qualities
                quality += line.rstrip().replace(' ', '') 
                #if len(quality) >= len(sequence): #bases and qualities line up
                    #break
                #else:
                line = file.readline()
                
    file.close() 
       
#This function parses the Phix sequencing reads into a string generator.       
def parseReadsPhiXString(filename):  
    readID, sequence, quality = '', '', ''
    file = open(filename, 'r')
    
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
        
        elif not sequence or not quality:
            sequenceLines = [] 
            while not line.startswith('+'): #not placeholder line (third line)
                #rstrip() - removes leading/trailing whitespace
                #replace() - removes whitespace from within string
                N = [pos for pos, char in enumerate(sequence) if char == 'N'] #positions of N in read
                line = line.rstrip().upper().replace(' ', '')
                sequenceLines.append(line) #no whitespace in string sequence
                line = file.readline()
            sequence = ''.join(sequenceLines) #merge lines to form original sequence
            sequenceNoNs = sequence.replace('N', '') #remove Ns
            temp = sequenceNoNs
        
            qualityLines = []
            qualityNoNs = ''
            
            while True: #collect base qualities
                line = line.rstrip().replace(' ', '')
                qualityLines.append(line) 
                quality = ''.join(qualityLines) #merge lines to form quality
                if len(quality) >= len(sequence): #bases and qualities line up
                    break
                else:
                    line = file.readline()
            
            for i in range(len(quality)): 
                if i not in N: #remove indices corresponding to Ns in read
                    qualityNoNs = qualityNoNs + quality[i]
                
            yield temp, qualityNoNs
        
    file.close() 
    