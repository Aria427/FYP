"""
This file includes functions for efficient and inefficient parsing of a genome and sequencing reads.
"""

from itertools import groupby
import gzip
import bz2
import io
import subprocess
import time

#To efficiently read the genome:
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
    #ignore boolean (x[0]) and hold header or sequence since they alternate
    iteration = (x[1] for x in groupby(filehandle, lambda line: line[0] == ">"))
    for header in iteration:
        header.next()[1:].strip() #drop '>'
        seq = ''.join(s.strip() for s in iteration.next()) #join all sequence lines
        yield seq           
        
def readGenome3(filename):
    genome = '' 
    #with open(filename, 'r') as file: #opening a file for reading
    with gzip.open(filename, 'r') as gzipFile:
        with io.BufferedReader(gzipFile) as file:
            for line in file:
                if line[0] != '>': #ignore header line with genome information
                    genome += line.rstrip() #add each line of bases to the string 
                    #rstrip() removes any trailing whitespace from the ends of the string (trim off new line/tab/space)
    return genome

#To efficiently read the sequencing reads:        
def readSequence1(filename): #fastQ 
    readID, sequence, quality = '', '', ''
    #file = open(filename, 'r')
    file = bz2.BZ2File(filename, 'r')
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
 
def readSequence2(filename):
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
             
#To inefficiently read the sequencing reads:          
def readSequenceInefficient(filename):
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
The following functions were tested but proved to be insufficient:   

from Bio import SeqIO
                   
def readGenomeSeqIO(filename): #keeps looping even for PhiX  
    genome = ''   
    #rU - open file for reading in universal readline mode (works across platforms due to differing newline characters)    
    filehandle = open(filename, 'rU')#, encoding="latin1")
    for record in SeqIO.parse(filehandle, "fasta"): #SeqRecord iterator - multiple records
    #record = SeqIO.read(filehandle, "fasta") #one record
        genome += record.seq #parses records one by one, without changing the file order
    filehandle.close()
    return genome
"""