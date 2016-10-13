"""
This file includes functions for efficient and inefficient parsing of a genome and sequencing reads.
"""

from Bio import SeqIO
from itertools import groupby
import csv
import pandas

#To efficiently read the genome:
def readGenome(filename):
    filehandle = open(filename, 'r')#, encoding="latin1")
    #ignore boolean (x[0]) and hold header or sequence since they alternate
    iteration = (x[1] for x in groupby(filehandle, lambda line: line[0] == ">"))
    for header in iteration:
        header.next()[1:].strip() #drop '>'
        seq = "".join(s.strip() for s in iteration.next()) #join all sequence lines
        yield seq        
        
#To read the genome using SeqIO:
def readGenomeSeqIO(filename):
    genome = ''   
    #rU - open file for reading in universal readline mode (works across platforms due to differing newline characters)    
    filehandle = open(filename, 'rU')#, encoding="latin1")
    for record in SeqIO.parse(filehandle, "fasta"): #SeqRecord iterator - multiple records
    #record = SeqIO.read(filehandle, "fasta") #one record
        genome += record.seq #parses records one by one, without changing the file order
    filehandle.close()
    return genome

#To inefficiently read the genome:
def readGenomeInefficient(filename):
    genome = '' 
    with open(filename, 'r') as file: #opening a file for reading
        for line in file:
            if line[0] != '>': #ignore header line with genome information
                genome += line.rstrip() #add each line of bases to the string 
                #rstrip() removes any trailing whitespace from the ends of the string (trim off new line/tab/space)
    return genome
  
#To efficiently read the sequencing reads:
def readSequence(filename):
    #filehandle = open(filename, 'rU', encoding="latin1")
    #tsvreader = csv.reader(filehandle, delimiter="\t")
    with open(filename) as file:     
        while True: #loops every 4 lines (each read is a set of 4) indefinitely until EOF
            next(file) #skip tag line 
            seq = file.readline().rstrip() #string of DNA bases
            next(file) #skip + line
            next(file) #skip quality sequence line          
            if len(seq) == 0: #reached EOF
                break
            yield seq
  
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
     