"""
This file includes functions for efficient and inefficient parsing of a genome and sequencing reads.
"""

#To efficiently read the genome:
def readGenome(filename): #fastA
    file = open(filename, 'r')
    line = file.readline()
    if not line: #reached EOF
        return None

    sequenceLines = []
    while True: #read  sequence lines up to blank line
        line = file.readline().rstrip()
        if line == "": #reached end of record or end of file
            break
        sequenceLines.append(line)     
    sequence = ''.join(sequenceLines) #merge lines to form sequence
    
    file.close()
    return sequence

#To inefficiently read the genome:
def readGenomeInefficient(filename):
    genome = '' 
    with open(filename, 'r') as file: #opening a file for reading
        for line in file:
            if line[0] != '>': #ignore header line with genome information
                genome += line.rstrip() #add each line of bases to the string 
                #rstrip() removes any trailing whitespace from the ends of the string (trim off new line/tab/space)
    return genome
  
def removeNonAscii(string): #removes non-ASCII characters
    return ''.join(i for i in string if ord(i) < 128)

    #To efficiently read the sequencing reads:        
def readSequence(filename): #fastQ
    readID, sequence, quality = None, None, None
    file = open(filename, 'r')
    while True: #runs until EOF
        line = file.readline()
        if not line: #reached EOF
            break

        if line.startswith('@'): #first line of read/record
            #sequence = removeNonAscii(sequence)
            #yield readID, sequence, quality #where each loop iteration ends
            yield sequence

            #reset to default values
            readID = line.rstrip()
            sequence = None
            quality = None   

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
        
        elif not quality:
            quality = []
            while True: #collect base qualities
                quality += line.rstrip().replace(' ', '') 
                if len(quality) >= len(sequence): #bases and qualities line up
                    break
                else:
                    line = file.readline()
    
    file.close()
    #yield readID, sequence, quality  
    yield sequence
      
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
from itertools import groupby 
    
def readGenomeTry(filename):
    filehandle = open(filename, 'r')#, encoding="latin1")
    #ignore boolean (x[0]) and hold header or sequence since they alternate
    iteration = (x[1] for x in groupby(filehandle, lambda line: line[0] == ">"))
    for header in iteration:
        header.next()[1:].strip() #drop '>'
        seq = "".join(s.strip() for s in iteration.next()) #join all sequence lines
        yield seq        
        
def readGenomeSeqIO(filename):  
    genome = ''   
    #rU - open file for reading in universal readline mode (works across platforms due to differing newline characters)    
    filehandle = open(filename, 'rU')#, encoding="latin1")
    for record in SeqIO.parse(filehandle, "fasta"): #SeqRecord iterator - multiple records
    #record = SeqIO.read(filehandle, "fasta") #one record
        genome += record.seq #parses records one by one, without changing the file order
    filehandle.close()
    return genome

def readSequenceTry(filename):
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
"""