"""

"""

from Bio import SeqIO
from itertools import groupby
import csv
import pandas

#To read the genome efficiently:
def readGenome(filename):
    filehandle = open(filename, 'r', encoding="latin1")
    faiter = (x[1] for x in groupby(filehandle, lambda line: line[0] == ">"))
    for header in faiter:
        header.__next__()[1:].strip() #drop '>'
        seq = "".join(s.strip() for s in faiter.__next__()) #join all sequence lines
        yield(seq)

#To read the genome using SeqIO:
def readGenomeSeqIO(filename):
    genome = ''   
    #rU - open file for reading in universal readline mode (works across platforms due to differing newline characters)    
    filehandle = open(filename, 'rU', encoding="latin1")
    for record in SeqIO.parse(filehandle, "fasta"): #SeqRecord iterator - multiple records
    #record = SeqIO.read(filehandle, "fasta") #one record
        genome += record.seq #parses records one by one, without changing the file order
    filehandle.close()
    return genome

#To inefficiently read the genome:
def readGenomeInefficient(filename):
    genome = '' 
    with open(filename, 'r', encoding="latin1") as file: #opening a file for reading
        for line in file:
            if line[0] != '>': #ignore header line with genome information
                genome += line.rstrip() #add each line of bases to the string 
                #rstrip() removes any trailing whitespace from the ends of the string (trim off new line/tab/space)
    return genome
  
#To read and parse the reads:
def readSequence(filename):
    sequences = []
    #filehandle = open(filename, 'rU', encoding="latin1")
    #tsvreader = csv.reader(filehandle, delimiter="\t")
    with open(filename, encoding="latin1") as file:     
        while True: #loops every 4 lines (each read is a set of 4) indefinitely until EOF
            file.readline() #read tag
            seq = file.readline().rstrip() #string of DNA bases
            file.readline() #+
            file.readline() #quality sequence - each quality score is Phredd33 encoded (converts to ASCII character)
            if len(seq) == 0: #reached EOF
                break
            sequences.append(seq)
    return sequences