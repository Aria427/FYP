import plotly.tools as tls
tls.set_credentials_file(username='Aria427', api_key='rois9sg5we')

import collections
import pdb
from path import path
import numpy
import matplotlib.pyplot as plt
import plotly.plotly as plty
import plotly.graph_objs as pltg
from PIL import Image, ImageDraw
import sys

#To read the genome:
def readGenome(filename):
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

#The following is a naive algorithm for exact matching where all occurrences are recorded:
def naiveExact(pattern, text):
    matchOffsets = []
    #loop through every position from where P could start without running past the end of T 
    for i in range(len(text) - len(pattern) + 1): #loop over all possible alignments of P in T from left to right	 
        match = True
        for j in range(len(pattern)): #loop over characters in P from left to right
            #i'th alignment and j'th character            
            if text[i+j] != pattern[j]: #compare characters of P with T
                match = False #mismatch; reject alignment
                break
        if match:
            matchOffsets.append(i)  #all characters matched; record
    return matchOffsets

#Hamming distance = minimum no of substitutions required to change one string into another.
#The following is also a naive algorithm but for approximate matching using the Hamming distance:
def naiveApproxHamming(pattern, text, maxHammingDist=1):
    matchOffsets = []
    for i in range(len(text) - len(pattern) + 1): 
        mismatches = 0
        for j in range(len(pattern)):
            if text[i+j] != pattern[j]:
                mismatches += 1     #mismatch
                if mismatches > maxHammingDist:
                    break           #exceeded maximum distance
        if mismatches <= maxHammingDist: #approximate match
            matchOffsets.append(i)
    return matchOffsets      

#The genome is double stranded and so the reads can come from one strand or the other.    
#To match both the read and the reverse complement of the read to the genome: 
def reverseComplement(read):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'} #each base is associated with its complementary base
    reverseRead = ''
    for base in read:
        reverseRead = complement[base] + reverseRead #complement added to beginning in order to reverse the read from end to start
    return reverseRead

#To align the reads against the genome to see how many match:
def align(sequence, genome):
    readsMatched = 0
    readsCount = 0
    readsOffsets = []
    for read in reads:
        read = read[:25] #prefix of read as all 100 bases have a smaller chance of matching
        matches = naiveApproxHamming(read, genome) #check if read matches in forward direction of genome
        matches.extend(naiveApproxHamming(reverseComplement(read), genome)) #add results of any matches in reverse complement of genome
        readsCount += 1
        if len(matches) > 0: #match - read aligned in at least one place
            readsMatched += 1
        readsOffsets.append(matches)
    return readsMatched, readsCount, readsOffsets

#To create a data visualisation of the matched reads against the genome:
def visualisation(readsOffsets, outputFile):
    readsOffsets.sort()                     #sort list
    readsOffsets = sum(readsOffsets, [])    #flatten list
    offsetsCount = collections.Counter(readsOffsets) #record count of each offset => length of match
    #print(len(offsetsCount)) #181
    file = open(outputFile, 'w')
    #arranges offsets with corresponding match length in ascending order
    for i in range(len(genome)):
        if offsetsCount[i] != 0:
            #f.write(i, ":", offsetsCount[i])
            file.write('-' * offsetsCount[i])
        elif (offsetsCount[i] == 0) & (offsetsCount[i+1] != 0):
            file.write('\n')
            file.write(' ' * (i+1)) #indentation  
        elif (offsetsCount[i] == 0) & (offsetsCount[i+1] == 0):
            file.write(' ')
    file.close()
    return outputFile
    
#pdb.set_trace()
#genome = readGenome(path('Data\HumanGenome.fa.gz').abspath())
#reads = readSequence(path('Data\HumanSequencingReads.tsv.bz2').abspath())
genome = readGenome(path('Data\PhixGenome.fa').abspath())
reads = readSequence(path('Data\PhiXSequencingReads1000.fastq').abspath())

#print("Length of the genome: ", len(genome))
#print("Frequency of each base in the genome: ", collections.Counter(genome))

#print("Length of the reads: ", len(reads))
#readsFreq = collections.Counter()
#for read in reads:
    #readsFreq.update(read)
#print("Frequency of each base in the reads: ", readsFreq)

matches, count, offsets = align(reads, genome)
print(matches, "/", count, " reads matched the genome")
#The result is not 100% but this is to be expected due to sequencing errors. 
    
file = path('Output Test Files\DataVisualisationTest.txt').abspath()
file = visualisation(offsets, file)

img = Image.new('RGBA', (5000, 5000), (255, 255, 255, 0)) 
draw = ImageDraw.Draw(img) 
offsets.sort()                     #sort list
offsets = sum(offsets, [])    #flatten list
offsetsCount = collections.Counter(offsets)
for i in range(len(genome)):
    if offsetsCount[i] != 0:
        draw.line(((i,1000), (offsetsCount[i],1000)), fill=128, width=3)
        #file.write('-' * offsetsCount[i])
    #elif (offsetsCount[i] == 0) & (offsetsCount[i+1] != 0):
        #file.write('\n')
        #file.write(' ' * (i+1)) #indentation  
    #elif (offsetsCount[i] == 0) & (offsetsCount[i+1] == 0):
        #file.write(' ')
draw.line((1000,2000, 1500,3000), fill=128, width=2)
img.show()

"""
x = []
y = []

for i in range(len(genome)):
    x.extend([i])
    if offsetsCount[i] != 0:
        y.extend([i+offsetsCount[i]])
    elif (offsetsCount[i] == 0) & (offsetsCount[i+1] != 0):
        y.extend([None]) #None
    elif (offsetsCount[i] == 0) & (offsetsCount[i+1] == 0):
        y.extend([None])

trace = pltg.Scatter(x, y)

fig = dict(data=[trace])
plty.iplot(fig, filename='simple-connectgaps')"""

#fig.savefig('yourfilename.png')
    
