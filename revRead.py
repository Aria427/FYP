#!/usr/bin/env python
#This file contains the optimisation for the reads by storing the reverse complement of each.

from path import path

#This function reads and outputs the sequencing reads one by one.      
def readInputReads(readsFile):
    flag, sequence, quality = '', '', ''
    with open(readsFile, 'r') as file:
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

#This function finds the reverse complement of a sequencing read.   
def reverseComplement(read):
    complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'} #each base is associated with its complementary base
    reverseRead = ''
    for base in read:
        reverseRead = complement[base] + reverseRead #complement added to beginning in order to reverse the read from end to start
    return reverseRead

"""
readsFile = path('Data\HumanSequencingReads.tsv').abspath() 
optimisedFile = path('Output Data\HumanReads.tsv').abspath()

readSeq = readInputReads(readsFile) #Human reads=28,094,847

with open(optimisedFile, 'w') as f:    
    for read, quality in readSeq: 
        #ignore flag and N's as not used
        f.write( '%s\t%s\n' % (read, quality) )
        print read
        reverseRead = reverseComplement(read)
        f.write( '%s\t%s\n' % (reverseRead, quality) )
"""       
#This function reads and outputs the sequencing reads one by one after optimisation.      
def readOptimisedReads(optimisedReadsFile):
    sequence, quality = '', ''
    with open(optimisedReadsFile, 'r') as file:
        while True: #runs until EOF
            line = file.readline() 
            if not line: #reached EOF
                break
            
            line = line.split()
            sequence = line[0]
            quality = line[1]

            yield sequence, quality

optimisedFile = path('Output Data\HumanReads.tsv').abspath()            
optReadSeq = readOptimisedReads(optimisedFile)   
 
with open(optimisedFile, 'r') as f:    
    for read, quality in optReadSeq:   
        print read, quality
        
        
    