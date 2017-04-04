#!/usr/bin/env python
#This file contains the map step for the aligment MapReduce implementation.

import sys
import gzip

GAP = 1
MISMATCH = 1
MATCH = 0

#This function reads the genome in chunks (groups of bytes) from the file stored in S3.
def readGenomeChunks(s3File, bytesNum=100):
    with gzip.open(s3File, 'r') as f:
        f.readline()
        for chunk in iter(lambda: f.read(bytesNum), ''):
            data = chunk.rstrip().upper().replace('N', '').replace('\n', '').replace(' ', '')
            yield data

def readGenome(s3File):
    genome = ''
    with gzip.open(s3File, 'r') as f: 
        for line in f:
            if line and line[0] != '>': #ignore header line with genome information
                l = line.rstrip().upper().replace('N', '').replace(' ', '')
                genome += l
        return genome
                  
#This function reads the sequencing reads as input to the mapper.        
def readInputReads(file):
    flag, sequence, quality = '', '', ''
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

#This function creates a score matrix of trial alignments of the two sequences.
def scoreMatrix(seq1, seq2):
    rows = len(seq1) + 1 #extra row for gap (-) => +1
    cols = len(seq2) + 1 #extra column for gap (-) => +1
    
    #2D matrix/graph of scores based on trial alignments of different base pairs
    scoreMatrix = [[0 for col in range(cols)] for row in range(rows)]

    maxScore = 0
    maxPos = None #(x, y) or (row, column) position of highest score in matrix
    
    #fill in scoring matrix
    for i in range(1, rows):
        for j in range(1, cols):
            score = calculateScore(seq1, seq2, scoreMatrix, i, j)
            if score > maxScore:
                #path with the highest cummulative score is the best alignment
                maxScore = score
                maxPos = (i, j) 
            scoreMatrix[i][j] = score

    assert maxPos is not None, 'The (x, y) position with the highest score was not found.'

    return scoreMatrix, maxPos

#This function calculates the score of a given position in the score matrix.
def calculateScore(seq1, seq2, scoreMatrix, x, y):
    similarity = MATCH if seq1[x - 1] == seq2[y - 1] else MISMATCH

    #score is based on upper-left(diagonal), upper and left neighbours
    diagonalScore = scoreMatrix[x - 1][y - 1] + similarity
    upScore = scoreMatrix[x - 1][y] + GAP
    leftScore = scoreMatrix[x][y - 1] + GAP

    return max(0, diagonalScore, upScore, leftScore)    
    
#This function determines the next move to be taken by the trace path.
def nextMove(scoreMatrix, x, y):
    diagonal = scoreMatrix[x - 1][y - 1]
    up   = scoreMatrix[x - 1][y]
    left = scoreMatrix[x][y - 1]
    
    
    if diagonal >= up and diagonal >= left: #tie goes to diagonal move
        #1 -> diagonal move, 0 -> end
        return 1 if diagonal != 0 else 0
    
    elif up > diagonal and up >= left: #tie goes to up move
        #2 -> up move, 0 -> end
        return 2 if up != 0 else 0
    
    elif left > diagonal and left > up: #tie goes to left move
        #3 -> left move, 0 -> end
        return 3 if left != 0 else 0         
    
#This function traces a path from the bottom-right to the top-left corner of the score matrix;
#It returns the optimal path.
def traceback(scoreMatrix, startPosition, seq1, seq2):
    #A move is determined by the score of the 3 adjacent squares.
    #possible moves
    diagonal = 1    #match/mismatch
    up = 2          #gap in sequence 1
    left = 3        #gap in sequence 2
    end = 0
        
    alignedSeq1 = []
    alignedSeq2 = []
    x, y = startPosition
    
    move = nextMove(scoreMatrix, x, y)
    
    while move != end:
        if move == diagonal:
            alignedSeq1.append(seq1[x - 1])
            alignedSeq2.append(seq2[y - 1])
            x -= 1
            y -= 1
        elif move == up:
            alignedSeq1.append(seq1[x - 1])
            alignedSeq2.append('-')
            x -= 1
        else: #if move == left
            alignedSeq1.append('-')
            alignedSeq2.append(seq2[y - 1])
            y -= 1

        move = nextMove(scoreMatrix, x, y)

    alignedSeq1.append(seq1[x - 1])
    alignedSeq2.append(seq2[y - 1])

    return ''.join(reversed(alignedSeq1)), ''.join(reversed(alignedSeq2))   

#This function generates a special string displaying identities (|), gaps (-), and mismatches (:).
def smithWatermanApproximateString(alignedSeq1, alignedSeq2, maxMismatches):
    identities, gaps, mismatches = 0, 0, 0
    alignmentString = []
    matchOffsets = []

    for i, (base1, base2) in enumerate(zip(alignedSeq1, alignedSeq2)):
        if base1 == base2:
            alignmentString.append('|')
            identities += 1
            if mismatches <= maxMismatches:
                matchOffsets.append(i)
        elif '-' in (base1, base2):
            alignmentString.append(' ')
            gaps += 1
        else:
            alignmentString.append(':')
            mismatches += 1

    return matchOffsets, ''.join(alignmentString), identities, gaps, mismatches

#This function implements the Smith Waterman approximate matching algorithm.
def smithWatermanApproximate(text, pattern, maxMismatches):
    sMatrix, startPos = scoreMatrix(text, pattern) #intialise score matrix

    #optimal path through the score matrix = optimal local sequence alignment
    refAligned, readAligned = traceback(sMatrix, startPos, text, pattern)
    assert len(refAligned) == len(readAligned), 'aligned strings are not the same size'
    
    matchOffsets, alignmentString, identities, gaps, mismatches =  \
            smithWatermanApproximateString(text, pattern, maxMismatches)
    
    return matchOffsets               
                
#This function finds the reverse complement of a sequencing read.   
def reverseComplement(read):
    complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'} #each base is associated with its complementary base
    reverseRead = ''
    for base in read:
        reverseRead = complement[base] + reverseRead #complement added to beginning in order to reverse the read from end to start
    return reverseRead

#This function takes the quality value Q (rounded integer) and converts it to its respective character. 
def QtoPhred33(Q):
    return chr(Q + 33) #converts integer to character according to ASCII table

#This function takes the Phred-33 encoded character and converts it back to Q.  
def phred33ToQ(qual):
    return ord(qual) - 33 #converts character to integer according to ASCII table

#This function aligns the reads to the genome using the Smith Waterman local alignment algorithm.
def alignSmithWaterman(read, quality, genome):
    readQualityDictionary = {} #key:read, value:list of quality integers
    
    #PhiX read line=100 & genome line=70, Human read line=60 & genome line=50
    read1 = read[:50] #first 50 characters 
    read2 = read[50:] #last 50 characters

    #maximum number of mismatches = 2
    reverseRead1 = reverseComplement(read1)
    matchOffset = smithWatermanApproximate(genome, read1, 2) #check if read matches in forward direction of genome
    matchOffset.extend(smithWatermanApproximate(genome, reverseRead1, 2)) #add results of any matches in reverse complement of genome
    reverseRead2 = reverseComplement(read2)
    matchOffset = smithWatermanApproximate(genome, read2, 2) #check if read matches in forward direction of genome
    matchOffset.extend(smithWatermanApproximate(genome, reverseRead2, 2))       

    if len(list(matchOffset)) > 0: #match - read aligned in at least one place
        qualityQ = []
        for q in quality:
            qualityQ.append(phred33ToQ(q))
        readQualityDictionary[read] = qualityQ

    return matchOffset, readQualityDictionary 
    
def main():
    #hard-coded reference genome stored in S3 via Amazon EMR
    #genomeFile = 's3://fyp-input-gen/HumanGenome_200000.fa.gz' #Frankfurt region doesn't work, Ireland does
    genomeFile = './human' #-cacheArchive s3://fyp-input/HumanGenome.fa.gz#human
    #g = readGenome(genomeFile)   
    readSeq = readInputReads(sys.stdin) #Human reads=28,094,847
    
    for read, quality in readSeq: 
        #Human genome=64,185,939 lines -> 3,273,481,150 bytes
        genome = readGenomeChunks(genomeFile, 250000) #250,000 bytes = 0.23842MB
        overlap = '' #size of read-1
        filesOffset = 0 #file is split in chunks so offset needs to change according to chunk
        
        for g in genome:
            g = overlap + g
            offset, rqDict = alignSmithWaterman(read, quality, g)
         
            #write results to STDOUT (standard output)
            for o in offset: #to remove empty list and '[' ']' characters
                #tab-delimited, key:offset of match with reads, value:<default count of 1, genome subsequence matched, read matched, corresponding quality> 
                print '%s\t%s\t%s\t%s\t%s' % (o+filesOffset, 1, g[o:o+len(read)], read, quality) 
                #The output here will be the input for the reduce step  
            
            overlap = g[-99:] #100-1 for PhiX, 60-1 for Human read => -1 as last 60 have already been read
            filesOffset += (len(g)-100) #store offset according to overlap as file is read in chunks  
 
if __name__ == '__main__':
    main()            
