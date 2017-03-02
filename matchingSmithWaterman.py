#!/usr/bin/env python
#This file includes functions for the Smith-Waterman algorithm.

#scores taken from Wikipedia
MATCH = 2
MISMATCH = -1
GAP = -1

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
def displayAlignment(alignedSeq1, alignedSeq2):
    identities, gaps, mismatches = 0, 0, 0
    alignmentString = []

    for base1, base2 in zip(alignedSeq1, alignedSeq2):
        if base1 == base2:
            alignmentString.append('|')
            identities += 1
        elif '-' in (base1, base2):
            alignmentString.append(' ')
            gaps += 1
        else:
            alignmentString.append(':')
            mismatches += 1

    return ''.join(alignmentString), identities, gaps, mismatches

#This function dispays the score matrix generated.
def displayMatrix(scoreMatrix):
    for row in scoreMatrix:
        for col in row:
            print('{0:>4}'.format(col))
        print()

#Smith-Waterman alignment
seq1 = 'ATAGACGACATACAGACAGCATACAGACAGCATACAGA'
seq2 = 'TTTAGCATGCGCATATCAGCAATACAGACAGATACG'

scoreMatrix, startPos = scoreMatrix(seq1, seq2) #intialise score matrix

#optimal path through the score matrix = optimal local sequence alignment
seq1Aligned, seq2Aligned = traceback(scoreMatrix, startPos, seq1, seq2)
#assert len(seq1Aligned) == len(seq2Aligned), 'aligned strings are not the same size'

#pretty print results
alignmentString, identities, gaps, mismatches = displayAlignment(seq1Aligned, seq2Aligned)

print 'Identities = %d/%d' % (identities, len(seq1Aligned))
print 'Gaps = %d/%d' % (gaps, len(seq1Aligned))

for i in range(0, len(seq1Aligned), 60):
    seq1Slice = seq1Aligned[i:i+60]
    print('Query  {0:<4}  {1}  {2:<4}'.format(i + 1, seq1Slice, i + len(seq1Slice)))
    print('             {0}'.format(alignmentString[i:i+60]))
    
    seq2Slice = seq2Aligned[i:i+60]
    print('Sbjct  {0:<4}  {1}  {2:<4}'.format(i + 1, seq2Slice, i + len(seq2Slice)))
