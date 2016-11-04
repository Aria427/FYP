#!/usr/bin/env python
#This file includes functions for Huffman coding to compress the DNA.

import struct

codes = {}

#This function calculates the frequency of each encountered base (A, C, G, T, N):
def frequency(sequence):
    freq = {} #dictionary
    for base in sequence:
        freq[base] = freq.get(base,0) + 1
    return freq
 
#This function constructs a list of tuples to input in the Huffman algorithm:
def sortFrequency(freq):
    bases = freq.keys()
    tuples = []
    for base in bases:
        tuples.append((freq[base], base))
    tuples.sort() #least frequency to the left; i.e. ascending order
    return tuples   
  
#This function constructs the Huffman coding tree:
def constructTree(tuples):
    while len(tuples) > 1:
        leastBases = tuple(tuples[0:2]) #combine the 2 bases with least frequencies                 
        restBases = tuples[2:]          #rest of bases
        branchFreq = leastBases[0][0] + leastBases[1][0] #branch frequency = sum of child frequencies
        tuples = restBases + [(branchFreq, leastBases)] #add branch to list end
        tuples.sort() #sort branch accordingly
    return tuples[0] #tree inside list

#This function cuts off the frequencies from the tree to nest the representation:
def cutTree(tree) :
    node = tree[1] #ignore root frequency                                  
    if type(node) == type(''): #leaf node 
        return node
    else: 
        left = cutTree(node[0]) #cut left subtree
        right = cutTree(node[1]) #cut right subtree
        return (left, right) #recombine

#This function assigns the traversed tree codes to the bases:
def assignCodes(node, path=''):
    global codes
    if type(node) == type(''): #leaf node
        codes[node] = path #set code to base
    else:                              
        assignCodes(node[0], path + '0') #branch to left
        assignCodes(node[1], path + '1') #branch to right

#This function generates the binary codes from the above functions:
def codeGeneration(sequence):
    freq = frequency(sequence)
    tuples = sortFrequency(freq)
    tree = constructTree(tuples)
    nestedTree = cutTree(tree)
    assignCodes(nestedTree)
    return nestedTree, codes

#This function encodes the DNA into binary codes and saves said encoding to file:       
def encode(sequence, output):
    global codes
    encoding = ''
    outFile = open(output, 'wb')
    for base in sequence: 
        encoding += codes[base]
        outFile.write(struct.pack('=s', codes[base]))
        #outFile.write(codes[base])
    outFile.close()
    return encoding
 
#This function decodes the binary back into its original DNA string:
def decode(tree, binary):
    decoding = ''
    node = tree
    for bit in binary:
        if bit == '0': #left branch 
            node = node[0]
        elif bit == '1': #right branch 
            node = node[1]
        else:
            raise TypeError('Value other than 0 or 1 was found in binary code.')
        if type(node) == type(''): #leaf node
            decoding += node #found base => add to decoding
            node = tree #loop for next base
    return decoding 
    