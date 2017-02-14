#!/usr/bin/env python
#This file includes functions for the LZW compression of the genome sequence.

#More practical adaptation of LZ78.
#Codeword (i,c) becomes (i) -> only the index is transmitted.
#The dictionary must be primed with all symbols of the source alphabet;
#also built and synchronized at both encoder and decoder.

#This function encodes the sequence using LZW dictionary-based coding/compression.
def encode(data):
    #build and synchronise dictionary at encoder
    dictSize = 256
    dictionary = dict((chr(i), i) for i in xrange(dictSize))
 
    encoding = []
    p = "" #pattern p

    for a in data:
        entry = p + a
        if entry in dictionary:
            p = entry #accumulate input symbols in p as long as p is in dictionary
        else: #when p+a is not in dictionary:
            dictionary[entry] = dictSize #new dictionary entry: p+a
            encoding.append(dictionary[p]) #new output: index of p
            dictSize += 1
            p = a #next p is a
 
    if p: #last output occurs without entry of new character
        encoding.append(dictionary[p]) #last output: last index of p
        
    return encoding #dictionary indices of encoded sequence
 
#This function decodes the sequence using LZW dictionary-based coding/compression.
def decode(data):
    #build and synchronise dictionary at decoder
    dictSize = 256
    dictionary = dict((chr(i), chr(i)) for i in xrange(dictSize))
    
    p = result = data.pop(0) #start with p=NULL

    for a in data: #decode input a
        if a in dictionary: #for first decoded input only p+a=a
            entry = dictionary[a]
        elif a == dictSize: #p+a not in dictionary
            entry = p + a
        else:
            raise ValueError('Bad compressed a: %s' % a)
        result.append(entry)

        dictionary[dictSize] = p + entry[0] #add p+a in dictionary
        dictSize += 1

        p = entry #next p=a (last decoded input)
        
    return result