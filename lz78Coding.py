#!/usr/bin/env python
#This file includes functions for the LZ78 coding of the genome sequence -> the reasoning behind the compression.

#The dictionary has to be built and synchronized at both encoder and decoder.
#Inputs are coded with a double (i, c) where:
#   i is an index of the dictionary entry containing the longest match;
#   c is the code of the next input character following the matched sequence.

#Encoding:
#Output = (0, C(x))
#   =>No match has been found in dictionary
#   =>Enter symbol x in dictionary
#Output = (i, C(y))
#   =>i-th dictionary entry has been matched
#   =>Enter i-th dictionary pattern + y in dictionary

#This function encodes the sequence using LZ78 dictionary-based coding/compression.
def encode(data):
    dictionary = {0: ''}
    index = 0
    
    #build and synchronise dictionary at encoder
    dynamicDictionary = lambda dictionary,                                 \
                        key: dictionary.get(key)                           \
                        or dictionary.__setitem__(key, len(dictionary))    \
                        or 0
    
    #output encoding of (i, c) list
    return [token                                                       
            for code in data                                  
                for token in [(index, code)]
                    for index in [dynamicDictionary(dictionary, token)]  
                        if not index                                    
           ] + [(index, '')]

#This function decodes the sequence using LZ78 dictionary-based coding/compression.
def decode(data):
    dictionary = {0: ''}
    j =  ''.join
    
    #build and synchronise dictionary at decoder
    dynamicDictionary = lambda dictionary,                                     \
                        value: dictionary.__setitem__(len(dictionary), value)  \
                        or value
    
    #output decoded seqeunce
    return j([dynamicDictionary(dictionary, dictionary[index] + code)   
                for (index, code) in data])

#print encode('ACGTACGTACGTACGT')
