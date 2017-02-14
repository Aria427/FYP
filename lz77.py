#!/usr/bin/env python
#This file includes functions for the LZ77 compression of the genome sequence.

#This function finds the lookback index for the suffix in the given window.
def lookback(window, suffix):
    return len(window) - (window + suffix).find(suffix)

#This function updates the window by adding x and keeping windowSize most recent symbols.
def grow(window, x, windowSize):
    window += x
    if len(window) > windowSize:
        window = window[-windowSize:] #keep windowSize last symbols
    return window

#This function encodes the sequence using L77 dictionary-based coding/compression.
def encode(data, windowSize):
    encoding = [] #list of codewords
    window = ''   #buffer of last symbols
    suffix = ''   #suffix of data to be coded
    index = 0     #data index
    
    while index < len(data):
        x = data[index]
        #print ('index=' + str(index) + ' | window=' + window + ' | suffix=' + suffix)
        #print ('READ: x=' + x)
       
        #if suffix + x could not be described by previously seen symbols
        if (window + suffix).find(suffix + x) < 0:
           if suffix == '': #no previously seen substring
               codeWord = (0, 0, x) #code x 
               encoding.append(codeWord)
               window = grow(window, x, windowSize) #add x to window
           else:
               i = lookback(window, suffix) #find #symbols back that suffix starts
               #step-size of 1 is hard-coded
               codeWord = (1, i, suffix[0]) #code #symbols
               encoding.append(codeWord)
               window = grow(window, suffix, windowSize) #add suffix to window

               suffix = '' #reset suffix
               index -= 1  #push back last symbol
        else:
           suffix += x #append lasy symbol to search suffix
       
        index += 1 #increment data index

    if suffix != '': #data finished => flush buffer
        i = lookback(window, suffix) 
        codeWord = (1, i, suffix[0])
        encoding.append(codeWord)
        
    return encoding
 
#This function decodes the sequence using L77 dictionary-based coding/compression.
def decode(data, windowSize):  
    decoding = []
    
    

    return decoding
    
#print encode('ACGTACGTACCCGGTT', 4)