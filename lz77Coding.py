#!/usr/bin/env python
#This file includes functions for the LZ77 encoding of the genome sequence -> the reasoning behind the compression.

import gzip

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
    encoding = [] #codewords
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
   
#This function encodes the sequence found in the file in chunks.
def encodeFile(sequenceFile, encodedFile):
    out = open(encodedFile, 'w')
    with gzip.open(sequenceFile) as f:
        #read genome into line-by-line generator     
        subseqs = (line.upper().replace('\n', '').replace('N', '') 
                    for line in f if line and line[0] != '>') #ignore header line
                
        encoding = ''
        lineCount = 0
        for s in subseqs:
            try:
                #encoding every two lines is more accurate than every one line
                #window will lose some characters if every one line is considered
                s = s + next(subseqs) 
            except StopIteration:
                pass
            if lineCount == 0: #normal encoding
                encoding = encoding + str(encode(s, 4)).strip('[]') + ', ' #format output
            else: 
                #ignore first codeword as it produces non-sensical output;
                #this is due to slight reading difficulties as every 2 lines are taken
                encoding = encoding + str(encode(s, 4)[1:]).strip('[]') + ', '
            lineCount += 1
        out.write(encoding)
    out.close()    
    return encoding
    
#This function encodes the sequence found in the file as a whole.
def encodeWholeFile(sequenceFile, encodedFile):
    out = open(encodedFile, 'w')
    with open(sequenceFile) as f:
        data = ''
        for line in f:
            if line and line[0] != '>':
                data = data + line.upper().replace('\n', '').replace('N', '')
        encoding = encode(data, 4)
        out.write(str(encoding)) 
    out.close()
    return encoding      
    