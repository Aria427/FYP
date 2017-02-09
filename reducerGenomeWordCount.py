#!/usr/bin/env python

import sys
from itertools import groupby
from operator import itemgetter

#This function reads in the output from the mapper using a generator.
def readMapperOutput(file, separator='\t'):
    for line in file:
        yield line.rstrip().split(separator, 1)

data = readMapperOutput(sys.stdin, separator='\t') #input comes from STDIN (standard input)
  
#groupby() groups word-count pairs by word and generates iterator 
# which returns consecutive keys and their respective group:
#   currentWord -> word (key) 
#   group -> iterator yielding all [currentWord, count] items
for currentWord, group in groupby(data, itemgetter(0)):
    try:
        totalCount = sum(int(count) for currentWord, count in group)
        print "%s\t%d" % (currentWord, totalCount) #write result to STDOUT
    except ValueError:
        pass #when count not a number

"""        
currentWord = None
currentCount = 0
word = None

for line in sys.stdin: #input comes from STDIN (standard input)
    line = line.strip() #remove whitespace
    word, count = line.split('\t', 1) #parse input from map step

    try:
        count = int(count) #convert count (currently a string) to int
    except ValueError:
        continue #when count not a number

    #This IF-switch only works because Hadoop sorts map output 
    # by key (word) before it is passed to the reducer.
    if currentWord == word:
        currentCount += count
    else:
        if currentWord:
            print '%s\t%s' % (currentWord, currentCount) #write result to STDOUT
        current_count = count
        current_word = word

if currentWord == word:
    print '%s\t%s' % (currentWord, currentCount) #output last word
"""
