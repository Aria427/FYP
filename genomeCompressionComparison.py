#!/usr/bin/env python
#This file include the function for the comparison of the various genome compression methods implemented. 

import matplotlib as plt

#This function generates a histogram of the different human genome sizes, before and after compression.
def compressionComparison():
    unzippedSize = 3273481150
    gzippedSize = 983,659,424
    intSize = 0  #fileParsing.parseGenomeInt() 
    longSize = 0 #fileParsing.parseGenomeLong()
    
    dictionary = {u'Unzipped': unzippedSize, u'G-zipped': gzippedSize,  \
                  u'Int Compression': intSize, u'Long Compression': longSize} #key:value

    plt.bar(range(len(dictionary)), dictionary.values(), align='center', color='purple')
    plt.title('Bar Chart of Human Genome Compression Sizes')
    plt.xticks(range(len(dictionary)), dictionary.keys())
    plt.ylabel('Size')
    plt.show()