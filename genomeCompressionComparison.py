#!/usr/bin/env python
#This file include the function for the comparison of the various genome compression methods implemented. 

import matplotlib.pyplot as plt
from collections import OrderedDict

#This function generates a histogram of the different human genome sizes, before and after compression.
def compressionComparison(output):
    unzippedSize = 3040 #3273481150
    gzipSize = 938 #983659424
    lzwSize = 988 #1036185345
    intSize = 930 #975782472  #fileParsing.parseGenomeInt() 
    longSize = 1810 #1951564944  #fileParsing.parseGenomeLong()
    
    dictionary = {u'Unzipped': unzippedSize, 
                  u'GZIP Compression': gzipSize,  
                  u'LZW Compression': lzwSize,
                  u'Int Compression': intSize, 
                  u'Long Compression': longSize} #key:value
    
    #sort dictionary sizes in descending order
    dictionary = OrderedDict(sorted(dictionary.items(), key=lambda x: x[1], reverse=True))
    
    plt.bar(range(len(dictionary)), dictionary.values(), align='center', color='purple')
    plt.title('Bar Chart of Human Genome Compression Sizes')
    plt.xticks(range(len(dictionary)), dictionary.keys(), rotation=45)
    plt.ylabel('Size / MB')
    plt.tight_layout() #ensures whole plot is saved to file
    plt.savefig(output)
    plt.show()