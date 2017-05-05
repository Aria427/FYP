#!/usr/bin/env python
#This file include the function for the time comparison of the various genome compression methods implemented. 

import matplotlib.pyplot as plt
from collections import OrderedDict
from path import path

#This function generates a histogram of the different human genome compression times.
def compressionComparison(output):
    gzipTime = 18.23
    lzwTime = 241.44
    intTime = 30.38
    longTime = 29.58
    
    dictionary = {u'GZIP': gzipTime,  
                  u'LZW': lzwTime,
                  u'Int': intTime, 
                  u'Long': longTime} #key:value
    
    #sort dictionary sizes in descending order
    dictionary = OrderedDict(sorted(dictionary.items(), key=lambda x: x[1], reverse=True))
    
    plt.bar(range(len(dictionary)), dictionary.values(), align='center', color='green', width=0.5)
    plt.title('Bar Chart of Human Genome Compression Run Times')
    plt.xticks(range(len(dictionary)), dictionary.keys(), rotation=45)
    plt.xlabel('Compression Methods')
    plt.ylabel('Time / Minutes')
    plt.tight_layout() #ensures whole plot is saved to file
    plt.savefig(output)
    plt.show()
    
outputFile = path('Documentation\genTimeComp.png').abspath()
compressionComparison(outputFile)