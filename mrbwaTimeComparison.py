#!/usr/bin/env python
#This file include the function for the time comparison of the various aligners implemented. 

import matplotlib.pyplot as plt
import numpy as np
from collections import OrderedDict
from path import path

#This function generates a histogram of the different human genome compression times.
def compressionComparison(output):    
    dictionary = {u'Hamming': 0,  
                  u'Boyer Moore': 0,
                  u'Kmer Index': 0, 
                  u'Burrows Wheeler': 0,
                  u'BWA': 0     }
             
    #sort in alphabetical order
    dictionary = OrderedDict(sorted(dictionary.items(), key=lambda x: x[0]))              
    
    aligners = ['Boyer Moore', 'Burrows Wheeler', 'BWA', 'Hamming', 'Kmer Index']
    n = len(aligners)
    index = np.arange(n)
    
    align1 = [0.01, 97.48, 65.45, 34.30, 2.46]
    align5 = [0.007, 33.50, 33.58, 17.47, 1.22]
    preProc = [0.14, 0, 7.03, 0, 0.15]
    
    align1Bar = plt.bar(index, align1, width=0.2, color='red', align='center', bottom=preProc)
    align5Bar = plt.bar(index+0.2, align5, width=0.2, color='purple', align='center', bottom=preProc)
    preBar1 = plt.bar(index, preProc, width=0.2, color='green', align='center')
    preBar5 = plt.bar(index+0.2, preProc, width=0.2, color='green', align='center')
    
    plt.title('Bar Chart of Aligners Run Times')
    plt.xticks(range(len(dictionary)), dictionary.keys(), rotation=45)
    plt.xlabel('Aligners')
    plt.ylabel('Time / Minutes')
    plt.legend([preBar1, align1Bar, align5Bar], ['Preprocessing', '1 Node', '5 Nodes'])
    plt.tight_layout() #ensures whole plot is saved to file
    plt.savefig(output)
    plt.show()
    
outputFile = path('Documentation\mrbwaTimeAlign.png').abspath()
compressionComparison(outputFile)

