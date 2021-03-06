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
                  u'FM Index': 0,
                  u'Burrows Wheeler': 0}
             
    #sort in alphabetical order
    dictionary = OrderedDict(sorted(dictionary.items(), key=lambda x: x[0]))              
    
    aligners = ['Boyer Moore', 'Burrows Wheeler', 'FM Index', 'Hamming', 'Kmer Index']
    n = len(aligners)
    index = np.arange(n)
    
    cat = [27.38, 26.19, 0.06, 31.56, 2.40]
    mr1 = [97.48, 65.45, 0.16, 34.30, 2.46]
    mr5 = [33.50, 33.58, 0.08, 17.47, 1.22]
    preProc = [0, 7.03, 3.18, 0, 0.15]
    
    catBar = plt.bar(index-0.2, cat, width=0.2, color='blue', align='center', bottom=preProc)
    mr1Bar = plt.bar(index, mr1, width=0.2, color='red', align='center', bottom=preProc)
    mr5Bar = plt.bar(index+0.2, mr5, width=0.2, color='purple', align='center', bottom=preProc)
    preBar1 = plt.bar(index-0.2, preProc, width=0.2, color='green', align='center')
    preBar2 = plt.bar(index, preProc, width=0.2, color='green', align='center')
    preBar3 = plt.bar(index+0.2, preProc, width=0.2, color='green', align='center')
    
    plt.title('Bar Chart of Map-Reduce Aligners Run Times')
    plt.xticks(range(len(dictionary)), dictionary.keys(), rotation=45)
    plt.xlabel('Map-Reduce Aligners')
    plt.ylabel('Time / Minutes')
    plt.legend([preBar1, catBar, mr1Bar, mr5Bar], ['Preprocessing', 'cat', 'Hadoop1', 'Hadoop5'])
    plt.tight_layout() #ensures whole plot is saved to file
    plt.savefig(output)
    plt.show()
    
outputFile = path('Documentation\genTimeAlign.png').abspath()
compressionComparison(outputFile)

