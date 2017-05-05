#!/usr/bin/env python
#This file include the function for the time comparison of the run BWA. 

import matplotlib.pyplot as plt
import numpy as np
from collections import OrderedDict
from path import path

#This function generates a histogram of the different human genome compression times.
def compressionComparison(output):    
    dictionary = {u'1 thread': 0,  
                  u'5 threads': 0 }
             
    #sort in alphabetical order
    dictionary = OrderedDict(sorted(dictionary.items(), key=lambda x: x[0]))              
    
    data = np.array([[0.14, 0.14],
                    [0.005, 0.002],
                    [0.005, 0.005]])
    
    colours = ['green', 'purple', 'red']
    labels = ['Preprocessing', 'Alignment', 'File Conversion']

    index = np.arange(data.shape[1])
    for i in range(data.shape[0]):
        plt.bar(index, data[i], width=0.2, color=colours[i % len(colours)],
                align='center', bottom=np.sum(data[:i], axis = 0),
                label=labels[i])
    
    plt.title('Bar Chart of BWA Run Times')
    plt.xticks(range(len(dictionary)), dictionary.keys(), rotation=45)
    plt.xlabel('Number of Threads')
    plt.ylabel('Time / Minutes')
    plt.legend(loc='best')
    plt.tight_layout() #ensures whole plot is saved to file
    plt.savefig(output)
    plt.show()
    
outputFile = path('Documentation/bwaTimeAlign.png').abspath()
compressionComparison(outputFile)


