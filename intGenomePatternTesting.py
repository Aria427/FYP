#!/usr/bin/env python
#This file contains some testing functions with regards to finding repetitive patterns in the compressed integer genome.

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as plticker

#This function returns the frequency of each integer found in the genome:
def countIntegers(genome):
    unique, counts = np.unique(genome, return_counts=True)
    frequencyArray = np.asarray((unique, counts)).T
    return unique, counts

#This function returns the frequency of each pair of integers found in the genome:
def countIntegerPairs(genome):
    unique, counts = np.unique((genome, genome[1:]), return_counts=True)
    frequencyArray = np.asarray((unique, counts)).T
    return frequencyArray    
    
#This function returns the frequency of each triple of integers found in the genome:
def countIntegerTriples(genome):
    unique, counts = np.unique((genome, genome[1:], genome[2:]), return_counts=True)
    frequencyArray = np.asarray((unique, counts)).T
    return frequencyArray
 
#This function creates a histogram displaying the frequencies of the relative integers:
def createHistogram(counts):
    fig, ax = plt.subplots()
    ax.hist(counts[0:35])
    loc = plticker.MultipleLocator(base=1.0) #this locator puts ticks at regular intervals
    ax.xaxis.set_major_locator(loc)
    plt.show()