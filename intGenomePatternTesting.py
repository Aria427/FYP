#!/usr/bin/env python
#This file contains some testing functions with regards to finding repetitive patterns in the compressed integer genome.

import numpy as np
from itertools import groupby
from operator import itemgetter
import matplotlib.pyplot as plt
import matplotlib.ticker as plticker

#This function returns the frequency of each integer found in the genome:
def countIntegers(genome):
    unique, counts = np.unique(genome, return_counts=True)
    #frequencyArray = np.asarray((unique, counts)).T
    return unique, counts

#This function identifies each integer and pair/triple/... of integers found in the genome:
def identifyConsecutiveIntegers(genome):
    for k, g in groupby(enumerate(genome), lambda (i, x): i-x):
        group = map(itemgetter(1), g) 
        yield group
    
#This function returns the frequency of each group of consecutive integers found in the genome:
def countConsecutiveIntegers(genome):
    pairs = []
    group = identifyConsecutiveIntegers(genome)
    for g in group:
        if len(g) > 1: #if pair of integers is observed
            pairs.append(g)
    unique, counts = np.unique(pairs, return_counts=True)
    #frequencyArray = np.asarray((unique, counts)).T
    return unique, counts  
 
#This function creates a histogram displaying the frequencies of the relative integers:
def createHistogram(counts):
    fig, ax = plt.subplots()
    ax.hist(counts[0:35])
    loc = plticker.MultipleLocator(base=1.0) #this locator puts ticks at regular intervals
    ax.xaxis.set_major_locator(loc)
    plt.show()