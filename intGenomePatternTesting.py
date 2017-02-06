#!/usr/bin/env python
#This file contains some testing functions with regards to finding repetitive patterns in the compressed integer genome.

import numpy as np
from itertools import groupby, combinations
from operator import itemgetter
from collections import Counter
import matplotlib.pyplot as plt
import matplotlib.ticker as plticker

#This function returns the frequency of each integer found in the genome.
#Integer 0 occurs most with frequency 319594 in human genome (total integers = 83203868).
def countIntegers(genome):
    uniqueKeys, counts = np.unique(genome, return_counts=True) 
    #frequencyArray = np.asarray((uniqueKeys, counts)).T
    return uniqueKeys, counts

def createHist(genome):
    hist, bins = np.histogram(genome, density=True)
    width = 0.7*(bins[1]-bins[0])
    centre = (bins[:-1]+bins[1:])/2
    plt.bar(centre, hist, align='center', width=width)
    plt.show() 
    
#This function identifies each integer and pair/triple/... of integers found in the genome.
def identifyConsecutiveIntegers(genome):
    for k, g in groupby(enumerate(genome), lambda (i, x): i-x):
        group = map(itemgetter(1), g) 
        yield group
    
#This function returns the frequency of each pair of consecutive integers found in the genome.
#Pair (0,1) occurs most with frequency 380 in human genome (total pairs = 178)
def countIntegerPairs(genome):
    pairs = []
    counter  = Counter()
    group = identifyConsecutiveIntegers(genome)
    for g in group:
        if len(g) > 1 & len(g) < 3: #if pair of integers is observed
            pairs.append(g)
    for p in pairs:
        p.sort()
        for c in combinations(p, 2):
            counter[c] += 1
    return counter

#This function returns the frequency of each triple of consecutive integers found in the genome:
#No triples exist in human genome so this function does not apply.
def countIntegerTriples(genome):
    triples = []
    counter  = Counter()
    group = identifyConsecutiveIntegers(genome)
    for g in group:
        if len(g) > 2 & len(g) < 4: #if triple of integers is observed
            triples.append(g)
    for t in triples:
        t.sort()
        for c in combinations(t, 2):
            counter[c] += 1
    return triples
 
#This function creates a histogram displaying the frequencies of the relative integers:
def createIntegerHistogram(counts):
    fig, ax = plt.subplots()
    ax.hist(counts[0:35])
    loc = plticker.MultipleLocator(base=1.0) #this locator puts ticks at regular intervals
    ax.xaxis.set_major_locator(loc)
    plt.show()

#This function creates a histogram displaying the frequencies of the relative integer pairs:
def createPairHistogram(counter):
    for pair, count in counter.items():
        print(pair, count)
    fig, ax = plt.subplots()
    ax.hist(count)
    loc = plticker.MultipleLocator(base=1.0)
    ax.xaxis.set_major_locator(loc)
    plt.show()  
