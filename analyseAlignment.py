#!/usr/bin/env python
#This file includes functions for calculating the time it takes to align versus the matches produced.

import alignmentString

from time import time
from matplotlib import pyplot as plt

def plotTimeVsMatches(pattern, text, output):
    timeHamming, matchesHamming = analyseHamming(pattern, text)
    timeEdit, matchesEdit = analyseEdit(pattern, text)
    timeKmer, matchesKmer = analyseKmer(pattern, text)
    timeFM, matchesFM = analyseFM(pattern, text)
    hamming, = plt.plot(timeHamming, matchesHamming, 'o', label='Hamming')
    edit, = plt.plot(timeEdit, matchesEdit, 'o', label='Edit')
    kmer, = plt.plot(timeKmer, matchesKmer, 'o', label='Kmer')
    fm, = plt.plot(timeFM, matchesFM, 'o', label='FM')
    plt.title('Time Duration vs Matches Generated')
    plt.xlabel('Time/seconds')
    plt.ylabel('Matches')
    plt.legend(handles=[hamming, edit, kmer, fm], loc=4) #4=lower-right
    plt.ylim([0, 110])
    plt.savefig(output)
    plt.show()
    
def analyseHamming(pattern, text):
    start = time()
    readsMatched, readsCount, _ = alignmentString.alignHamming(pattern, text)
    end = time()
    timeLength = (end - start)
    print 'Approximate Hamming lasted %.2f seconds to produce %d/%d matches.' % (timeLength, readsMatched, readsCount)
    return timeLength, readsMatched
    
def analyseEdit(pattern, text):
    start = time()
    readsMatched, readsCount, _ = alignmentString.alignEdit(pattern, text)
    end = time()
    timeLength = (end - start)
    print 'Approximate Edit lasted %.2f seconds to produce %d/%d matches.' % (timeLength, readsMatched, readsCount)
    return timeLength, readsMatched    
    
def analyseKmer(pattern, text):
    start = time()
    readsMatched, readsCount, _ = alignmentString.alignKmer(pattern, text)
    end = time()
    timeLength = (end - start)
    print 'Kmer indexing lasted %.2f seconds to produce %d/%d matches.' % (timeLength, readsMatched, readsCount)
    return timeLength, readsMatched
    
def analyseFM(pattern, text):
    start = time()
    readsMatched, readsCount, _ = alignmentString.alignFM(pattern, text)
    end = time()
    timeLength = (end - start)
    print 'FM indexing lasted %.2f seconds to produce %d/%d matches.' % (timeLength, readsMatched, readsCount)
    return timeLength, readsMatched
    