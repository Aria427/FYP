#!/usr/bin/env python
#This file contains the file splitting for the humnan sequencing reads.

from path import path
from itertools import izip_longest

readsFile = path('Data\HumanSequencingReads.tsv').abspath() #28,094,847 human reads in total, each read=60

#This function groups the data into fixed-length chunks of size n.
def grouper(n, iterable, fillvalue=None):
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)

n = 1000001 #=>30 files
with open(readsFile) as fin:
    #split file by lines to not cut files in middle of read 
    for i, group in enumerate(grouper(n, fin, fillvalue=''), 1):
        with open('Data\HumanReads_{0}'.format(i * n), 'w') as fout:
            fout.writelines(group)            