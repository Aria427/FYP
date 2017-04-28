#!/usr/bin/env python
#This file contains the preprocessing of the genome as done by the Burrows Wheeler aligner.

from path import path
import mapperAlignBurrowsWheeler as burrows

genomeFile = path('Data/PhiXGenome.fa.gz').abspath()
genome = burrows.readGenome(genomeFile)

preprocFile = path('Preprocessed/BurrowsWheelerGenome').abspath()

with open(preprocFile, 'w') as f:
    bw = burrows.bwt(genome) 
    #print bw
    bwr = burrows.bwt(genome[::-1]) 
    #print bwr
    f.write(bw + '\n' + bwr)

with open(preprocFile, 'r') as f:
    bw, bwr = f.read().rstrip().split('\n', 1)
    #print bw
    #print bwr
    