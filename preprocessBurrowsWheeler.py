#!/usr/bin/env python
#This file contains the preprocessing of the genome as done by the Burrows Wheeler aligner.

from path import path
import mapperAlignBurrowsWheeler as burrows

genomeFile = path('Data/HumanGenome_Part100Update.gz').abspath()
genome = burrows.readGenome(genomeFile)

preprocFile = path('Preprocessed/BurrowsWheelerGenome').abspath()
with open(preprocFile, 'w') as f:
    bw = burrows.bwt(genome) 
    bwr = burrows.bwt(genome[::-1]) 
    f.write(bw + '\n' + bwr)
