#!/usr/bin/env python
#This file contains the preprocessing of the genome as done by the FM Index aligner.

from path import path
import mapperAlignFmIndex as fm
import pickle

genomeFile = path('Data/PhiXGenome.fa.gz').abspath()
genome = fm.readGenome(genomeFile)

preprocFile = path('Preprocessed/FmIndexGenome').abspath()
with open(preprocFile, 'w') as f:
    index = fm.fmIndex(genome)
    pickle.dump(index, f)
 
with open(preprocFile, 'r') as f:
    index = pickle.load(f)
    