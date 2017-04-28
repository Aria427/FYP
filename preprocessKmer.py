#!/usr/bin/env python
#This file contains the preprocessing of the genome as done by the Kmer Index aligner.

from path import path
import mapperAlignKmerIndex as kmer
import pickle

genomeFile = path('Data/PhiXGenome.fa.gz').abspath()
genome = kmer.readGenome(genomeFile)

preprocFile = path('Preprocessed/KmerIndexGenome').abspath()
with open(preprocFile, 'w') as f:
    index = kmer.KmerIndex(genome, 10)
    pickle.dump(index, f)
    
with open(preprocFile, 'r') as f:
    index = pickle.load(f)    
    