#!/usr/bin/env python
#This file contains the preprocessing of the genome as done by the Kmer Index aligner.

from path import path
import mapperAlignKmerIndex as kmer

genomeFile = path('Data/HumanGenome_Part100Update.gz').abspath()
genome = kmer.readGenome(genomeFile)

preprocFile = path('Preprocessed/KmerIndexGenome').abspath()
with open(preprocFile, 'w') as f:
    index = kmer.KmerIndex(genome, 10)
    f.write(index)
    