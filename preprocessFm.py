#!/usr/bin/env python
#This file contains the preprocessing of the genome as done by the FM Index aligner.

from path import path
import mapperAlignFmIndex as fm

genomeFile = path('Data/HumanGenome_Part100Update.gz').abspath()
genome = fm.readGenome(genomeFile)

preprocFile = path('Preprocessed/FmIndexGenome').abspath()
with open(preprocFile, 'w') as f:
    index = fm.fmIndex(genome)
    f.write(index)
    