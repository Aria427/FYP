#!/usr/bin/env python
#This file contains the main functionality of the program.

import fileParsing
import alignmentInt
#import analyseAlignment
#import visualisation
import intGenomePatternTesting
import genomePatternTesting
import lz77
import lz78
import lzw

from path import path
import numpy as np
import gzip

import sys
import pdb
import time

#pdb.set_trace()
genomeFile = path('Data\HumanGenome.fa.gz').abspath()
#binaryGenomeFile = path('Output Data\HumanGenomeZip.bin').abspath()
#fileParsing.parseGenomeInt(genomeFile, binaryGenomeFile)
#binaryGenomeFile = path('Output Data\PhixGenome.bin').abspath() 
#fileParsing.parseGenomeInt(path('Data\PhixGenome.fa').abspath(), binaryGenomeFile)

#with open(binaryGenomeFile , 'rb') as f:
    #decodedGenome = np.fromfile(f, dtype=np.int)
    #print decodedGenome #length = 161856144
#d = reduce(lambda x,y: x+str(y), decodedGenome, '')
#decodedGenome = int(d)

#with gzip.open(genomeFile, 'r') as f:
    #data = f.read().rstrip().upper().replace('N', '').replace('\n', '').replace(' ', '')
    
#file read entirely    
from collections import Counter   
with open(path('Data\OneMB.fa').abspath()) as f:
    data = f.read().rstrip().upper().replace('N', '').replace('\n', '').replace(' ', '')
    l = genomePatternTesting.kmerList(data, 4)
    c = Counter(l) 
"""Counter({'TGCT': 12661, 'TGGT': 10773, 'TTTT': 10394, 'GTTT': 9638, 'GCTG': 9637, 'TTAT': 9071, 'TTCT': 8883, 'GCTT': 8881, 'TTTC': 8694, 'AAAA': 8689, 'TATG': 8504, 'TGAT': 8315, 'TTGC': 8127, 'TTTG': 8126, 'TATT': 8126, 'ATTT': 7937, 'TGTT': 7749, 'TTGA': 7558, 'ATGA': 7557, 'CTGA': 7371, 'TAAA': 7370, 'TGAC': 7370, 'CGTT': 7369, 'GGTT': 7181, 'ATGC': 7181, 'TGGC': 7180, 'TTGG': 7180, 'TGAG': 7179, 'ATGG': 6993, 'AATG': 6992, 'GATG': 6991, 'TTCA': 6803, 'CTTC': 6803, 'AAAG': 6802, 'TTAA': 6614, 'TCTG': 6614, 'GATT': 6613, 'GCCG': 6613, 'CTGG': 6613, 'CTTG': 6426, 'TTTA': 6425, 'AAAT': 6423, 'ATTG': 6423, 'CTAA': 6237, 'TCAA': 6235, 'CTGC': 6235, 'TCTT': 6048, 'TACT': 6047, 'GTTG': 5859, 'TAAT': 5858, 'CGCT': 5858, 'CTTT': 5670, 'ACTT': 5669, 'GATA': 5668, 'GTTA': 5668, 'TTAC': 5667, 'ATTA': 5481, 'AAGG': 5292, 'TGGA': 5292, 'GTTC': 5291, 'CGCC': 5291, 'CAAA': 5289, 'GCTA': 5103, 'ACTG': 5102, 'GCGT': 5102, 'GGCT': 5102, 'AAGC': 5101, 'CTTA': 5101, 'GACT': 4914, 'TCGT': 4914, 'ATGT': 4914, 'TTCG': 4913, 'AATC': 4912, 'CCGC': 4725, 'GGTA': 4725, 'TGCC': 4725, 'GTAT': 4724, 'AAGA': 4724, 'AGAA': 4724, 'AATT': 4724, 'ACGC': 4722, 'CTGT': 4722, 'TCGC': 4536, 'CGTC': 4536, 'CATT': 4536, 'GCTC': 4535, 'AGCT': 4535, 'TTCC': 4347, 'TCAG': 4347, 'TCAT': 4347, 'CTAT': 4347, 'ATAA': 4346, 'AGGA': 4346, 'AGGC': 4346, 'GGAT': 4346, 'AGAT': 4344, 'AGGT': 4158, 'GACC': 4158, 'ATTC': 4158, 'GTCA': 4158, 'GGTG': 4158, 'CATG': 4158, 'AGTT': 4157, 'GAAG': 4157, 'CAAG': 4157, 'GAGA': 4155, 'GACG': 4154, 'TCTA': 3969, 'CTCT': 3969, 'AACA': 3969, 'TTGT': 3969, 'GTGA': 3968, 'GAGG': 3968, 'ACCA': 3967, 'GCAA': 3967, 'CCAA': 3966, 'CAAC': 3965, 'CAAT': 3780, 'GTGG': 3780, 'CCGT': 3779, 'GAAA': 3779, 'ATAT': 3779, 'CCGA': 3779, 'TGCG': 3779, 'CTAC': 3591, 'GTCT': 3591, 'ACGT': 3590, 'TGAA': 3590, 'CCTG': 3589, 'GCAG': 3589, 'CGCA': 3589, 'GGCG': 3588, 'CGAC': 3587, 'CCAT': 3402, 'TAAC': 3402, 'ACTA': 3402, 'CTCC': 3402, 'CCTT': 3402, 'ATCA': 3402, 'CAGG': 3401, 'GAGT': 3401, 'TGCA': 3401, 'TACC': 3401, 'CAGA': 3400, 'GCGC': 3399, 'AAAC': 3213, 'GAAT': 3213, 'CGTG': 3213, 'CACT': 3213, 'TCCT': 3213, 'ATCT': 3212, 'GGAA': 3212, 'AACG': 3212, 'TATC': 3212, 'TCTC': 3024, 'ACCG': 3024, 'CAGT': 3024, 'GTCG': 3024, 'CCCT': 3024, 'GTGC': 3023, 'ACAA': 3023, 'GGCA': 2835, 'GCCA': 2835, 'GCCT': 2835, 'GGTC': 2835, 'TGTC': 2834, 'CTCA': 2834, 'TGTG': 2834, 'AACC': 2832, 'AATA': 2646, 'ACTC': 2646, 'AGCG': 2646, 'TCCG': 2646, 'AACT': 2645, 'TACG': 2645, 'ACCT': 2644, 'CGCG': 2644, 'TCCA': 2643, 'GAGC': 2457, 'TAAG': 2457, 'CGGT': 2457, 'GGAC': 2457, 'CCTC': 2457, 'CTCG': 2457, 'CATC': 2457, 'TCCC': 2457, 'AAGT': 2456, 'AGAG': 2455, 'CGAG': 2268, 'ATAC': 2268, 'TCAC': 2268, 'ATCG': 2268, 'CCTA': 2268, 'GCAT': 2268, 'AGTG': 2267, 'CGAT': 2267, 'GACA': 2267, 'CGTA': 2267, 'GTAC': 2079, 'AGTC': 2079, 'GGCC': 2079, 'AGAC': 2079, 'ACCC': 2079, 'ACAC': 2079, 'ACAT': 2079, 'TATA': 2079, 'CCAG': 2078, 'GCGG': 2078, 'GAAC': 2078, 'GGAG': 2078, 'GCGA': 2078, 'ACGA': 2077, 'GTAA': 1890, 'CAGC': 1890, 'CACG': 1889, 'TGTA': 1889, 'TCGA': 1889, 'CGAA': 1701, 'ACAG': 1701, 'GTGT': 1701, 'CCCA': 1701, 'CCAC': 1700, 'ATCC': 1699, 'GCAC': 1512, 'TAGG': 1512, 'CCCC': 1512, 'TACA': 1512, 'GTAG': 1511, 'CGGA': 1511, 'ACGG': 1511, 'GGGT': 1511, 'GTCC': 1511, 'AGCA': 1511, 'AGTA': 1323, 'TTAG': 1323, 'AGGG': 1323, 'TCGG': 1323, 'TAGA': 1323, 'CACC': 1323, 'CGGC': 1322, 'AGCC': 1322, 'TGGG': 1321, 'GCCC': 1134, 'CACA': 1134, 'ATAG': 945, 'TAGT': 945, 'CCGG': 945, 'CCCG': 945, 'CATA': 945, 'GGGA': 944, 'GGGC': 756, 'GGGG': 567, 'CTAG': 567, 'CGGG': 567, 'TAGC': 566})"""  

#file read in chunks so missing patterns on boundary
#print genomePatternTesting.countInt(path('Data\OneMB.fa').abspath())
"""Counter({'TGCT': 12627, 'TGGT': 10742, 'TTTT': 10357, 'GTTT': 9612, 'GCTG': 9608, 'TTAT': 9040, 'GCTT': 8855, 'TTCT': 8853, 'AAAA': 8666, 'TTTC': 8663, 'TATG': 8479, 'TGAT': 8290, 'TATT': 8106, 'TTGC': 8099, 'TTTG': 8098, 'ATTT': 7912, 'TGTT': 7723, 'ATGA': 7539, 'TTGA': 7532, 'CTGA': 7351, 'TGAC': 7348, 'CGTT': 7346, 'TAAA': 7341, 'ATGC': 7164, 'TTGG': 7162, 'TGGC': 7161, 'GGTT': 7160, 'TGAG': 7159, 'GATG': 6975, 'AATG': 6972, 'ATGG': 6969, 'TTCA': 6782, 'CTTC': 6781, 'AAAG': 6776, 'TTAA': 6599, 'TCTG': 6598, 'GATT': 6598, 'CTGG': 6593, 'GCCG': 6590, 'ATTG': 6413, 'CTTG': 6410, 'TTTA': 6407, 'AAAT': 6403, 'TCAA': 6223, 'CTGC': 6219, 'CTAA': 6214, 'TCTT': 6036, 'TACT': 6030, 'GTTG': 5842, 'TAAT': 5841, 'CGCT': 5841, 'CTTT': 5657, 'TTAC': 5652, 'ACTT': 5652, 'GATA': 5650, 'GTTA': 5650, 'ATTA': 5466, 'TGGA': 5279, 'CGCC': 5279, 'AAGG': 5275, 'CAAA': 5275, 'GTTC': 5273, 'ACTG': 5090, 'GCGT': 5088, 'AAGC': 5087, 'CTTA': 5087, 'GCTA': 5085, 'GGCT': 5082, 'ATGT': 4903, 'TCGT': 4902, 'AATC': 4898, 'GACT': 4895, 'TTCG': 4895, 'CCGC': 4711, 'AGAA': 4711, 'GTAT': 4710, 'AATT': 4709, 'AAGA': 4706, 'CTGT': 4706, 'GGTA': 4706, 'ACGC': 4705, 'TGCC': 4705, 'AGCT': 4524, 'TCGC': 4523, 'CATT': 4523, 'GCTC': 4519, 'CGTC': 4516, 'GGAT': 4339, 'AGGA': 4338, 'TCAG': 4337, 'ATAA': 4335, 'AGAT': 4335, 'TCAT': 4335, 'TTCC': 4334, 'CTAT': 4334, 'AGGC': 4332, 'GTCA': 4150, 'CAAG': 4148, 'GGTG': 4148, 'GACC': 4147, 'CATG': 4147, 'AGTT': 4146, 'AGGT': 4146, 'ATTC': 4146, 'GACG': 4144, 'GAAG': 4143, 'GAGA': 4142, 'TCTA': 3960, 'CTCT': 3959, 'GTGA': 3957, 'CCAA': 3956, 'ACCA': 3956, 'AACA': 3956, 'CAAC': 3955, 'TTGT': 3955, 'GCAA': 3955, 'GAGG': 3954, 'CAAT': 3772, 'GTGG': 3769, 'CCGA': 3769, 'GAAA': 3767, 'ATAT': 3767, 'TGCG': 3767, 'CCGT': 3764, 'CTAC': 3581, 'GCAG': 3581, 'GTCT': 3579, 'CGAC': 3577, 'ACGT': 3577, 'TGAA': 3577, 'CGCA': 3577, 'GGCG': 3576, 'CCTG': 3576, 'ATCA': 3395, 'CAGG': 3395, 'CCTT': 3393, 'TAAC': 3392, 'CAGA': 3391, 'GCGC': 3391, 'TACC': 3391, 'CTCC': 3390, 'GAGT': 3390, 'CCAT': 3389, 'ACTA': 3388, 'TGCA': 3387, 'TCCT': 3207, 'CGTG': 3206, 'AAAC': 3204, 'GAAT': 3204, 'ATCT': 3204, 'GGAA': 3204, 'CACT': 3204, 'AACG': 3200, 'TATC': 3200, 'CAGT': 3017, 'GTGC': 3016, 'CCCT': 3015, 'GTCG': 3014, 'TCTC': 3013, 'ACAA': 3013, 'ACCG': 3011, 'TGTG': 2828, 'GGTC': 2826, 'GGCA': 2825, 'GCCA': 2825, 'GCCT': 2825, 'TGTC': 2823, 'AACC': 2823, 'CTCA': 2822, 'AGCG': 2639, 'AATA': 2638, 'ACTC': 2638, 'AACT': 2638, 'CGCG': 2637, 'TCCG': 2637, 'TCCA': 2637, 'TACG': 2636, 'ACCT': 2633, 'GAGC': 2453, 'GGAC': 2451, 'CGGT': 2450, 'CTCG': 2450, 'TAAG': 2449, 'CCTC': 2449, 'AGAG': 2448, 'CATC': 2448, 'AAGT': 2447, 'TCCC': 2445, 'CCTA': 2265, 'CGAG': 2262, 'CGAT': 2261, 'ATCG': 2261, 'CGTA': 2261, 'TCAC': 2260, 'GACA': 2259, 'GCAT': 2259, 'ATAC': 2258, 'AGTG': 2257, 'ACCC': 2075, 'AGAC': 2074, 'ACAC': 2074, 'TATA': 2074, 'AGTC': 2073, 'ACGA': 2072, 'GGCC': 2072, 'GGAG': 2072, 'GTAC': 2071, 'CCAG': 2071, 'ACAT': 2071, 'GAAC': 2071, 'GCGG': 2070, 'GCGA': 2070, 'TGTA': 1885, 'CAGC': 1884, 'CACG': 1883, 'GTAA': 1882, 'TCGA': 1882, 'GTGT': 1697, 'CCAC': 1696, 'ACAG': 1695, 'CGAA': 1694, 'ATCC': 1693, 'CCCA': 1692, 'AGCA': 1509, 'GTAG': 1508, 'TAGG': 1508, 'TACA': 1508, 'GCAC': 1507, 'GGGT': 1507, 'CGGA': 1506, 'ACGG': 1505, 'GTCC': 1505, 'CCCC': 1505, 'AGTA': 1321, 'TTAG': 1321, 'TAGA': 1321, 'CACC': 1320, 'AGGG': 1319, 'TCGG': 1319, 'CGGC': 1319, 'TGGG': 1317, 'AGCC': 1316, 'GCCC': 1130, 'CACA': 1130, 'ATAG': 943, 'TAGT': 943, 'GGGA': 943, 'CCGG': 942, 'CCCG': 942, 'CATA': 940, 'GGGC': 754, 'GGGG': 566, 'CGGG': 566, 'CTAG': 565, 'TAGC': 564})"""

print genomePatternTesting.countInt(genomeFile)
#has bug when run on human genome -> Counter({'CHR1': 1, '>CHR': 1}) 
#=> only reads header line


reads = fileParsing.parseReadsInt(path('Data\HumanSequencingReads.tsv.bz2').abspath()) 
#reads = fileParsing.parseReadsPhiXInt(path('Data\PhiXSequencingReads1000.fastq').abspath())


#analyseAlignment.plotTimeVsMatches(reads, decodedGenome, path('Output Test Files\AlignmentAnalysis.png').abspath())


#matchesCount, totalCount, offsets = alignmentInt.alignHamming(reads, decodedGenome)
#print '%d/%d reads matched the genome.' % (matchesCount, totalCount) #The result is not 100% but this is to be expected due to sequencing errors. 
#print offsets


"""
textFile = path('Output Test Files\DataVisualisationTest.txt').abspath()
textFile = visualisation.visualisationText(genome, offsets, textFile)

jpgFile = path('Output Test Files\DataVisualisationTest.jpg').abspath()
jpgFile = visualisation.visualisationJPG(genome, offsets, jpgFile) 

pngFile = path('Output Test Files\DataVisualisationTest.png').abspath()
pngFile = visualisation.visualisationGD(genome, offsets, pngFile)
"""
#visualisation.visualisationTkinter(decodedGenome, offsets)
