#!/usr/bin/env python
#This file includes functions for compressing and decompressing sequences using the inbuilt LZW module.

import lzw

#This function compresses the sequence into a file, using LZW.
def compress(sequenceFile, binaryFile):
    seqBytes = lzw.readbytes(sequenceFile, 100000)
    compressedBytes = lzw.compress(seqBytes)
    lzw.writebytes(binaryFile, compressedBytes)
    
#This function decompresses the sequence into a file, using LZW.
def decompress(binaryFile, sequenceFile):
    compressedBytes = lzw.readbytes(binaryFile, 100000)
    seqBytes = lzw.decompress(compressedBytes)
    lzw.writebytes(sequenceFile, seqBytes)
