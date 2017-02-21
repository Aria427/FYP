#!/usr/bin/env python
#This file includes functions for the LZ77 compression of the genome sequence.

from bitarray import bitarray
import gzip

class LZ77Compressor:
    #This function initialises the class with a default window size of 4.
    def __init__(self, windowSize=4):
        self.windowSize = min(windowSize, 100) 
        #look-ahead buffer contains the next characters to be transmitted
        self.lookaheadSize = 15 #length of match is at most 4 bits
        
    #This function compresses the sequence given.
    def compress(self, data):#inputFile, outputFile):
        encodeBuffer = bitarray(endian='little') #explicitly set endianness for different systems
        
        #with open(inputFile, 'rb') as f:
        #    data = f.read().upper().replace('\n', '').replace('N', '')
    
        i = 0
        while i < len(data):
            match = self.findLongestMatch(data, i)
            
            #Compression format:
            #   No previous matches in window => 0 bit + 8 bits 
            #   Match => 1 bit  
            #               + 12 bits pointer (distance from current position to start of match)
            #               + 4 bits (length of match) 
            
            if match: 	
                (bestMatchDistance, bestMatchLength) = match

                #add 1 bit flag, followed by 12 bit for distance, and 4 bit for length of match 
                encodeBuffer.append(True) #True = 1 bit
                encodeBuffer.frombytes(chr(bestMatchDistance >> 4)) #12 bit for distance
                encodeBuffer.frombytes(chr(((bestMatchDistance & 0xf) << 4) | bestMatchLength)) #4 bit for length of match

                #print "(1, %i, %i)" % (bestMatchDistance, bestMatchLength),
                i += bestMatchLength

            else: #no match
                #add 0 bit flag, followed by 8 bit for character
                encodeBuffer.append(False) #False = 0 bit
                encodeBuffer.frombytes(data[i]) #8 bit for character
				
                #print "(0, %s)" % data[i],
                i += 1

        #bitarray size = 8*n;
        #so if number of bits is not a multiple of 8, fill buffer with 0s
        #encodeBuffer.fill() 
        
        #with open(outputFile, 'wb') as f:
        #    f.write(encodeBuffer.tobytes())

        return encodeBuffer
    
    #This function compresses the sequence found in the input file using compress();
    #The compressed data is written to the output file line by line.   
    def compressFile(self, inputFile, outputFile):
        with gzip.open(inputFile, 'rb') as f:
            #data = ''
            #for line in f:
            #    if line and line[0] != '>':
            #        data = data + line.upper().replace('\n', '').replace('N', '') 
            #encoding = self.compress(data)
            
            #read genome into line-by-line generator     
            subseqs = (line.upper().replace('\n', '').replace('N', '') 
                        for line in f if line and line[0] != '>') #ignore header line
                    
            encoding = bitarray(endian='little')
            for s in subseqs:
                encoding.append(self.compress(s))
            encoding.fill()
            
            with open(outputFile, 'wb') as out:
                out.write(encoding.tobytes())
        
        return encoding
        
    #This function decompresses the sequence found in the input binary file line by line;  
    #The decompressed data is written to the output file.
    def decompress(self, inputFile, outputFile):
        #data = bitarray(endian='little')
        decodeBuffer = []

        with open(inputFile, 'rb') as f:
            for line in f:
                data = bitarray(endian='little')
                data.fromfile(line)

                while len(data) >= 9:
                    flag = data.pop(0)
        
                    if not flag:
                        byte = data[0:8].tobytes()
                        decodeBuffer.append(byte)
                        del data[0:8]
                    else:
                        byte1 = ord(data[0:8].tobytes())
                        byte2 = ord(data[8:16].tobytes())
        
                        del data[0:16]
                        distance = (byte1 << 4) | (byte2 >> 4)
                        length = (byte2 & 0xf)
                        
                        for i in range(length):
                            decodeBuffer.append(decodeBuffer[-distance])
        
        decoding =  ''.join(decodeBuffer)

        with open(outputFile, 'wb') as f:
            f.write(decoding)
	
        return decoding

    #This function finds the longest match of a sequence within a window size.
    def findLongestMatch(self, data, currentPosition):
        bufferEnd = min(currentPosition + self.lookaheadSize, len(data) + 1)
        bestMatchDistance = -1
        bestMatchLength = -1
        
        #output any subsequence of length 1;
        #8 bits uncompressed is better than the 13 bits used for flag compression
        for j in range(currentPosition + 2, bufferEnd): #subsequences of length 2 or greater are considered
            startIndex = max(0, currentPosition - self.windowSize)
            subsequence = data[currentPosition:j]
            
            for i in range(startIndex, currentPosition):
                repetitions = len(subsequence) / (currentPosition - i)
                last = len(subsequence) % (currentPosition - i)

                matchedSequence = data[i:currentPosition] * repetitions + data[i:i+last]

                if matchedSequence == subsequence and len(subsequence) > bestMatchLength:
                     bestMatchDistance = currentPosition - i 
                     bestMatchLength = len(subsequence)

        if bestMatchDistance > 0 and bestMatchLength > 0:
		return (bestMatchDistance, bestMatchLength)
        
        return None