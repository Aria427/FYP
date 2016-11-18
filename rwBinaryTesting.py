#!/usr/bin/env python
#This file contains some testing functions with regards to writing and reading a binary file.

import fileParsing

import sys
import struct
from bitarray import bitarray
import numpy

print("Native byteorder: ", sys.byteorder)

with open('Output Test Files\longlong.bin', 'wb') as f:
    sequence = 'GAGTTTTATCGCTTCCATGACGCAGAAGTTA'
    print 'Length of long long sequence: %d' % len(sequence)
    binary = fileParsing.baseToBinary(sequence)
    print 'Length of binary converted long long sequence: %d' % len(binary)
    x = int(binary, 2)
    print 'Byte sequence: %d' % x
    print 'Bit length of byte sequence: %d' % x.bit_length()
    f.write(struct.pack('q', x)) #long long - 8 bytes

with open('Output Test Files\int.bin', 'wb') as f:
    sequence = 'GAGTTTTATCGCTTC'
    print 'Length of int sequence: %d' % len(sequence)
    binary = fileParsing.baseToBinary(sequence)
    print 'Length of binary converted int sequence: %d' % len(binary)
    x = int(binary, 2)
    print 'Byte sequence: %d' % x
    print 'Bit length of byte sequence: %d' % x.bit_length()
    f.write(struct.pack('i', x)) #int - 4 bytes

with open('Output Test Files\Bitarray.bin', 'wb') as f:
    sequence = 'GAGTTTTATCGCTTC'
    ba = bitarray()
    ba.encode(fileParsing.bases, sequence)
    ba.tofile(f)

print '\n'

with open('Output Test Files\int.bin' , 'rb') as f:
    fileBytes = f.read(4) #read first 4 bytes 
    data = fileBytes
    while fileBytes != '':
        fileBytes = f.read(4) #read the rest 4 bytes by 4 bytes
        data += fileBytes 
    print 'Encoded data: %s' % data
    print 'Representation of encoded data: %s' % repr(data)
    print 'Hexadecimal format of encoded data: %s' % data.encode('hex')
    print 'Integer format of encoded data: %s' % ' '.join([str(ord(a)) for a in data])
    for d in data:
        i = ord(d) #integer value of byte
        print 'Integer value of byte: %s' % i
        bin = '{0:b}'.format(i) 
        print 'Binary decoding of byte: %s' % bin
        hex = '{0:x}'.format(i) 
        print 'Hexadecimal decoding of byte: %s' % hex
        oct = '{0:o}'.format(i) 
        print 'Octal decoding of byte: %s' % oct
    for i in range(0, len(data), 4):
        fileBytes = struct.unpack('i', data[i:i+4])[0] #[0] to remove unnecessary syntax  
        print 'Unpacked data range: %s' % fileBytes

with open('Output Test Files\int.bin' , 'rb') as f:
    fileBytes = struct.unpack('i', f.read(4))[0]
    print 'Unpacked data direct: %s' % fileBytes 
        
with open('Output Test Files\int.bin', 'rb') as f:
    for chunk in iter(lambda: f.read(4), ''):
        fileBytes = struct.unpack('i', chunk)[0]        
        print 'Unpacked data lambda: %s' % fileBytes
    
with open('Output Test Files\int.bin', 'rb') as f:
    fileBytes = numpy.fromfile(f, dtype=numpy.int)
    print 'Unpacked data numpy: %s' % fileBytes
    
ba = bitarray()
with open('Output Test Files\Bitarray.bin', 'rb') as f:
    ba.fromfile(f)
    print ba
    