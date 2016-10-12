"""
This file includes functions for data visualisation in the form of a text file and a jpeg file.
"""

import collections
from PIL import Image, ImageDraw

#To create a data visualisation of the matched reads against the genome in a text file:
def visualisationText(genome, readsOffsets, outputFile):
    readsOffsets.sort()                     #sort list
    readsOffsets = sum(readsOffsets, [])    #flatten list
    offsetsCount = collections.Counter(readsOffsets) #record count of each offset => length of match
    file = open(outputFile, 'w')
    for i in range(len(list(genome))):
        if offsetsCount[i] != 0:
            file.write('-' * offsetsCount[i])
        elif (offsetsCount[i] == 0) & (offsetsCount[i+1] != 0):
            file.write('\n')
            file.write(' ' * (i+1)) #indentation  
        elif (offsetsCount[i] == 0) & (offsetsCount[i+1] == 0):
            file.write(' ')
    file.close()
    return outputFile
  
#To create a data visualisation of the matched reads against the genome in a jpg file:
def visualisationJPG(genome, readsOffsets, outputFile):
    readsOffsets.sort()                     #sort list
    readsOffsets = sum(readsOffsets, [])    #flatten list
    offsetsCount = collections.Counter(readsOffsets) #record count of each offset => length of match
    img = Image.new('RGBA', (len(list(genome)), 2500), (255, 255, 255, 0)) 
    draw = ImageDraw.Draw(img)
    draw.line(((0, 10), (len(list(genome)), 10)), fill='blue', width=5)
    j = 0
    for i in range(len(list(genome))):
        if offsetsCount[i] != 0:
            draw.line(((i, 20+j), (i+offsetsCount[i], 20+j)), fill='purple', width=5)
            j += 10
    img.show() #display generated image
    img.save(outputFile, 'JPEG', quality=80, optimize=True, progressive=True)
    return outputFile
