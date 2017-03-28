#!/usr/bin/env python
#This file includes functions for data visualisation in the form of files and a Tkinter GUI.

from PIL import Image, ImageDraw
from Bio.Graphics import GenomeDiagram
from Bio.SeqFeature import SeqFeature, FeatureLocation
import Tkinter
import numpy as np
from itertools import groupby
from operator import itemgetter
from collections import Counter
from analyseMapReduceOutput import readReducerOutputs

#This function creates a data visualisation of the matched reads against the genome in a text file.
def visualisationText(genome, readsOffsets, outputFile):
    readsOffsets.sort()                     #sort list
    readsOffsets = sum(readsOffsets, [])    #flatten list
    offsetsCount = Counter(readsOffsets) #record count of each offset => length of match
    file = open(outputFile, 'w')
    for i in range(len(genome)):
        if offsetsCount[i] != 0:
            file.write('-' * offsetsCount[i])
        elif (offsetsCount[i] == 0) & (offsetsCount[i+1] != 0):
            file.write('\n')
            file.write(' ' * (i+1)) #indentation  
        elif (offsetsCount[i] == 0) & (offsetsCount[i+1] == 0):
            file.write(' ')
    file.close()
    return outputFile
  
#This function creates a data visualisation of the matched reads against the genome in a jpg file.
def visualisationJPG(genome, readsOffsets, outputFile):
    readsOffsets.sort()                     #sort list
    readsOffsets = sum(readsOffsets, [])    #flatten list
    offsetsCount = Counter(readsOffsets) #record count of each offset => length of match
    img = Image.new('RGBA', (len(genome), 2500), (255, 255, 255, 0)) 
    draw = ImageDraw.Draw(img)
    draw.line(((0, 10), (len(genome), 10)), fill='blue', width=5)
    j = 0
    for i in range(len(genome)):
        if offsetsCount[i] != 0:
            #draw.line(((i, 20+j), (i+offsetsCount[i], 20+j)), fill='purple', width=5)
            draw.line(((i, 20+j), (i+100, 20+j)), fill='purple', width=5) #100 for PhiX, 60 for Human read
            draw.text((i, 20+j), "%d" % i, fill=0)
            j += 10
    img.save(outputFile, 'JPEG', quality=80, optimize=True, progressive=True)
    return outputFile

#This function creates a data visualisation of the matched reads against the genome using 'GenomeDiagram'.    
def visualisationGD(genome, readsOffsets, outputFile):
    readsOffsets.sort()                     #sort list
    readsOffsets = sum(readsOffsets, [])    #flatten list
    offsetsCount = Counter(readsOffsets) #record count of each offset => length of match
    dia = GenomeDiagram.Diagram(outputFile)
    featuresTrack = dia.new_track(1, greytrack=False)
    featuresSet = featuresTrack.new_set()
    for i in range(len(genome)):
        if offsetsCount[i] != 0:
            #feature = SeqFeature(FeatureLocation(i, i+offsetsCount[i]), strand=+1)
            feature = SeqFeature(FeatureLocation(i, i+100), strand=+1) #100 for PhiX, 60 for Human read
            featuresSet.add_feature(feature, name="%d" % i, label=True, label_size=8) 
    dia.draw(format='linear', pagesize=(len(genome), 2500), fragments=10, start=0, end=len(genome))
    dia.write(outputFile, "PNG")
    return outputFile  

#This function creates a data visualisation of the matched reads against the genome using Tkinter.
def visualisationTkinter(genome, readsOffsets):
    readsOffsets.sort()                     #sort list
    #readsOffsets = sum(readsOffsets, [])    #flatten list
    offsetsCount = Counter(readsOffsets) #record count of each offset => length of match
    window = Tkinter.Tk()
    canvas = Tkinter.Canvas(window, bg='white', width=len(genome), height=500)
    
    def onObjectClick(event): 
        canvas.create_text(50, event.y, text='Clicked at offset %d' % event.x, tag='%d' % event.x)                 
        print 'Clicked at offset', event.x  
    
    canvas.create_line(0, 5, len(genome), 5, fill='red', width=2.5)
    j = 0
    for i in range(len(genome)):
        if offsetsCount[i] != 0:
            #line = canvas.create_line(i, 10+j, i+offsetsCount[i], 10+j, fill='purple', width=2.5) 
            line = canvas.create_line(i, 10+j, i+100, 10+j, fill='purple', width=2.5) #100 for PhiX, 60 for Human read
            canvas.tag_bind(line, '<Button-1>', onObjectClick)
            j += 5
    canvas.pack()
    window.mainloop()

#This function visualises the MapReduce output results using the Tkinter function defined above.   
def visualise(binaryGenomeFile, mapReduceResultFiles):
    with open(binaryGenomeFile , 'rb') as f:
        decodedGenome = np.fromfile(f, dtype=np.int)
        data = readReducerOutputs(mapReduceResultFiles)

        offsets = []
        for currentOffset, group in groupby(data, itemgetter(0)):
            try:
                for offset, count, rq in group:	
                    offsets.append(int(offset))
            except ValueError:
                pass      
    
        visualisationTkinter(decodedGenome, offsets)
    