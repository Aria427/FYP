#!/usr/bin/env python
#This file includes functions for calculating the algorithmic time complexity of functions.

import timeit
from matplotlib import pyplot
from functools import partial

def plotTC(function, pattern, text, numTests):
    x = []
    y = []
    for i in range(0, 2):
        testTimer = timeit.Timer(partial(function, pattern, text))
        t = testTimer.timeit(number=numTests)
        x.append(i)
        y.append(t)
    return pyplot.plot(x, y, 'o')