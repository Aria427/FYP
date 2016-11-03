#!/usr/bin/env python
#This file includes functions for calculating the algorithmic time complexity of functions.

import timeit
from matplotlib import pyplot as plt
from functools import partial

def plotTC(function, pattern, text, numTests):
    testTimer = timeit.Timer(partial(function, pattern, text))
    t = testTimer.timeit(number=numTests)
    return plt.plot(0, t, 'o')