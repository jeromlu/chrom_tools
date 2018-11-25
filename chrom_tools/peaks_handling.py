# -*- coding: utf-8 -*-
"""
Created on Fri Oct 26 18:52:57 2018

@author: jeromlu2
"""

import numpy as np
import peakutils
from peakutils.plot import plot as pplot

from chrom_tools import Chromatogram


class PeakIntegrator(object):
    pass
    

class PeakFinder(object):
    '''class that can find peaks in a curve'''
    def __init__(self, chrom = None):
    