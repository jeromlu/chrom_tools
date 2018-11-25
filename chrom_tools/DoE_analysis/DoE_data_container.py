# -*- coding: utf-8 -*-
"""
Created on Thu Oct 11 21:52:44 2018

@author: JEROMLU2
"""

class TargetResultsData(object):
    '''
    Class to contain one target's fit results together with fitting history
    '''
    
    def __init__(self, parent = None):
        super(TargetResultsData, self).__init__(parent)
        
        #all results
        self.results = None
        
        #name of fitted target
        self.target_name = None
        
        #list of history, if there was stepvise fit
        self.history = None
        
        
class TargetsResultsContainer(object):
    
    def __init__(self):
        
        #pd.DataFrame to contain more TargetResultsData objects
        self.
        
        
        