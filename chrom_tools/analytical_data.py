# -*- coding: utf-8 -*-
"""
Created on Tue Sep 11 08:17:33 2018

@author: jeromlu2
"""

import logging
logger = logging.getLogger(name = 'analytical_data.py')

#import of python packages
import os
import sys
import ntpath
import re

#import of third party packages
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms

#import of my or modified packages


class AnalyticalData(object):
    
    
    HEADER_ROW = 1
    INDEX_COL = 2
    
    def __init__(self, fname = None):
        
        
        
        self.loads = pd.DataFrame([])
        
        self.eluates = pd.DataFrame([])
        
        if fname is not None:
            self.fname = fname
            self.load_data()
            
            
    def load_data(self):
        '''Loads the excel file into the Analytical data data structure
        '''
        
        na_values = ['n.a.']
        header_row = 1
        index_col = 2
        
        try:
            self.eluates = pd.read_excel(self.fname,
                                       sheet_name = 'Eluates',
                                       header = header_row,
                                       index_col = index_col,
                                       na_values = na_values)
            
            self.loads = pd.read_excel(self.fname,
                                       sheet_name = 'Loads',
                                       header = header_row,
                                       index_col = index_col,
                                       na_values = na_values)   
            if self.loads.index.is_unique == False:
                print(self.loads[self.loads.index.duplicated()])
                raise IndexError('Index of loads is not unique.')
            #if self.eluates.index.is_unique == False:
            #    print(self.eluates.index[self.eluates.index.duplicated()])
            #    raise IndexError('Index of eluates is not unique.')
                
        except Exception as e:
            self.print_err()
            
    #funkcija da doda oznake nad stolpce, imam dve opciji (yield ali pa visino stolpca)
    def autolabel(self, ax, rects, values = False):
        '''Function that attaches some text labels to the rectangles
        '''
        span = max(ax.get_ylim())-min(ax.get_ylim())
        if len(values) <1:
            for rect in rects:
                height = rect.get_height()
                ax.text(rect.get_x()+rect.get_width()/2., height + 0.03 * span , '{0:0.0f}'.format(height),
                        ha='center', va='bottom',fontsize = 12) 
        else:
            for rect,value in zip(rects,values):
                if type(value) is str:
                    text = value
                else:
                    text = '{0:0.0f} %'.format(value)
                #text = str(value)
                height = rect.get_height()
                ax.text(rect.get_x()+rect.get_width()/2., height + 0.03 * span, text,
                        ha='center', va='bottom', fontsize = 12,
                        bbox={'facecolor':'white', 'alpha' :0.85, 'edgecolor':'none', 'boxstyle' :'round,pad=0.2'})   
              
            
    def plot_load_eluate_CQA(self):
        '''Plots selected samples (eluates) in comparison to the starting material
        '''
        pass
    
    def plot_sequence_linked(self, seq_el_IDs, quality_attribute, step_names = None):
        '''plots linked results, uses PyEnergyDiagrams
        '''
        from pyEnergyDiagram.energydiagram import ED
        diagram = ED()
        
        if quality_attribute not in self.loads.columns:
            print('Qulity attribute not in the loads!!')
            print('You may select:')
            print(self.loads.columns)
            return
        
        if quality_attribute not in self.eluates.columns:
            print('Qulity attribute not in the eluates!!')
            print('You may select:')
            print(self.eluates.columns)
            return
        
        el_data = self.eluates.loc[seq_el_IDs]
        self.eluates['el_start_val'] = el_data.loadID.map(self.loads[quality_attribute])
        el_data = self.eluates.loc[seq_el_IDs]
        print(el_data['el_start_val'])
        
        #create levels
        for i, ser_iterator in enumerate(el_data.iterrows()):
            idx, row = ser_iterator
            load = row['el_start_val']
            diagram.add_level(load, row.loadID)
            eluate = row[quality_attribute]
            diagram.add_level(eluate, idx, 'last')
            print(i,i+1)
            
            #add links and arrows
            if (i > 0):
                if((i % 2) == 1):
                    diagram.add_link(i, i+1)
                    diagram.add_arrow(i-1, i)

                    
            
            if step_names is not None:
                print('TODO')

        
        diagram.plot(show_IDs=True)
        return diagram
        
        
                    
    def print_err(self):
            exc_type, exc_obj, exc_tb = sys.exc_info()           
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            err_msg = '{0}:\n{1}\nError occurred in file: {2}'.format(exc_type, exc_obj, fname)
            print(err_msg)
            logger.error(err_msg)     
            
            


class OriginatorData(object):
    '''Class for keeping the data bout the reference biologics
    '''
    
    def __init__(self, fname = None, *args, **kwargs):
        
        self.fname = fname
        
        self.originator = pd.DataFrame([])

    
    def load_data(self):
        
        try:
        
            raw_orig = pd.read_excel(self.fname, skiprows = [0])
            
            self.originator['average'] = raw_orig.mean()
            self.originator['sigma'] = raw_orig.std()
            self.originator['lower'] = self.originator.average - 3 * self.originator.sigma
            self.originator['upper'] = self.originator.average + 3 * self.originator.sigma
        except Exception as e:
            print(e)

        


            
            
if __name__ == '__main__':
    
    folder = 'C:/Users/jeromlu2/LukaFiles/04_Programiranje/01_Python/03_IPythonNotebooks/temp/IPI/03_DSP Development/04_DSPD/01_Akta_runs/'
    fname = 'GNT904_data_akta_runs.xlsx'
    
    data = AnalyticalData(folder + fname)
    diagram = data.plot_sequence_linked(['1Z918','1Z920'],'SEC_HMWs')
    diagram2 = data.plot_sequence_linked(['1Z918'],'SEC_HMWs')
    
    
                       