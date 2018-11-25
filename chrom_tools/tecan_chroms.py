# -*- coding: utf-8 -*-
"""
Created on Tue Sep 11 20:51:48 2018

@author: JEROMLU2

Defines class for handling Trinean measurements.
"""

#import of python packages
import bisect
import ntpath
import re
import os
import sys

#import of third party packages
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms

#my packages
import chrom_tools

CODEC = 'utf-16'
FRAC_SIZE = 0.3 #mL


#how to fill chromatogram from Tecan data....

class TecanChrom(object):
    
    def __init__(self, name = '', units = None, curves = pd.DataFrame([]), cv = np.nan,
                 fractions = None, phases = None, injection_marks = None, run_log = '',
                 flow = 1):
        
        self.name = name
        
        #TODO: list of units for each curve
        self.units = units
        
        self.curves = curves
        
        self.signals = [col_name for col_name in curves.columns if 'vol' not in col_name]
      
        self.cv = cv
        
        self.fractions = fractions
        
        self.phases = phases
        
        self.injection_marks = injection_marks
        
        self.run_log = run_log
        
        #sampling data and analytical data
        self.sampling_data = {}
        
    def get_x_values(self, signal, offset = 0, cv = 1):
        '''returns x values of the signal (pd.Series)\n
        erases nan values
        '''
        prefix = 'vol'
        return (self.curves[prefix + signal].dropna() - offset) / cv
    
    def get_y_values(self, signal):
        '''returns signal values (pd.Series)
        erases nan values
        '''
        return self.curves[signal].dropna()
        
class TecanChroms(object):
    '''Tecan chromatograms container 
    '''
    
    def __init__(self):
        self.__fname = ""
        self.__chromatograms = []
        self.__chromFromID = {}
        self.__dirty = False
        self.__currentCurve = np.array([]) #tole je treba se nastudirat
        
    def __len__(self):
        return len(self.__chromatograms)
    
    def __iter__(self):
        for pair in self.__chromatograms:
            yield pair[1]
            
    def __getitem__(self, k):
        if type(k) == str:
            names = self.get_chrom_names()
            if k in names:
                idx = bisect.bisect_left(names, k)
                return self.__chromatograms[idx][1]
            else:
                return None
        elif (type(k) == list) and (type(k[0]) == str):
            return [[name, self[name]] for name in k]
        elif type(k) == int:
            return self.__chromatograms[k][1]
        else:
            return self.__chromatograms[k]        
    
    def __contains__(self, name):
        '''
        Checks if chromatogram is already in the chromatograms container.
        '''
        if name in self.get_chrom_names():
            return True
        else:
            return False
        
    def get_chrom_names(self):
        '''Return list of chrom names'''
        return [el[0] for el in self.__chromatograms]
    
    def chromatograms(self):
        return self.__chromatograms
    
    def filename(self):
        return self.__fname  
      
    def add_chromatogram(self, chrom):
        '''add chromatogram object to TecanChroms
        '''
        
        if id(chrom) in self.__chromFromID:
            return False
        key = chrom.name
        bisect.insort_left(self.__chromatograms, [key, chrom])
        self.__chromFromID[id(chrom)] = chrom
        self.__dirty = True
        return True  

      
    def load_from_Trinean_xls(self, fname, start_salts = None, step_lengths = None,
                              end_salts = None):
        '''
        Loads Concnetration curves from Trinean exported excel file.
        Sample names have to be specific... Created with this application.
        It loads up to eight curves, depending on how many columns were used
        during the Tecan screening.
        
        Parameters
        ----------
        fname : str
            name of the file. It has to be full path with file extension at the end 
            or it can be relative path with file extansion at the end.
        start_salts : list or None
            list of floats, each float represents the starting salt molarity
            during the each step
        step_lengths : list or None
            list of the step lengths, number of fractions in each step
        end_salts : list or None
            list of salt concentration at the end of the step. Could be None
            if there is no gradient.
            
        Returns
        -------
        succes : bool
            True if succesfuly loaded
            False if there were some problems, or file was already loaded
        msg : str
            descriptive message about the problem or succes
        '''        
        
        error = None
        try:
            nan_val = ['-']
            folder, file = os.path.split(fname)
            data_raw = pd.read_excel(fname, na_values= nan_val, index_col = 0,
                                 sheet_name = 'Trinean Report 1', header = 19)
            data = data_raw.dropna(axis = 0, subset = ['Avg\nA280 Concentration (mg/ml)'])
            data['Sample name'] = data.index
            data = data[~data['Sample name'].str.contains('blank')]
            cols = ['Sample name','Avg\nA280 Concentration (mg/ml)']
            print(data.shape)
    
            #razdelim po ploscicah
            split_df = pd.DataFrame(data['Sample name'].str.split().tolist(), 
                            columns = ['Fractions', 'x', 'Column', 'Sample_mark'], index = data['Sample name'])
            mask = split_df['Column'] == '1'
        
            chromatograms = pd.DataFrame(columns = split_df[mask]['Fractions'].values)
            for i in range(1,9):
                mask = split_df['Column'] ==str(i)
                print(mask.shape)
                chromatograms.loc['Column_' + str(i)] = data[mask]['Avg\nA280 Concentration (mg/ml)'].values
            #sortiranje po indexu, mogoce ne rabim
            chromatograms.loc['num',:] = pd.to_numeric(pd.DataFrame(chromatograms.columns.str.split('_').tolist())[1]).values
            chromatograms.sort_values('num', axis = 1, inplace = True)
            chromatograms.drop('num', inplace = True)
            
            for idx, curve in chromatograms.iterrows():
                name = file.split('.')[0]
                nr_points = len(curve)
                curve_vol = FRAC_SIZE * pd.Series(np.arange(nr_points))
                if start_salts:
                    salt_x, salt_y = self.salt_conc(start_salts, step_lengths, end_salts) #[0.1, 0.2, 0.5], [4,7,4] , [None, 0.5, None])
                    print(salt_x)
                    print(np.arange(nr_points))
                    print(len(salt_x),len(salt_y), len(curve_vol), len(curve))
                    curves = pd.DataFrame([],columns = ['volConc', 'Conc', 'volCond', 'Cond'])
                    curves['volCond'] = salt_x * FRAC_SIZE
                    curves['Cond'] = salt_y
                    fill = pd.Series([np.nan]*(len(salt_y) - nr_points))
                    curves['volConc'] = curve_vol.append(fill).values
                    curves['Conc'] = curve.append(fill).values
                else:
                    curves = pd.DataFrame(np.column_stack([curve_vol, curve]),
                                          columns = ['volConc', 'Conc'])
                chrom = TecanChrom(name + '_' + idx, curves = curves)
                print(chrom.name)
                self.add_chromatogram(chrom)
            return True, chromatograms
        
        except Exception as error:
            self.print_err()
            
        finally:
            if error is not None:
                return False, error
            
            msg = 'Loaded additional chromatogram: {0}'.format(fname)
            return True, msg
        
    def salt_conc_step(self, start_salt_conc, start_fr, end_fr, end_salt_conc = None):
        '''
        Calculates the salt concentration in a step
        
        Parameters
        ---------
        start_salt_conc : float
            concentration at the beginning of the step
        start_fr : int
            at which fraction the step begins
        end_fr : int
            at which fraction step ends
        end_salt_conc : float or None
            if None ther is no gradient, if float gradient goes up to this value
        '''
        #creates independant values   
        x = np.arange(start_fr, end_fr + 1)
        
        #if end_salt_conc is not None it calculates gradient parameters
        if end_salt_conc:
            k = (end_salt_conc - start_salt_conc) / (end_fr - start_fr)
            n = (start_salt_conc * end_fr - start_fr * end_salt_conc) / (end_fr - start_fr)
            y = n + k * x
        #else there is no gradient only constant value
        else:
            y = np.ones(len(x)) * start_salt_conc
        return x, y
    
    def salt_conc(self, start_salts, step_lengths, end_salts):
        '''
        Calculates and creates salt concentrations
        Uses salt_conc_step function
        
        Parameters
        ----------
        start_salts : list
            list of floats, each float represents the starting salt molarity
            during the each step
        step_lengths : list
            list of the step lengths, number of fractions in each step
        end_salts : list
            list of salt concentration at the end of the step. Could be None
            if there is no gradient.
        '''
      
        xs = []
        ys = []
        
        #I do some initial checking to see if the lengths of the list are OK
        if len(start_salts) != len(end_salts) or len(step_lengths) != len(start_salts):
            print('lists should be same size')
            return None
        #points where change, next step occures
        changes = np.r_[0,np.cumsum(step_lengths)]

        
        #for each step I calculate the salt concentration and add it to appropriate list
        for i, start in enumerate(start_salts):
            x_i, y_i = self.salt_conc_step(start, changes[i], changes[i+1], end_salts[i])
            xs.append(x_i)
            ys.append(y_i)
        
        #to get final curve I concatienate all the lists
        x = np.concatenate(xs)
        y = np.concatenate(ys)
        #plt.plot(x, y)
        
        return x, y

    @staticmethod
    def formats():
        return "*.xlsx"
    
    def print_err(self):
            exc_type, exc_obj, exc_tb = sys.exc_info()           
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            err_msg = '{0}:\n{1}\nError occurred in file: {2}\n at line: {3}'.format(exc_type, exc_obj, fname, exc_tb.tb_lineno)
            print(err_msg)    
    
    
if __name__ == '__main__':
    
    folder = 'C:/Users/jeromlu2/LukaFiles/04_Programiranje/01_Python/03_IPythonNotebooks/temp/IPI/03_DSP Development/04_DSPD/02_Tecan_screenings/04_CEX/01_podatki/Basel_test/TS1456/'
    file = 'TS1456.xlsx'
    
    t_chroms = TecanChroms()
    chromatograms = t_chroms.load_from_Trinean_xls(folder + file)
    
    
    