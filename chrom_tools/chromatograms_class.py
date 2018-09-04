# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 21:12:55 2018

@author: JEROMLU2
"""
#import of python packages
import bisect
import ntpath
import re

#import of third party packages
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms

CODEC = 'utf-16'

def colorAtNumber(n):
    return list(mpl.rcParams['axes.prop_cycle'])[n]['color']


class Chromatogram(object):
    
    def __init__(self, name = '', units = None, curves = pd.DataFrame([]), cv = np.nan,
                 fractions = None, phases = None, injection_marks = None, run_log = '',
                 flow = 1):
        
        self.name = name
        
        self.units = units
        
        self.curves = curves
        
        self.signals = [col_name for col_name in curves.columns if 'vol' not in col_name]
      
        self.cv = cv
        
        self.fractions = fractions
        
        self.phases = phases
        
        self.injection_marks = injection_marks
        
        self.run_log = run_log
        
        #preverim ce jih imam
        self.all_alignments = self.get_all_alignments()
        
        self.flow = flow
        
    def get_all_alignments(self):
        '''finds all possible alignments
        
        merges injections, phases, and fractions
        '''        
        #print(pd.concat([self.fractions, self.phases ,self.injection_marks]))
        if (self.fractions is None) and (self.phases is None) and (self.injection_marks is None):
            return None
        else:
            return pd.concat([self.fractions, self.phases ,self.injection_marks])

    
    def get_offset(self, name):
        '''returns the volume at which specific (name) alignment was created
        '''
        offsets = self.get_all_alignments()
        if name in offsets:
            return offsets.loc[name]
        else:
            print('Na izbiro imas:\n', offsets) 
    
    def get_run_ID(self):
        '''Returns the D-number and scouting number of chromatogram
        '''
        match = re.match('\w+-(\w+)-\w+-\w+\s+(\d\d\d)', self.name)
        return match.group(1), match.group(2)
    
    def find_nearest_idx(self, array, value):
        '''vrne indeks vrednosti v tabeli array, ki je najblizje isakni vrednosti (value)
        Primer:
        tab = np.arange(15)
        find_nearest(tab, 6.6)
        result = 7
        '''
        idx = (np.abs(array-value)).argmin()
        return idx
    
    def nearest_at_vol(self, signal, value):
        '''Returns index in column of selected signal that is closest to the given value
        '''
        x = self.get_x_values(signal)
        idx = self.find_nearest_idx(x, value)
        y = self.get_y_values(signal)
        return y[idx]
    
    def get_x_values(self, signal):
        '''returns x values of the signal (pd.Series)\n
        erases nan values
        '''
        prefix = 'vol'
        return self.curves[prefix + signal].dropna()
    
    def get_y_values(self, signal):
        '''returns signal values (pd.Series)
        erases nan values
        '''
        return self.curves[signal].dropna()
    
    def get_curve_cut(self, signal, cut):
        '''returns only part of the chromatogram signal, based on the volume limits
        one can also provide string as a name for a limit
        
        signal....  ime signala [npr. Cond, UV_1 280, UV_2 300] (string)
        cut.......  info about part of the curve that you are interested in (list)
                    which may contain strings or numbers
        '''
        limits = []
        if len(cut) > 2:
            limits_raw = cut[:2]
        else:
            limits_raw = cut
        
        for value in limits_raw:
            if type(value) == str:
                if (self.all_alignments is not None):
                    if value in self.all_alignments:
                        limits.append(self.all_alignments[value])
                    else:
                        print('You may only select: ', self.all_alignments.index)
                else:
                    print('There are no named alignments possible')
            else:
                limits.append(value)
                
        x = self.get_x_values(signal)
        mask =  (x > min(limits)) & (x < max(limits))
        
        x_out = x[mask]
        y_out = self.get_y_values(signal)[mask]
        
        return x_out, y_out
    
    def maxima(self, signal, n , x_limit = None, treshold = None):
        '''Returns maxima (np.array) of chromatogram\n
        white noise is a problem:
        - change n
        - only maxima above tresh
        - have to implement slope detection
        '''
        x = self.get_x_values(signal).values
        y = self.get_y_values(signal).values
        if x_limit:
            mask = (x < max(x_limit)) & (x > min(x_limit))
            x = x[mask]
            y = y[mask]
        
        mask= np.r_[True, y[1:] > y[:-1]] & \
                                    np.r_[y[:-1] > y[1:], True]
        iteration_list = range(n+1)
        for i in iteration_list[2:]:
            a=np.ones(i,dtype=bool)
            new_mask= np.r_[a, y[i:] > y[:-i]] & \
                                        np.r_[y[:-i] > y[i:], a]
            mask=mask & new_mask
        if treshold:
            tresh_mask = y > treshold
            mask = mask & tresh_mask
            
        peaks = np.array(list(zip(x[mask],y[mask])))
        return peaks
    
    def minima(self, signal, n, x_limit = None, treshold = None):
        '''Returns minima (np.array) of chromatogram\n
        white noise is a problem:
        - change n
        - only maxima above tresh
        - have to implement slope detection
        '''
        
        x = self.get_x_values(signal).values
        y = self.get_y_values(signal).values
        if x_limit:
            mask = (x < max(x_limit)) & (x > min(x_limit))
            x = x[mask]
            y = y[mask]
        
        
        mask= np.r_[True, y[1:] < y[:-1]] & \
                                    np.r_[y[:-1] < y[1:], True]
        iteration_list = range(n+1)
        for i in iteration_list[2:]:
            a=np.ones(i,dtype=bool)
            new_mask= np.r_[a, y[i:] < y[:-i]] & \
                                        np.r_[y[:-i] < y[i:], a]
            mask=mask & new_mask
            
        if treshold:
            tresh_mask = y > treshold
            mask = mask & tresh_mask
            
        peaks = np.array(list(zip(x[mask],y[mask])))
        return peaks
    
    def remove_baseline(self, signal):
        '''deterines the base line of a signal
        
        not implemented yet'''
        pass
    
    def curve_stat_moments(self, signal, x_limit, base_line):
        '''calculates first four moments of curve\n
        by means of numerical integration\n
        x_limit..... limits the area of integration'''
        
        x = self.get_x_values(signal).values
        y = self.get_y_values(signal).values
        if x_limit:
            mask = (x < max(x_limit)) & (x > min(x_limit))
            x = x[mask]
            y = y[mask]
            
        y = y - base_line
        mu_0 = np.trapz(y, x)
        mu_1 = (np.trapz(y*x, x)) / mu_0
        mu_2 = (np.trapz(y * (x-mu_1)**2, x)) / mu_0
        skew = (np.trapz(y * (x-mu_1)**3, x)) / (mu_0 * mu_2**3)
        print('stat moments: mu_0: {0:.1f}, mu_1: {1:.1f}, mu_2: {2:.2f}'.format(mu_0, mu_1, mu_2))
        return mu_0, mu_1, mu_2, skew
    
    def HETP_calc(self, signal, x_limit, L, base_line = None, plot = False):
        '''Calculates HETP using three different methods (Jungbauer book) of selected signal\n
        signal........selceted signal (str)\n
        x_limit.......upper and lower volume limit of the curve (list)\n
        L.............column length in cm (float)\n
        base_line.....base line of the curve (float)\n
        plot..........plots curve that was used for determination of the HETP
        '''    
        x = self.get_x_values(signal).values
        y = self.get_y_values(signal).values
        if x_limit:
            mask = (x < max(x_limit)) & (x > min(x_limit))
            x = x[mask]
            y = y[mask]
        if base_line is None:    
            base_line = y.min()
        mu_0, mu_1, mu_2, skew = self.curve_stat_moments(signal, x_limit, base_line)
        
        
        HETP_1 = mu_2 * L / mu_1**2
        
        y = y - base_line
        t_max = x[y.argmax()]
    
        
        FWHM_low = x[abs(y[x < t_max] - (y.max() / 2)).argmin()]
        FWHM_high= x[y.argmax() + abs(y[x > t_max] - (y.max() / 2)).argmin()]
    
        FWHM = FWHM_high - FWHM_low
        print('delta: ', FWHM,' t_max: ', t_max,' L ', L)
        HETP_2 = L / 5.54 * (FWHM / t_max)**2
        
        C_max = y.max() + base_line
        HETP_4 = 2 * np.pi * (mu_0 / (C_max * t_max) )**2 * L
        
        if plot == True:
            fig, ax = plt.subplots(figsize = (9, 5))
            trans = transforms.blended_transform_factory(
                    ax.transData, ax.transAxes)
            ax.plot(x,y)
            ax.vlines(t_max, 0, 1, transform = trans)
            ax.annotate('', 
                        xy = (FWHM_low, y.max() / 2), 
                        xytext = (FWHM_high, y.max() / 2), 
                        arrowprops = dict(arrowstyle='<->'))
            ax.text(t_max, y.max()/2, 'FWHM = {0:.2f}'.format(FWHM), ha = 'center', va = 'bottom')
            ax.vlines(mu_1, 0, 1, color = 'g', transform = trans)
            ax.axvspan(mu_1 - mu_2/2, mu_1 + mu_2/2, facecolor='g', alpha=0.3)      
            text ='HETP_stat: {0:.3f}\nHETP_Uni {1:.3f}\nHETP_4: {2:.3f}'.format(HETP_1, HETP_2, HETP_4)
            ax.text(0.1, 0.99, text, transform = ax.transAxes, va = 'top')
            text ='Plates/meter:\n\nStat: {0:.3f}\nUnicorn: {1:.3f}\nCetrti: {2:.3f}'
            text = text.format(100/HETP_1, 100/HETP_2, 100/HETP_4)
            ax.text(0.1, 0.5, text, transform = ax.transAxes, va = 'top')
            ax.set_xlabel('Volume [mL]')
            ax.set_ylabel(signal)
            fig.tight_layout()
        print('skew: ', skew)
        return HETP_1, HETP_2, HETP_4     

    
    def plot_signals(self, signals = None, fig = None, lab = ''):
        '''narise signale
        signal na drugem mestu da na desno os
        vsak signal ma svojo y skalo a vsi imajo isto x skalo'''
        
        all_available_sig = [name for name in self.curves.columns if not name.startswith('vol')]
        
        if signals is None:
            signals = all_available_sig
        
        if fig is None:
            fig = plt.figure()
        axes = {}
        lns = []
        
        ax_0 = fig.add_subplot(111)
        ax_current = ax_0
        wrong_flag = False
        
        
        for i, signal in enumerate(signals):
            
            if signal in self.curves.columns:
                if i > 0:
                    ax_current = ax_0.twinx()
                #set x and y data       
                x = self.curves['vol' + signal].dropna().values
                y = self.curves[signal].dropna().values
                
                #plot data on last axes
                current_color = colorAtNumber(i)
                ln = ax_current.plot(x, y, label = signal + '_' + lab, 
                                     linewidth = 3, color = current_color)
                lns.append(ln[0])
                
                if i == 0:
                    ax_current.set_ylabel(signal + ' [{0}]'.format(self.units[signal]))  
                elif i == 1:
                    ax_current.set_ylabel(signal + ' [{0}]'.format(self.units[signal]))
                else:
                    ax_current.set_yticks([])
                 
                    
                #add matplotlib axes (mnozina je za eno) to the axes dictionary
                axes[signal] = ax_current
                
            else:
                wrong_flag= True
                print(signal + ' is not in your data')
        
        if wrong_flag:
            print('You can select only: {0}'.format(all_available_sig))
        ax_0.set_xlabel('Volume [{0}]'.format(self.units[0]))
        labs = [l.get_label() for l in lns]
        ax_0.legend(lns, labs, loc='upper center')
        ax_0.set_title(self.name)
        return ax_0, fig
        
       
                
                
        
        
class ChromatogramsContainer(object):
    '''class to hold chromatogram objects defined above'''
    
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
                
            
    def get_chrom_names(self):
        '''Return list of chrom names'''
        return [el[0] for el in self.__chromatograms]
    
    def chromatograms(self):
        return self.__chromatograms
    
    def filename(self):
        return self.__fname        
        
    def add_chromatogram(self, chrom):
        '''add chromatogram object to ChromatogramsContainer'''
        
        if id(chrom) in self.__chromFromID:
            return False
        key = chrom.name
        bisect.insort_left(self.__chromatograms, [key, chrom])
        self.__chromFromID[id(chrom)] = chrom
        self.__dirty = True
        return True
    
    def remove_chromatogram(self, chrom):
        
        if id(chrom) not in self.__chromFromID:
            return False
        del self.__chromFromID[id(chrom)]
        name = chrom.name
        i = bisect.bisect_left(self.__chromatograms, [name, chrom])
        del self.__chromatograms[i]
        del chrom
        self.__dirty = True
        return True
        
    def clear(self):
        
        self.__chromFromID = {}
        self.__chromatograms = []
        
    def create_sim_chromatogram(self):
        '''Solves differential equation for PFR and creates data'''
        pass
    
    def load_chromatogram(self, fname):
        '''
        from Unicorn exported file
        '''
        
        error = None
        try:
            name = ntpath.basename(fname).split('.')[0]
            #read in the exported ascii, csv
            temp_df = pd.read_table(fname, skiprows=[0], low_memory=False,
             sep='\t',encoding= CODEC)
        
            #rename the column names (glupa struktura exporta je)
            #prefix = 'vol'
            #if units[0 == 'ml']
            col = np.array([ ('vol'+name, name) for name in temp_df.columns if 'Unnamed' not in name]).flatten()
            
            #set the columns of curves to new columns
            temp_df.columns = col
            
            #set the units of the columns
            units = temp_df.iloc[0]
            temp_df = temp_df.drop(temp_df.index[0])
            

            fractions = None
            #check if there are fractions, if there are set the fractions data
            if 'Fraction' in col:
                temp_df, fractions = self.get_fractions(temp_df)
            
            injections = None
            if 'Injection' in col:
                temp_df, injections = self.get_injections(temp_df)
                
            
            phases = None
            #check if there are phases defined (only in Unicorn 7...)
            if 'Run Log' in col:
                temp_df, phases = self.get_phases(temp_df)
            
            curves = temp_df.astype(np.float64)
            
            chrom_metadata = self.get_chrom_metadata(fname)
            
            loaded_chrom = Chromatogram(name, units, curves, chrom_metadata,
                                        fractions, phases, injections)
            
            self.add_chromatogram(loaded_chrom)
            
        except Exception as e:
            error = 'During the load_chromatogram we got an error: {}'.format(e)
        finally:
            if error is not None:
                return False, error
            
            msg = 'Loaded additional chromatogram: {0}'.format(fname)
            return True, msg

    def get_chrom_metadata(self, fname):
        '''look at the folder for chrom_add_info.txt, and get additional info about the run'''
        error = None
        fh = None
        
        try:

            fh = open(ntpath.dirname(fname)+ '/chrom_add_info.txt', 'r')
            for line in fh:
                items_list = line.split('\t')
                if items_list[0] == ntpath.splitext(ntpath.basename(fname))[0]:
                    return True, float(items_list[-1])
        except (IOError, OSError, ValueError) as e:
            error = 'Failed to load metadata: {0}'.format(e)
        finally:
            if fh is not None:
                fh.close()
            if error is not None:
                return False, error

            
    def get_fractions(self, temp_df):
        '''iz dataFrame dobim frakcije
        specificno za Unicorn 7'''
        
        fractions = temp_df['volFraction'].dropna().astype(np.float64)
        fraction_labels = temp_df['Fraction'].dropna().values
        
        temp_df = temp_df.drop('volFraction',axis = 1)
        temp_df = temp_df.drop('Fraction', axis = 1)
        
        fractions.index = fraction_labels
                                    
        return temp_df, fractions
    
    def get_injections(self, temp_df):
        '''iz dataFrame dobim injection_markse 
        specificno za Unicorn 7'''
    
        injection_marks = temp_df['volInjection'].dropna().astype(np.float64)
        injection_marks.index = injection_marks.index = 'Injection_' + \
                                        injection_marks.index.astype(str)
      
        temp_df = temp_df.drop('volInjection',axis = 1)
        temp_df = temp_df.drop('Injection', axis = 1)
                                  
        return temp_df, injection_marks
    
    def get_phases(self, temp_df):
        '''v Unicornu 7 vrne imena in pozicij zacetkov in konca faze'''
        
        phases = temp_df['volRun Log'].dropna().astype(np.float64)
        phase_labels = temp_df['Run Log'].dropna().values
        
        temp_df = temp_df.drop('volRun Log',axis = 1)
        temp_df = temp_df.drop('Run Log', axis = 1)
        
        phases.index = phase_labels
        
        return temp_df, phases

    def get_possible_alignments(self, intersection = False):
        '''returns all possible alignements, by name,
        if intersection, it returns only indices present in all loaded chromatograms
        '''

        idx = self.__chromatograms[0][1].all_alignments.index.unique()
        if len(self.__chromatograms) < 2:
            return idx
        
        if intersection:
            for name, chrom in self.__chromatograms[1:]:
                idx = idx.intersection(chrom.all_alignments.index.unique())        
        else:
            for name, chrom in self.__chromatograms[1:]:
                idx = idx.union(chrom.all_alignments.index.unique())
        return list(idx)
        
        
    def load_multiple_chrom(self):
        '''Loads all chromatograms from a folder'''
        pass

    @staticmethod
    def formats():
        return "*.csv *.asc"



if __name__ == '__main__':
    
    print('\n\n\n\n')
    folder = '../../../01_Akta_files/'
    c = ChromatogramsContainer()
    ok, msg = c.load_chromatogram(folder + 'LAG525 AEC VCS MuLV Run1 001.csv')
    print(msg)
    ok, msg = c.load_chromatogram(folder + 'LAG525-D024-17-ALC 001.csv')
    print(msg)
    ok, msg = c.load_chromatogram(r'../../../01_Akta_Files/LAG525 AEC VCS MuLV Run2 001.csv')
    print(msg)
    ok, msg = c.load_chromatogram(folder + 'Testi/Test kolone LAG525(42) MabSelect Sure 001.csv')
    print(len(c))
    chrom = c[-1]
    HETPs = chrom.HETP_calc('Cond', [0,20], 19.8, plot = False)
    text ='\nHETP_stat: {0:.3f}, HETP_Uni {1:.3f}, HETP_4: {2:.3f}'.format(*HETPs)
    print(text)
    #c.plot_signals(['UV 1_300', 'Cond'])
    limits = [0,20]
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    