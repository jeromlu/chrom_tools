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

#my packages
import chrom_tools

CODEC = 'utf-16'

def color_at_number(n):
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
        
        #sampling data and analytical data
        self.sampling_data = {}
        
    def get_all_alignments(self):
        '''Concantinates all possible alignments, if present
        merges injections, phases, and fractions.
        
        Returns: None or concatinated pd.DataFrame
        '''        
        #print(pd.concat([self.fractions, self.phases ,self.injection_marks]))
        if (self.fractions is None) and (self.phases is None) and (self.injection_marks is None):
            return None
        else:
            return pd.concat([self.fractions, self.phases ,self.injection_marks], sort = True)

    
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
    
    @staticmethod
    def find_nearest_idx(array, value):
        '''Finds neareast index of value in array that is closest to the passed
        
        Parameters:
        -----------
        array : numpy.array
            array of values
        value : float
            value to be searched for in passed array
        
        Returns:
        --------
        idx : index of the found value
        
        Usage:
        ------
        >>> tab = np.arange(15)
        >>> find_nearest(tab, 6.6)
        >>> result = 7
        '''
        idx = (np.abs(array - value)).argmin()
        return idx
    
    def nearest_at_vol(self, signal, value):
        '''Returns index in column of selected signal that is closest to the given value
        '''
        x = self.get_x_values(signal)
        idx = self.find_nearest_idx(x, value)
        y = self.get_y_values(signal)
        return y[idx]
    
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
    
    def get_meta_data(self, fname = None):
        '''Reads the data about sampling and loads
        do this with analytical data class
        Parameters
        ----------
        fname : str or None
            file where analytical data are saved
        Returns
        -------
        None
        '''
        if fname:
            flag = self.sampling_data = chrom_tools.AnalyticalData(fname)
            return flag
        else:
            return False


    def mark_sampled_fractions(self, ax, fname, offset = 0, cv = 1, sample_ID = None):
        
        #initial settings
        myColors = ['green', 'red', 'blue', 'orange', 'grey']
        i = 0 # color counter
        trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)
        bbox={'facecolor':'white', 'alpha' :0.6, 'edgecolor':'none', 'boxstyle' :'round,pad=0.2'}
        
        if not self.get_meta_data(fname):
            print('No data about sampling information!')
            return
        data_elu = self.sampling_data.eluates
        data_loads = self.sampling_data.loads
        eks, scout = self.get_run_ID()
        mask = (data_elu.eks == eks) & (data_elu.scoutNum == int(scout))
        current_run = data_elu[mask]
        
        #create for coloring
        x_lim = list(ax.get_xlim())
        fr_min = (self.fractions.min() - offset) / cv
        fr_max = (self.fractions.max() - offset) / cv
        
        if fr_min > x_lim[0]:
            x_lim[0] = fr_min
        if fr_max < x_lim[1]:
            x_lim[1] = fr_max
        x = np.linspace(x_lim[0], x_lim[1], 1000)

        #x_positions of sample fractions and sample names
        names_X_pos = []
        sample_names = current_run.index
        
        fractions = self.fractions
        for s_name in sample_names:
            
            #if called with sample_ID paramater it marks only that sample
            if sample_ID and (s_name != sample_ID):
                print(sample_ID)
                continue
            upper_label = current_run.loc[s_name].fracUpper
            lower_label = current_run.loc[s_name].fracLower
            
            if (upper_label) is None or (lower_label is None):
                err_txt = 'For sample {0} there is no upper or lower limit!'
                print(err_txt.format(s_name))
                return
            if (upper_label not in fractions.index) or (lower_label not in fractions.index):
                err_txt = 'For sample {0} there is no upper or lower limit!'
                print(err_txt.format(s_name))
                print('Possible fractions are: ', fractions.index)
                return                
            
            upper = fractions.index.get_loc(upper_label) + 1
            lower = fractions.index.get_loc(lower_label)
            
            names_X_pos.append((fractions[upper]-fractions[lower])/2+ fractions[lower]-offset)
            ax.fill_between(x, 0, 0.7, 
                            where = ((fractions[lower]- offset)/cv < x) & \
                                    (x < (fractions[upper]- offset)/cv), 
                            facecolor=myColors[i], alpha=0.3, transform = trans)
            
            if i >= (len(myColors)-1):
                i = 0
            else:
                i=i+1
        
        for x_pos, s_name in zip(names_X_pos, sample_names):
            ax.text(x_pos / cv, 0.15, s_name,
                    ha = 'center', rotation = 'vertical', fontsize = 12,
                    transform=trans, zorder=5, bbox = bbox)

    def add_curve_to_axis(self, ax_0, signal, lab, cv = 1, offset = 0, 
                          twinx = False,  **kwargs):
        '''
        Adds curve to passed axis (matplotlib)
        Parameters
        ----------
        ax_0 : mpl.axes
        curve : list
            list of x, y values of the curve
        twinx : bool, default False
            if True new axes is created and returned
        cv : float
            column volume to normalize independant data (volume)
        offset : float
            translation of the curve
        **kwargs : plt.plot keyword options
            
        Returns
        -------
        line object and mpl.axes
        '''
        #set x and y data       
        x = self.get_x_values(signal, offset, cv)
        y = self.get_y_values(signal)
        
        if twinx:
            ax = ax_0.twinx()
        else:
            ax = ax_0
        
        ln, = ax.plot(x, y, label = signal + '_' + lab, 
                                     linewidth = 3, **kwargs)
        return ln, ax

    
    def plot_signals(self, signals = None, x_lim = None, fig = None, lab = '',
                     mark_samples = False, sample_ID = None, fname = None, 
                     align_to = None, cv = 1):
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
        if x_lim is not None:
            ax_0.set_xlim(x_lim)
        ax_current = ax_0
        wrong_flag = False
        
        if isinstance(align_to, str):
            if align_to in self.all_alignments.index:
                offset = self.all_alignments[align_to]
            else:
                offset = 0          
                print('You can only align to following:')
                print(self.all_alignments.index)
        elif isinstance(align_to, (int, float)):
            offset = align_to
        else:
            offset = 0
        
        
        for i, signal in enumerate(signals):
            
            if signal in self.curves.columns:
                if i > 0:
                    ax_current = ax_0.twinx()

                
                #plot data on last axes
                current_color = color_at_number(i)
                ln = self.add_curve_to_axis(ax_current, signal, lab, 
                                            color = current_color, offset = offset,
                                            cv = cv)
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
        if mark_samples:
            print(fname)
            self.mark_sampled_fractions(ax_0, fname, offset, cv, sample_ID)
        if wrong_flag:
            print('You can select only: {0}'.format(all_available_sig))
        ax_0.set_xlabel('Volume [{0}]'.format(self.units[0]))
        ax_0.set_zorder(20) #20 upam da je dovolj visoka stevilka
        ax_0.patch.set_alpha(0)
        labs = [l.get_label() for l in lns]
        l = ax_0.legend(lns, labs, loc='upper center')
        l.set_zorder(21)
        ax_0.set_title(self.name)
        return ax_0, fig
        
       
                
                
        
        
class ChromatogramsContainer(object):
    '''Class to hold chromatogram objects defined above.
    
    Most important method would be load_chromatogram, which loads additional chromatogram.
    
    It also defines method for convinient chromatogram selecting(by name or index).
    
    TODO: create chromatogram from simulated data.
    TODO: compare chromatograms.
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
        from Unicorn exported csv file, load the data
        
        Parameters
        ----------
        fname : str
            name of the file. It has to be full path with file extension at the end 
            or it can be relative path with file extansion at the end.
            
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
            name = ntpath.basename(fname).split('.')[0]
            if name in self:
                msg = 'Chromatogram {0} already loaded.'.format(name)
                raise ValueError(msg)
            #read in the exported ascii, csv
            temp_df = pd.read_table(fname, skiprows=[0], low_memory=False,
             sep='\t', encoding= CODEC)
        
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
            
            
            #convert to floating point numbers
            curves = temp_df.astype(np.float64)
            
            #get metadata, this should be changed
            chrom_metadata = self.get_chrom_metadata(fname)
            
            #instantiate the Chromatogram object with all parameters set
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
        
        
    def load_multiple_chroms(self, folder):
        '''
        Loads all chromatograms from a folder
        Parameters
        ----------
        folder : string
            path to the folder where Akta files should be
        Returns
        -------
        '''
        #all names in given folder
        names = os.listdir(folder)
        
        csv_files = [el for el in names if 
                     os.path.isfile(os.path.join(folder, el)) and
                     (os.path.splitext(el)[-1] == '.csv')]
        
        if len(csv_files) == 0:
            print('No csv files to load.')
        
        for name in csv_files:
            match = re.match('\w+-(\w+)-\w+-\w+\s+(\d\d\d).csv', name)
            if match:
                fname = os.path.join(folder, name)
                self.load_chromatogram(folder + fname)
                
    def compare_chroms(self, chrom_lst, signal, align_to = None, cv = False):
        '''
        Compare one signal of different chromatograms.
        
        Parameters
        ----------
        chrom_lst : list
            :list of integers or list of chrom names
        signal : str
            name of a signal, it should be present in at least one of the
            chromatograms
            
        Returns
        -------
        Axes of the comparison plot.
        mpl.axes
        '''
        
        if (type(chrom_lst) != list) and (not chrom_lst):
            msg = 'First argument should be list with chrom names.'
            print(msg)
            return None
        
        fig, ax_0 = plt.subplots(figsize = (10, 7))
        lns = []
        for name in chrom_lst:
            chrom = self[name]
            if signal in chrom.curves.columns:
                #offset, cv = chrom.get_offset(aign_to, ) 
                ln = chrom.add_curve_to_axis(ax_0, signal, lab = chrom.name, offset = 0, cv =1)
                lns.append(ln[0])
            else:
                print('Signal {0} not in the file {1}'.format(signal, chrom.name))
            
        labs = [l.get_label() for l in lns]
        l = ax_0.legend(lns, labs, loc = 'upper center')
        l.set_zorder(21)
            
            
        

    @staticmethod
    def formats():
        return "*.csv *.asc"



if __name__ == '__main__':
    import os
    print('\n\n\n\n')
    folder = './tests/01_Akta_files/'
    os.listdir(folder)
    c = ChromatogramsContainer()
    ok, msg = c.load_chromatogram(folder + 'GNT904A1-D243-18-CEX 001.csv')
    print(msg)
    ok, msg = c.load_chromatogram(folder + 'GNT904A1-D243-18-CEX 001.csv')
    print(msg)
    ok, msg = c.load_chromatogram(folder + 'GNT904A1-D242-18-CEX 001.csv')
    print(msg)
    ok, msg = c.load_chromatogram(folder + 'GNT904A1-D242-18-CEX 002.csv')
    print(msg)
    ok, msg = c.load_chromatogram(folder + 'GNT904A1-D243-18-CEX 002.csv')
    print(len(c))
    chrom = c[1]
    #HETPs = chrom.HETP_calc('Cond', [0,20], 19.8, plot = False)
    #text ='\nHETP_stat: {0:.3f}, HETP_Uni {1:.3f}, HETP_4: {2:.3f}'.format(*HETPs)
    #print(text)
    #c.plot_signals(['UV 1_300', 'Cond'])
    limits = [0,20]
    '''
    for chrom in c:
        signal = 'UV 1_280'
        if 'D243' in chrom.name:
            signal = 'UV 1_300'
        chrom.plot_signals([signal, 'Cond'], x_lim = [100, 400],  
                           fname = './tests/GNT904_data_akta_runs.xlsx', 
                           mark_samples=True)
    '''
    
    ok, msg = c.load_chromatogram(folder + 'GNT904A1-D269-18-CEX 001.csv')
    
    anal_fname = 'C:\\Users\\jeromlu2\\LukaFiles\\04_Programiranje\\01_Python\\03_IPythonNotebooks\\temp\\IPI\\03_DSP Development\\04_DSPD\\01_Akta_runs\\'
    chrom_2 = c['GNT904A1-D269-18-CEX 001']
    chrom_2.plot_signals(['UV 1_280', 'Cond'], x_lim = [0, 250],  
                                                       fname = anal_fname + 'GNT904_data_akta_runs.xlsx', 
                                                       mark_samples = True,
                                                       align_to = 'Elution')
    #c.load_multiple_chroms(folder)
    
    #c.compare_chroms([0,1,2, 3], 'UV 1_280')
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    