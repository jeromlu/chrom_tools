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
from matplotlib.transforms import blended_transform_factory

print(os.getcwd(), 'first', __name__)
#import of my or modified packages
#sys.path.append(os.path.join(os.path.dirname(__file__), ''))
from chrom_tools.pyEnergyDiagram.energydiagram import ED
import chrom_tools

#*********GLOBAL VARIABLES
HEADER_ROW = 1
INDEX_COL = 2
#descrpitors for auto labeling of bars or xtick_labels, should not be sensitive to the case
ALLOWED_LABEL_SELECTORS = ('values', 'rel_diff', 'yield')

class AnalyticalData(object):
    '''
    Data object for handling of analytical data. It keeps data about your samples.
    The excel file has to have two sheets: "Eluates" and "Loads".
    Second data file has to have data about the refereence biologics.
    All the sheets have to have same columns.
    Required columns are: sampleID, eks, scout, loadID ....
    
    Instantiation:
        AnalyticalData(fname, orig_fname)
    
    Parameters
        ----------
        fname : str 
            path and filename of the excel file where analytical data is stored
        orig_fname : str
            path and filename of the excel file where analytical data of originator is held
        
        Returns
        -------
        fig (plt.figure) and ax (fig.add_subplot()) 

    
    '''
    def __init__(self, fname = None, orig_fname = None, filter_dict = None):
        
        
        self.loads = pd.DataFrame([])
        
        self.eluates = pd.DataFrame([])
        
        if orig_fname is not None:
            self.originator = OriginatorData(orig_fname)
        else:
            self.originator = None
        
        if fname is not None:
            self.fname = fname
            self.load_data(filter_dict)
            
            
            
    def load_data(self, filter_dict):
        '''Loads the excel file into the Analytical data data structure
        '''
        
        na_values = ['n.a.', 'na']
        zeros = ['BPQL', 'BDL', 'LOQ']
        header_row = 1
        index_col = 2
        
        try:
            self.eluates = pd.read_excel(self.fname,
                                       sheet_name = 'Eluates',
                                       header = header_row,
                                       index_col = index_col,
                                       na_values = na_values)
            self.eluates.replace(zeros, 0, inplace = True)
            if (filter_dict is not None) and (type(filter_dict) == dict):
                self.filter_eluate_samples(filter_dict)
                
            self.loads = pd.read_excel(self.fname,
                                       sheet_name = 'Loads',
                                       header = header_row,
                                       index_col = index_col,
                                       na_values = na_values)   
            if self.loads.index.is_unique == False:
                print(self.loads[self.loads.index.duplicated()])
                raise IndexError('Index of loads is not unique.')
            self.loads.replace(zeros, 0, inplace = True)
            #if self.eluates.index.is_unique == False:
            #    print(self.eluates.index[self.eluates.index.duplicated()])
            #    raise IndexError('Index of eluates is not unique.')
            return True    
        except Exception as e:
            self.print_err()
            return False
    
    def filter_eluate_samples(self, filter_dict):
        '''
        Filters the data according to the data filter
        
        Parameters
        -----------
        
        filter_dict: dictionary
            dictionary of filters, dictionar key muxt be in columns, otherwise
            nothing happens
        
        Returns
        -----------
        
        sel_data (pd.DataFrame)
        '''
        mask = np.ones(len(self.eluates), dtype = 'bool')
        for key in filter_dict:
            if key in self.eluates.columns:
                mask_key = np.zeros(len(self.eluates), dtype = 'bool')
                for value in filter_dict[key]:
                    new_mask_value = self.eluates[key] == value
                    mask_key = mask_key | new_mask_value
                mask = mask & mask_key
            else:
                print(key, ' not in column names.')
                print('You may select:')
                print(self.eluates.columns)
        
        self.eluates = self.eluates[mask]
        
    @staticmethod
    def autolabel(ax, rects, values = None, fs = 12):
        '''Function that attaches some text labels to the rectangles. 
        It tries to format the text in some sensible way.
        
        Parameters:
        -----------
        ax : mpl.axes
            axes of bar plot
        rects : mpl.rectangles
            bars to which to attache the value labels
        values : None or list like object
            if None the labels are height of bars. If listlike object then 
            the passed values are attached to the bars. 
            Number of elements in list should be the same as number of rectangles.
        fs : int
            font size of the label text.
        '''
        span = max(ax.get_ylim())-min(ax.get_ylim())
        bbox = {'facecolor':'white', 'alpha' : 0.85, 'edgecolor':'none',
                'boxstyle' :'round, pad=0.2'}
        
        if values is None:
            for rect in rects:
                height = rect.get_height()
                ax.text(rect.get_x()+rect.get_width()/2., height + 0.03 * span ,
                        '{0:.4g}'.format(height),
                        ha='center', va='bottom', fontsize = fs,
                        bbox = bbox)

        else:
            for rect, value in zip(rects,values):
                if type(value) is str:
                    text = value
                elif value:
                    text = '{0:0.4g}'.format(value)
                height = rect.get_height()
                ax.text(rect.get_x()+rect.get_width()/2., height + 0.03 * span, text,
                        ha='center', va='bottom', fontsize = fs,
                        bbox = bbox)
                
    def set_ax_size(self, wid, hei, ax=None, onlyprint = False):
        """ w, h: width, height in inches """
        if not ax: ax=plt.gca()
        l = ax.figure.subplotpars.left
        r = ax.figure.subplotpars.right
        t = ax.figure.subplotpars.top
        b = ax.figure.subplotpars.bottom
        fig = ax.figure
        w, h = fig.get_size_inches() * fig.dpi
        #print('w:{0:.3g} h: {1:.3g}'.format(w,h))
        subplots = np.array([l, r, t, b])
        sub_abs = np.array([l*w, r*w, t*h, b*h])
        #print('l {0:.2f} r {1:.2f} t {2:.2f} b {3:.2f}'.format(*subplots))
        #print('Absolut_values:\nl {0:.2f} r {1:.2f} t {2:.2f} b {3:.2f}\n\n'.format(*sub_abs))
        #print('ax_w: {0:0.2f}'.format(sub_abs[1] - sub_abs[0]))
        if onlyprint:
            return
        figw = float(wid)/(r-l)
        #print(figw)
        figh = float(hei)/(t-b)
        #print(figh)
        ax.figure.set_size_inches(figw, figh)
    
    @staticmethod    
    def create_labels(array, sig_places = 3):
        labels =['{0:.' + str(sig_places) + 'g}'.format(name) for name in array]
        return labels
        
        
    def __plot_QA(self, ax, QA, sel_data, fixed_ylim = None, show_orig = True,
                  sample_labels = True, annotate_orig = True, sep_samples = None,
                  x_ticks = 'description',
                  bar_labels = {'Eluates' : None, 'Loads' : None},
                  rotation = 90, fs = 12, load_id_column = 'loadID',
                  add_x_axes = None):
        '''
        Internal function that plots eluate and load data to the axis. 
        It is used by two functions
        
        Parameters
        ----------
        ax : mpl.axis
            axis to which the data should be plotted
        QA : str 
            the name of qualtity attribute to show
        sel_data : pandas.DataFrame
            the data where selected samples residue, it has to include column
            with name same as QA
        fixed_ylim : list (len = 2)
            fixed y limit of plot
        show_orig : bool
            toggles if we want to show originator three sigma range or not
        bar_labels : dict
            dictionary of pair {'Eluates' : val1, 'Loads' : val2}, if any of 
            val1 or val2 are None there will be no labels on the specified bars,
            othervise val1 or val2 are descriptive strings (columns) what label
            to put above the plot's bars, not sensitive to the case
        annotate_orig : bool
            defines if the originator value should be annotated
        x_ticks : str or None
            controls the xtick labels, "description" means it will write
            the sampel description on x tick labels. If None no ticks are shown.
        add_x_axes : None or list
            list of additional columns to be visible as x axis, 
            only first two elements of list are used
            
        Returns
        -------
        None
        '''
        
        samples = sel_data.index
        #width of a bar on the plot
        bar_width = 0.4
        bar_shift = bar_width / 1.5 
        
        
        #**********************PLOTTING*****************************************
        #set x tick positions
        x = np.arange(len(samples))
        ax.set_xticks(x)
        
        #set the eluates data
        y = np.nan_to_num(sel_data[QA].values)
        
        #set the appropriate load data
        load_y = np.nan_to_num(sel_data[load_id_column].map(self.loads[QA]))

        #creating plots
        rectsL = ax.bar(x, load_y, width = bar_width, label = 'Load')
        rectsE = ax.bar(x + bar_shift, y, width = bar_width, label = 'Eluate')
        
        #set x limits
        ax.set_xlim(-bar_width, len(x) - 1 + bar_shift + bar_width)
        ax.set_title(QA)
        
        #***********************create x tick labels****************************
        if x_ticks == 'description':
            x_labels = sel_data.description
        elif x_ticks == 'sample_names':
            x_labels = samples
        elif x_ticks is None:
            #this is very specific if you are fractionating your elution peak or more peaks
            if sep_samples is not None:
                tab1=['Frac_' + str(i) for i in range(0, sep_samples)]
                tab2=['Frac_' + str(i) for i in range(0, len(samples) - sep_samples)]
                x_labels = pd.Series(tab1 + tab2)
        if sample_labels and (x_ticks is not None):
            x_labels = x_labels + '\n' + sel_data.index
        
        #set x tick labels
        if x_ticks is None:
            ax.set_xticklabels('')
        else:
            ax.set_xticklabels(x_labels, rotation = rotation, 
                               ha = 'center', fontsize = fs)
            
        #**********************adding additional x axes************************
        if add_x_axes:
            counter = 0
            pos_x_ax = - 0.23
            for column in add_x_axes:
                if column in sel_data.columns:

                    
                    ax_new = ax.twiny()
                    
                    # Move twinned axis ticks and label from top to bottom
                    ax_new.xaxis.set_ticks_position("bottom")
                    ax_new.xaxis.set_label_position("bottom")
                    
                    # Offset the twin axis below the host
                    ax_new.spines["bottom"].set_position(("axes", pos_x_ax))
                    
                    # Turn on the frame for the twin axis, but then hide all 
                    # but the bottom spine
                    ax_new.set_frame_on(True)
                    ax_new.patch.set_visible(False)
                    for sp in ax_new.spines.values():
                        sp.set_visible(False)
                    ax_new.spines["bottom"].set_visible(True)
                    
                    #set same ticks positions
                    ax_new.set_xlim(-bar_width, len(x) - 1 + bar_shift + bar_width)
                    ax_new.set_xticks(x)
                    
                    #set the labels 
                    x_labels = self.create_labels(sel_data[column])
                    
                    ax_new.set_xticklabels(x_labels)
                    ax_new.set_xlabel(column)
                    
                    ax_new.xaxis.set_label_coords(1.01, pos_x_ax)
                    ax_new.grid(None)
                    
                    
                    #update the position and counter
                    pos_x_ax = pos_x_ax - 0.07
                    counter = counter + 1
                    
                    if counter == 3:
                        break
        
        #**********************ADDING ORIGINATOR********************************
        if self.originator.is_loaded():
            originator = self.originator.originator
        
        if QA in originator.index:
            orig_y = originator.average[QA]
        else:
            orig_y = np.nan

        trans = blended_transform_factory(ax.transAxes, ax.transData)
        
        if ~np.isnan(orig_y):
            ax.hlines(orig_y, x[0] - bar_width, len(x) + bar_width, color = 'green')
            if annotate_orig:
                ax.annotate('Orig\n{0:.2f}'.format(orig_y), xy=(1.01, orig_y), 
                            xycoords=trans, clip_on=False, va='center', 
                            fontSize=12, color='green')
            sigma = originator.sigma[QA]
    
            x = np.arange(x[0]-1,len(x)+1, 0.1)
            ax.fill_between(x, orig_y - 3*sigma, orig_y + 3*sigma, 
                            where =(x > x[0]) & (x < len(x)), 
                            facecolor='green', alpha=0.3,transform = ax.transData)
    
            lower = min(np.r_[y,load_y,orig_y,orig_y-3*sigma, orig_y+3*sigma])
            upper = max(np.r_[y,load_y,orig_y,orig_y-3*sigma, orig_y+3*sigma])
        else:
            lower = min(np.r_[y,load_y])
            upper = max(np.r_[y,load_y])

        #************SET UP y limits*******************************************
        #set the y plot limits, rel_margin tells you how much margin should be 
        #at the top and bottom of the plot
        rel_marg=0.1
        if fixed_ylim is None:
            ax.set_ylim(lower - (upper - lower) * rel_marg,
                        upper + (upper - lower) * rel_marg)
        else:
            ax.set_ylim(*fixed_ylim)
            
        trans = blended_transform_factory(ax.transData, ax.transAxes)
        #separate samples
        if sep_samples is not None:
            ax.vlines(np.array(sep_samples) - bar_width, 0, 1, transform = trans, lw = 2.5)    
        
        #************ adding labels to columns *********************************
        for key, selector in bar_labels.items():
            
            possible_selectors = np.r_[self.eluates.columns, ALLOWED_LABEL_SELECTORS]
            if (selector is not None) and (selector in possible_selectors):
                #create list of values that I would like to show on the graph
                values = None                            
                if selector.lower() == 'values':
                    values = None
                elif selector.lower() == 'rel_diff':
                    values = (y - load_y) / load_y * 100
                elif selector.lower() == 'yield':
                    values = sel_data['yield'] * 100
                else:
                    values = sel_data[selector]
  
                #put labels above load bars or eluates bar 
                #***********ELUATES********************
                if key.lower() == 'eluates':
                    self.autolabel(ax, rectsE, values, fs = fs)
                    
                #***********LOADS**********************
                elif key.lower() == 'loads':
                    self.autolabel(ax, rectsL, values, fs = fs)
            elif selector is None:
                pass
            else:
                print('You can give only following parameters (case insesitive):')
                print(ALLOWED_LABEL_SELECTORS)
                print('Or you may give name of a column in eluates data (case sensitive)')
        
            
    def plot_load_eluate_QA(self, QA, sel_data, fixed_ylim = None, save_lab = '', save = False,
                            sample_labels = False, show_orig = True, sep_samples = None, fs = 12,
                            x_ticks = 'description', rotation = 20,
                            bar_labels = {'Eluates' : 'values', 'Loads' : None},
                            load_id_column = 'loadID', add_x_axes = None):
        '''
        Plots selected samples (eluates) with its starting material (loads), 
        optionally adds the originator three sigma range and different variaties
        of labels (see Parameters)
        
        Parameters
        ----------
        QA : str 
            the name of qualtity attribute to show
        sel_data : pandas.DataFrame
            the data where selected samples residue, it has to include column
            with name same as QA
        fixed_ylim : list (len = 2)
            fixed y limit of plot
        save : bool
            if True the created figure is saved in a created subfolder TODO
        save_lab : str
            string that is added to the save file name
        show_orig : bool
            toggles if we want to show originator three sigma range or not
        sample_labels : bool
            toggles if we want to show elutaes labels on xtick labels
        sep_samples: int
            separates samples by vertical line, number of sep_samples are on 
            the right side of the line
        x_ticks : str or None
            controls the xtick labels, "description" means it will write
            the sampel description on x tick labels. If None no ticks are shown.
        bar_labels : dict
            dictionary of pair {'Eluates' : val1, 'Loads' : val2}, if any of 
            val1 or val2 are None there will be no labels on the specified bars,
            othervise val1 or val2 are descriptive strings (columns) what label
            to put above the plot's bars, not sensitive to the case
        add_x_axes : None or list
            list of additional columns to be visible as x axis, 
            only first two elements of list are used
        
        Returns
        -------
        fig (plt.figure) and ax (fig.add_subplot()) 

        '''
        samples = sel_data.index
        if (len(samples) < 2) or (len (samples) > 15):
            print("\n\nError:\nYou've passed inapropriate number of samples.\
                  \nYou should pass at least 2 but not more than 15 samples.\n")
            return
            
        #calculates the width of the figure appropriate to the number of samples
        #TODO: improve this
        wid = int(len(samples)) + 1
        #fig, ax = plt.subplots(figsize = (wid, 6))
        fig, ax = plt.subplots()
        #print(load_id_column)
        #print(sel_data[load_id_column])
        self.__plot_QA(ax, QA, sel_data, fixed_ylim = fixed_ylim,
                       bar_labels = bar_labels, x_ticks = x_ticks, 
                       sep_samples = sep_samples,
                       sample_labels = sample_labels,
                       fs = fs, rotation = rotation,
                       load_id_column = load_id_column,
                       add_x_axes = add_x_axes)
                           
        
        ax.legend(bbox_to_anchor = (0.98, 1), loc=2, borderaxespad=0.)
        
        #****************setup proper size of the plot**************************
        #TODO: imporve this, currently try and error way (try using fisxed size axes..)
        fig.tight_layout()
        scale = 0.5
        wid =  scale * (4 + int(len(samples)))
        self.set_ax_size(wid, 4, ax)
        fig.tight_layout()

        if save:
            fig.savefig('./03_slike/sample_groups/PCE_testing/load_eluate' + QA + save_lab  + '.jpg' )
        return fig, ax
         
    def plot_load_eluate_QAs(self, QAs, sel_data, fixed_ylim = None, save_lab = '', save = False,
                            sample_labels = False, show_orig = True, sep_samples = None,
                            x_ticks = 'description', rotation = 20,
                            bar_labels = {'Eluates' : 'values', 'Loads' : None}):
        '''
        Plots many QAs (load aluate bars) within one figure. It has similar
        options as the method which plots only one QA. Here is no legend.
        
        Parameters
        ----------
        QA : str 
            the name of qualtity attribute to show
        sel_data : pandas.DataFrame
            the data where selected samples residue, it has to include column
            with name same as QA
        fixed_ylim : list (len = 2)
            fixed y limit of plot
        save : bool
            if True the created figure is saved in a created subfolder TODO
        save_lab : str
            string that is added to the save file name
        show_orig : bool
            toggles if we want to show originator three sigma range or not
        sample_labels : bool
            toggles if we want to show elutaes labels on xtick labels
        sep_samples: int
            separates samples by vertical line, sep_samples are on 
            the right side of the line
        x_ticks : str
            controls the xtick labels, "description" means it will write
            the sampel description on x tick labels
        bar_labels : dict
            dictionary of pair {'Eluates' : val1, 'Loads' : val2}, if any of 
            val1 or val2 are None there will be no labels on the specified bars,
            othervise val1 or val2 are descriptive strings (columns) what label
            to put above the plot's bars, not sensitive to the case
        
        Returns
        -------
        fig (plt.figure) and ax (fig.add_subplot()) 
        '''
        samples = sel_data.index
        if (len(samples) < 2) or (len (samples) > 6):
            print("\n\nError:\nYou've passed inapropriate number of samples.\
                  \nYou should pass at least 2 but not more than 6 samples.\n")
            return
        n_QAs = len(QAs)
        wid = (len(samples) * n_QAs)
        if n_QAs <= 5:
            fig, ax_list = plt.subplots(1, 5, figsize = (wid, 4))
        elif n_QAs > 5 and n_QAs <=10:
            fig, ax_list = plt.subplots(2, 5, figsize = (wid, 8))
        elif n_QAs > 10 and n_QAs <=15:
            fig, ax_list = plt.subplots(3, 5, figsize = (wid, 12))
        else:
            print('You may not plot more than 15 QAs')
            return False
        
        last_row = n_QAs - 5
        for i, (ax, QA) in enumerate(zip(ax_list.ravel(), QAs)):            
            xticks = None
            if i > last_row:
                xticks = 'description'
            
            self.__plot_QA(ax, QA, sel_data, fixed_ylim = fixed_ylim,
                       bar_labels = bar_labels, annotate_orig = False,
                       x_ticks = xticks)
            
        fig.tight_layout()
        if save:
            fig.savefig('./03_slike/sample_groups/PCE_testing/load_eluate' + QA + save_lab  + '.jpg' )
        return fig, ax
    
    def plot_sequence_linked(self, seq_el_IDs, QA, step_names = None,
                             fs = 12, **kwargs):
        '''Plots linked results of passed sequence of samples, uses PyEnergyDiagrams. 
        Got idea from physics energy diagrams.
        TODO: one way arrow
        
        Parameters:
        -----------
        seq_el_ID : list like object
            list of sample IDs. List of strings, the strings should be in data index.
        QA : string
            quality attribute. Name of a numeric column in data (self).
        step_names : None or list
            if list of strings, these are used as the names of each step
        fs : int
            font size passed to diagram plot
        kwargs : dict
            key word arguments passed to diagram plot see: chrom_tools.py_energy_diagram
            
        Usage:
        ------
        
        >>> PIL3_seq = ['DC80','DC81', 'DC83','DC84', '1Z600']
        >>> step_names = ['ALC', 'VIN', 'CEX', 'MMC', 'UF/DF']
        >>>
        >>> data.plot_sequence_linked(PIL3_seq, 'SEC_HMWs', step_names)
        
        '''
        
        diagram = ED()
        
        if QA not in self.loads.columns:
            print(QA)
            print('Quality attribute not in the loads!!')
            print('You may select:')
            print(self.loads.columns)
            return
        
        if QA not in self.eluates.columns:
            print(QA)
            print('Quality attribute not in the eluates!!')
            print('You may select:')
            print(self.eluates.columns)
            return
       
        n_steps = len(seq_el_IDs)
        el_data = self.eluates.loc[seq_el_IDs]
        self.eluates['el_start_val'] = el_data.loadID.map(self.loads[QA])
        el_data = self.eluates.loc[seq_el_IDs]
        
        #create levels
        for i, ser_iterator in enumerate(el_data.iterrows()):
            idx, row = ser_iterator
            load = row['el_start_val']
            diagram.add_level(load, row.loadID)
            eluate = row[QA]
            diagram.add_level(eluate, idx, 'last')
            if i + 1 <= n_steps:
                if i + 1 < n_steps:
                    diagram.add_link(2*i+1, 2*i+2)
                diagram.add_arrow(2*i, 2*i+1)

        
        diagram.plot(show_IDs = False, fs = fs, **kwargs)
        
        if step_names is not None:
            for i, step in zip(set(diagram.positions), step_names):
                start = i * (diagram.dimension + diagram.space)
                pos = start + diagram.dimension / 2
                trans = blended_transform_factory(diagram.ax.transData, diagram.ax.transAxes)
                diagram.ax.text(pos, 1.1, step, ha = 'center', fontsize = fs + 3, transform = trans)
        
        return diagram
        
    def plot_sample_chrom(self, sampleID, signals = None, chrom_folder = '.', 
                          QA = None, align_to = None, **kwargs):
        '''
        Loads Unicorn exported chromatogram, and plots UV and Conductivity curve
        
        Parameters:
        -----------
        
        Returns:
        --------
        
        Usage:
        ------
        '''
        if sampleID in self.eluates.index:
            run_ID = self.eluates.loc[sampleID, 'eks']
            scout_num = self.eluates.loc[sampleID, 'scoutNum']
        else:
            print('{0} not in eluates analytical data'.format(sampleID))
            return
        
        #all names in given folder
        names = os.listdir(chrom_folder)
        
        csv_files = [el for el in names if 
                     os.path.isfile(os.path.join(chrom_folder, el)) and
                     (os.path.splitext(el)[-1] == '.csv')]
        
        if len(csv_files) == 0:
            print('No csv files to load in {0}.\nChange folder?'.format(chrom_folder))
        print('{0:03g}'.format(scout_num), run_ID)
        
        chrom_name = None
        for name in csv_files:
            if (run_ID in name) and ('{0:03g}'.format(scout_num) in name):
                chrom_name = name
                break

            
        if chrom_name:
            container = chrom_tools.ChromatogramsContainer()
            success, msg = container.load_chromatogram(os.path.join(chrom_folder, chrom_name))
            print(msg)
            if success:
                chrom = container[0]
                ax, fig = chrom.plot_signals(signals, mark_samples = True, fname = self.fname,
                                             align_to = align_to, sample_ID = sampleID, **kwargs)
                if QA:
                    QA_value = self.eluates.loc[sampleID, QA]
                    txt = QA + ' = ' + str(QA_value)
                    ax.text(0.6, 0.9, txt, transform = ax.transAxes, fontsize = 12)
            else:
                print(msg)
        else:
            print('No file with same name found!!! Run ID: {0}, scout num: {1:03g}'.format(run_ID, scout_num))
            return
                
            
        
        return chrom, fig, ax
                    
    def print_err(self):
            exc_type, exc_obj, exc_tb = sys.exc_info()           
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            err_msg = '{0}:\n{1}\nError occurred in file: {2}'.format(exc_type, exc_obj, fname)
            print(err_msg)
            logger.error(err_msg)     
            
            


class OriginatorData(object):
    '''
    Class for keeping the data about the reference biologics (originator)
    The class posseses the convenient methods for:
        - loading
    '''
    
    def __init__(self, fname = None, *args, **kwargs):
        
        self.fname = fname
        self.loaded = False
        
        self.originator = pd.DataFrame([])
        if fname is not None:
            self.load_data()

    
    def load_data(self):
        
        try:
            self.loaded = False
            raw_orig = pd.read_excel(self.fname, skiprows = [0])

            self.originator['average'] = raw_orig.mean()
            self.originator['sigma'] = raw_orig.std()
            self.originator['lower'] = self.originator.average - 3 * self.originator.sigma
            self.originator['upper'] = self.originator.average + 3 * self.originator.sigma
            self.loaded = True

        except Exception as e:
            self.loaded = False
            print(e)
            
    def is_loaded(self):
        return self.loaded
        

        


            
            
if __name__ == '__main__':
    
    folder = 'C:/Users/jeromlu2/LukaFiles/04_Programiranje/01_Python/03_IPythonNotebooks/temp/IPI/03_DSP Development/04_DSPD/'
    fname = '01_Akta_runs/GNT904_data_akta_runs.xlsx'
    orig_fname = '03_Originator/QTPP data table original.xlsx'
    chrom_folder = '01_Akta_runs/02_AktaFiles/CEX_optimization'
    
    #or logic between values, and logic between keys
    filter_dict = {'class' : ['optimization', 'development_seq_1'],
                   'step' : ['CEX', 'MMC']}
    
    data = AnalyticalData(folder + fname, folder + orig_fname)
    #diagram = data.plot_sequence_linked(['1Z918','1Z920'],'SEC_HMWs')
    #diagram2 = data.plot_sequence_linked(['1Z918'],'SEC_HMWs')
    
    mask = (data.eluates.eks == 'D243') | (data.eluates.eks == 'D242')
    eks_data = data.eluates[mask]
    
    bar_labels = {'Eluates' : 'values', 'Loads' : 'SEC_HMWs'}
    data.plot_load_eluate_QA('HCP', eks_data.iloc[:10], fixed_ylim=[0], bar_labels = bar_labels, 
                             add_x_axes = ['pH', 'concEl'])
    #data.plot_load_eluate_QA('HCP', eks_data.iloc[:15], fixed_ylim=[0], 
    #                         bar_labels = bar_labels)
    #data.plot_load_eluate_QA('HCP', eks_data.iloc[:8], fixed_ylim=[0], 
    #                         bar_labels = bar_labels)
    #data.plot_load_eluate_QA('HCP', eks_data.iloc[:6], fixed_ylim=[0], 
    #                         bar_labels = bar_labels)
    #fig, ax = data.plot_load_eluate_QA('HCP', eks_data.iloc[:4], fixed_ylim=[0], 
    #                                   bar_labels = bar_labels)
    #data.plot_load_eluate_QA('HCP', eks_data.iloc[:2], fixed_ylim=[0], 
    #                         bar_labels = bar_labels)
    
    #QAs = ['HCP', 'SEC_HMWs',  'nrCE_SDS_main', 'nrCE-SDS_HHL', 'iCE_sum_AP', 'iCE_sum_BP',
    #       'Afucosylation', 'High_Mannose', 'Galactosylation']
    #data.plot_load_eluate_QAs(QAs, eks_data.iloc[:3])
    #data.plot_load_eluate_QAs(QAs, eks_data.iloc[:4])
    #data.plot_load_eluate_QAs(QAs, eks_data.iloc[:5])
    
    #bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    #width, height = bbox.width, bbox.height
    
    #width *= fig.dpi
    #height *= fig.dpi
    
    #chrom = data.plot_sample_chrom('2A337', ['UV 1_280', 'Cond'], chrom_folder = folder + chrom_folder)
    
    '''
    others = ['SEC_HMWs', 'BP3', 'HCP']
    for CQA in others:
        for name, data in eks_data.groupby(['eks']):
            text = 'Eks {0}'.format(name)
            sep_samples = 6
            if 'D243' in name:
                sep_samples = 4
            fig, ax = plotAllSamplesOneCQA(CQA, data.sort_index(), userLab = text, save=True, 
                                           x_ticks = 'kr_neki', sep_samples = sep_samples, sample_labels = True)
            '''
                           