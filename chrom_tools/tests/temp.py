# -*- coding: utf-8 -*-
"""
Created on Mon Sep 10 18:19:11 2018

@author: JEROMLU2
"""

def markSampledFractions(ax, curves, eks, scout, dataElu, cv, offset, fromExcel = True):
    '''Doda oznake vzorcev in obarva, kako smo frakcionirali, prebere iz excela, ta more biti pravilno napisan'''
    #initial settings
    trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)
    
    myColors = ['green', 'red', 'blue', 'orange', 'grey']
    i=0
    

    
    mask = (dataElu.eks == eks) & (dataElu.scoutNum == int(scout))
    currentRun = dataElu[mask]
    fracLower = currentRun.fracLower.dropna().values
    fracUpper = currentRun.fracUpper.dropna().values

    
    fractions, fractionLabels = get_fractions(curves) 
    #print(fracLower, fracUpper)
    x = np.arange(ax.get_xlim()[0],ax.get_xlim()[1], 0.01)
    

    for spLab,zgLab in zip(fracLower, fracUpper):     
        


        #print(zgLab, zg, spLab, sp)


        imenaX.append((fractions[zg]-fractions[sp])/2+ fractions[sp]-offset)
        ax.fill_between(x, 0, 0.7, where =((fractions[sp]-offset)/cv < x) & (x < (fractions[zg]-offset)/cv), 
                        facecolor=myColors[i], alpha=0.3, transform = trans)

        if i >= (len(myColors)-1):
            i = 0
        else:
            i=i+1

    #napisi zdruzenih frakcij
    imena = currentRun.sampleID.values
    
    bbox={'facecolor':'white', 'alpha' :0.6, 'edgecolor':'none', 'boxstyle' :'round,pad=0.2'}
    
    for x,s in zip(imenaX,imena):
        ax.text(x/cv, 0.15, s,ha = 'center', rotation = 'vertical', fontsize =12, transform=trans, zorder=5, bbox = bbox)

    diff_height = 0.07
    i = 1
    conc = currentRun.concEl.values
    HMWs = currentRun.SEC_HMWs.values
    HCPs = currentRun.HCP.values
    for i, x in enumerate(imenaX):
        if (i % 2) == 0:
            rel_height = 0.92
        else:
            rel_height = 0.72
        s = 'Conc: {0:.1f}'.format(conc[i])
        ax.text(x/cv, rel_height, s,ha = 'center', fontsize =12, transform=trans,zorder=5, bbox = bbox)
        s = 'HMWs: {0:.1f}'.format(HMWs[i])
        ax.text(x/cv, rel_height - diff_height, s,ha = 'center', fontsize =12, transform=trans,zorder=5, bbox = bbox)
        s = 'HCP: {0:.0f}'.format(HCPs[i])
        ax.text(x/cv, rel_height - 2*diff_height, s, ha = 'center', fontsize =12, transform=trans,zorder=5, bbox = bbox)
        

#    ax.fill_between(X1, 0.86, 1.05, where =(frakcije[13]/cv > X1) & (X1 > frakcije[1]/cv), 
#                    facecolor='green', alpha=0.3, transform=trans)
#    ax.text(((frakcije[13]-frakcije[1])/2 + frakcije[1])/cv,3000, imena[-1],fontsize=12, ha = 'center')
        
        
        
def plotAllSamplesOneCQA(CQA, selectedDataElu, userLab = '', fixed_ylim = None, save = False, sample_labels = False,
                         show_orig = True, sep_samples = None, x_ticks = 'description', **kwargs):
    '''specificna funkcija narise vse vzorce v selectedData dataframu za izbrani CQA'''
    
    #izracuna sirino gldee na stevilo vzorcev
    samples = selectedDataElu.sampleID.values
    wid = int(len(samples))+1
    fig = plt.figure(figsize=(wid,6))
    ax = fig.add_subplot(111)

    
    barWidth = 0.4


    y = np.nan_to_num(selectedDataElu[CQA].values)
    x = np.arange(len(y))
    loadY = np.nan_to_num(selectedDataElu.loadID_new.map(dataLoad[CQA]))

    rectsL = ax.bar(x, loadY, width = barWidth, label = 'Load')
    rectsE = ax.bar(x+barWidth/1.5, y, width = barWidth, label = 'Eluate')
    
    yieldY = selectedDataElu['yield'].values

    #loadDescr = selectedDataElu.loadID.map(dataLoad.description)
    #xTicks = loadDescr.value_counts().cumsum()-loadDescr.value_counts().iloc[0] +loadDescr.value_counts()/2
    ax.set_xticks(x)
    
    #create x tick labels
    if x_ticks == 'description':
        x_labels = selectedDataElu.description
    else:
        if sep_samples is not None:
            tab1=['Frac_' + str(i) for i in range(0,sep_samples)]
            tab2=['Frac_' + str(i) for i in range(0, len(y)-sep_samples)]
            x_labels = pd.Series(tab1 + tab2)
    if sample_labels:
        x_labels = x_labels + '_' + selectedDataElu.index
    #ax.set_xticklabels(selectedDataElu.index, rotation=90, ha = 'center', fontsize = 12)
    #ax.set_xticklabels(selectedDataElu.Peak_cut_end.astype(str)+'\n'+selectedDataElu.index, rotation=90, ha = 'center', fontsize = 12)
    ax.set_xticklabels(x_labels, rotation = 20, ha = 'right', fontsize = 12)
    ax.set_xlim(-barWidth,len(x))


    trans = transforms.blended_transform_factory(ax.transAxes,ax.transData)

    if CQA in originator.index:
        origY = originator.average[CQA]
    else:
        origY = np.nan
    
    if ~np.isnan(origY):
        ax.hlines(origY,x[0]-barWidth,len(x)+barWidth,color = 'green')
        ax.annotate('Orig\n{0:.2f}'.format(origY), xy=(1.01, origY), 
                    xycoords=trans, clip_on=False, va='center', fontSize=12, color='green')
        sigma = originator.sigma[CQA]

        x = np.arange(x[0]-1,len(x)+1,0.1)
        ax.fill_between(x, origY-3*sigma, origY+3*sigma, where =(x > x[0]) & (x < len(x)), 
                        facecolor='green', alpha=0.3,transform = ax.transData)

        lowerVal = min(np.r_[y,loadY,origY,origY-3*sigma, origY+3*sigma])
        upperVal = max(np.r_[y,loadY,origY,origY-3*sigma, origY+3*sigma])
    else:
        lowerVal = min(np.r_[y,loadY])
        upperVal = max(np.r_[y,loadY])
        

    relMarg=0.1
    if fixed_ylim is None:
        ax.set_ylim(lowerVal - (upperVal - lowerVal) * relMarg,upperVal + (upperVal - lowerVal) * relMarg)
    else:
        ax.set_ylim(fixed_ylim[0], fixed_ylim[-1])
    ax.set_title(CQA)
    
    #************ adding labels to columns *********************************
    
    #***********ELUATES
    #val = selectedDataElu['yield']*100 #sek selectedDataElu.USPmat
    #val = selectedDataElu['Analiza']
    val = (y - loadY) / loadY * 100
    autolabel(ax, rectsE, values = val)
    
    #***********LOADS
    #val = 'Pr_'+selectedDataElu.USP_pro.apply(str) + '_' + selectedDataElu.USP_rep.apply(str)
    #val = selectedDataElu.loadID_new.map(dataLoad[CQA])
    #autolabel(ax,rectsL,values = val)
    
    trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)
    #separate samples
    if sep_samples is not None:
        ax.vlines(sep_samples - barWidth, 0, 1, transform = trans, lw = 2.5)
    
    ax.set_ylabel('Relative amount [%]')
    ax.legend(bbox_to_anchor = (0.98, 1), loc=2, borderaxespad=0.)#loc='best')
    fig.subplots_adjust(0.1, 0.18, 0.9 , 0.96 )
    #fig.tight_layout()
    plt.show()
    if save:
        fig.savefig('./03_slike/sample_groups/PCE_testing/load_eluate' + CQA + userLab  + '.jpg' )
    return fig, ax



#funkcija da doda oznake nad stolpce, imam dve opciji (yield ali pa visino stolpca)
def autolabel(ax,rects, values = False):
    # attach some text labels
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
  
    

def plotAllSamplesCQAs(CQAs, selectedDataElu, userLab = '', fixed_ylim = None):
    '''specificna funkcija narise vse vzorce v selectedData dataframu za izbrane CQAje'''
    
    #izracuna sirino gldee na stevilo vzorcev
    samples = selectedDataElu.sampleID.values
    wid = int(len(samples) * 5)+1
    barWidth = 0.4
    
    n_CQAs = len(CQAs)
    if n_CQAs > 5 and n_CQAs <=10:
        fig, ax_list = plt.subplots(2,5,figsize = (wid,8))
    elif n_CQAs <= 5:
        fig, ax_list = plt.subplots(1,5,figsize = (wid,4))
    elif n_CQAs > 10:
        fig, ax_list = plt.subplots(3,5,figsize = (wid,12))
    else:
        print('Prevec CQAjev hoces naenkrat narisat.')
        return False
    
    for i, (ax, CQA) in enumerate(zip(ax_list.ravel(),CQAs)):
        y = np.nan_to_num(selectedDataElu[CQA].values)
        x = np.arange(len(y))
        loadY = np.nan_to_num(selectedDataElu.loadID_MMC.map(dataLoad[CQA]))

        rectsL = ax.bar(x, loadY, width = barWidth, label = 'Load')
        rectsE = ax.bar(x+barWidth/1.5, y, width = barWidth, label = 'Eluate')
    
        yieldY = selectedDataElu['yield'].values

    #loadDescr = selectedDataElu.loadID.map(dataLoad.description)
    #xTicks = loadDescr.value_counts().cumsum()-loadDescr.value_counts().iloc[0] +loadDescr.value_counts()/2
        last_row = n_CQAs - 5
        if i >= last_row:
            ax.set_xticks(x)
            if i == last_row:
                x_label = selectedDataElu.index
            elif i == last_row + 1:
                x_label = selectedDataElu.load.round().apply(str)
            elif i == last_row +2:
                x_label = selectedDataElu.peakCutEnd.apply(str)
            else:
                x_label = selectedDataElu.index
            ax.set_xticklabels(x_label, rotation=90, ha = 'center', fontsize = 12)
        else:
            ax.set_xticklabels('')
    #ax.set_xticklabels(chrom, rotation=90, ha = 'center', fontsize = 12)
        
        ax.set_xlim(-barWidth,len(x))


        trans = transforms.blended_transform_factory(ax.transAxes,ax.transData)

        if CQA in originator.index:
            origY = originator.average[CQA]
        else:
            origY = np.nan

        if ~np.isnan(origY):
            ax.hlines(origY,x[0]-barWidth,len(x)+barWidth,color = 'green')
            #ax.annotate('Orig\n{0:.0f}'.format(origY), xy=(1.01, origY), 
            #            xycoords=trans, clip_on=False, va='center', fontSize=12, color='green')
            sigma = originator.sigma[CQA]

            x = np.arange(x[0]-1,len(x)+1,0.1)
            ax.fill_between(x, origY-3*sigma, origY+3*sigma, where =(x > x[0]) & (x < len(x)), 
                            facecolor='green', alpha=0.3,transform = ax.transData)

            lowerVal = min(np.r_[y,loadY,origY,origY-3*sigma, origY+3*sigma])
            upperVal = max(np.r_[y,loadY,origY,origY-3*sigma, origY+3*sigma])
        else:
            lowerVal = min(np.r_[y,loadY])
            upperVal = max(np.r_[y,loadY])


        relMarg=0.1
        if fixed_ylim is None:
            ax.set_ylim(lowerVal - (upperVal - lowerVal) * relMarg,upperVal + (upperVal - lowerVal) * relMarg)
        else:
            ax.set_ylim(fixed_ylim[0], fixed_ylim[-1])
        ax.set_title(CQA,fontsize = 17)
    val = selectedDataElu['yield']*100 #sek selectedDataElu.USPmat
        #val = selectedDataElu['Analiza']
    autolabel(ax,rectsE,values = val)
        #val = 'Pr_'+selectedDataElu.USP_pro.apply(str) + '_' + selectedDataElu.USP_rep.apply(str)
        #val = selectedDataElu.loadID_MMC.map(dataLoad['Analiza'])
        #autolabel(ax,rectsL,values = val)
    #ax_list[0].set_ylabel('Relative amount [%]')
    #ax_list[0,0].legend(bbox_to_anchor=(0.98, 1), loc=2, borderaxespad=0.)#loc='best')
    #fig.subplots_adjust(0.06,0.3,0.94,0.91,wspace = 0,hspace = 0)
    fig.tight_layout()
    plt.show()
    return fig    