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