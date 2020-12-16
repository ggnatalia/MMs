#!/usr/bin/env python3

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import math
import scipy.stats as sp

########################### PLOT FUNCTIONS
def plot_heatmap(data, outputDir, title, vmin, vmax, center,  legendtitle = 'z', text = None, symmetric = None, figsize = (20, 20)):
    """ Plot a heatmap """
    palette = sns.diverging_palette(240, 10, n = 20, sep = 10, center = 'light') # n =21 [-1:1, 0.1]=> 21 breaks
    mask = np.ones(data.shape, dtype = bool)
    if symmetric: # plot triangular matrix
        mask[np.tril_indices_from(mask, k = 0)] = False
    else:
        mask = False
    plt.figure(figsize = figsize) # select figure size
    with sns.axes_style("white"):
        if text: # add text values per cell
            p = sns.heatmap(data, annot = True, fmt = '.2f' , cmap = palette, square = True, vmin = float(vmin), vmax = float(vmax), center = float(center), mask = mask,  cbar_kws = {'label': legendtitle}) 
        else:
            p = sns.heatmap(data, annot = False, fmt = '.2f' , cmap = palette, square = True, vmin = float(vmin), vmax = float(vmax), center = float(center), mask = mask,  cbar_kws = {'label': legendtitle}) 
    p.tick_params(right = False, bottom = True, labelright = False, labelbottom = True, top = False, labeltop = False, left = True, labelleft = True)
    p.set_title(title, fontsize = 18)
    plt.legend(fontsize = 'large', title_fontsize = '16')
    plt.tight_layout()
    #plt.show()
    plt.savefig('{}/{}.png'.format(outputDir, title))
 
    
def barplot(df, title, outputDir, rowfig = None, colfig = None, figsize = (20, 20), T = False, ylab = 'Axis y', xlab = 'Axis x', textSize = 8):
    " Barplot by rows per default. If you want to plot it by columns, use T (transpose) option."
    #fig = plt.figure()
    fig = plt.figure(figsize = figsize)
    if T:
        df = df.T
    if not rowfig:
        rowfig = math.floor(math.sqrt(len(df.index))) #sqrt(nSpecies)
    if not colfig:
        if len(df.index)%rowfig == 0:
            colfig = math.floor(len(df.index)/rowfig)  #Dividendo = Divisor * Cociente + Resto -> Cociente = (Dividendo - R)/Divisor; colfig = (nSpecies - nSpecies%sqrt(nSpecies))/sqrt(nSpecies)
        else:
            colfig = math.floor(len(df.index)/rowfig + 1 )
    for row in range(len(df.index)):
        xlabels = list(df.columns.values)
        x = range(len(xlabels))
        y = df.iloc[row]
        ax = fig.add_subplot(rowfig, colfig, row + 1)
        ax.bar(x, y, color = 'blue')
        ax.set_title('{}'.format(df.index.values[row]), fontsize = textSize + 2)
        ax.set_ylabel(ylab, fontsize = textSize)
        ax.set_xlabel(xlab, fontsize = textSize)
        ax.set_xticks(np.arange(len(xlabels)))
        ax.set_xticklabels(xlabels)
        ax.tick_params(axis="x", labelsize = textSize - 2,rotation=70 )
    plt.tight_layout()
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.4)
    #plt.show()
    plt.savefig('{}/{}.png'.format(outputDir, title))
    plt.close()


def barplotpc(df, title, outputDir, figsize = (20, 20), ylab = 'Axis y', xlab = 'Axis x', textSize = 8):
    stacked_data = df.apply(lambda x: x*100/sum(x), axis=1)
    # plot
    stacked_data.plot(kind="bar", stacked=True).legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.title(title)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.tight_layout()
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.4)
    #plt.show()
    plt.savefig('{}/{}.png'.format(outputDir, title))
    plt.close()  
      
   
      
      
def plot_entropy(array, title, outputDir ):
    nrow, ncol = np.shape(array)
    # SIMPLY CODE:
    entropies = {}
    for pos in range(ncol):
        entropies[pos] = sp.entropy(list(map(float,{nt:i for nt, i in zip(*np.unique(array[:,pos], return_counts=True)) if nt != '.'}.values())))
    # entropies = {pos: sp.entropy(list(map(float,{nt:i for nt, i in zip(*np.unique(array[:,pos], return_counts=True)) if nt != '.'}.values()))) for pos in range(ncol)}
    # Impossible to see 50000 positions, need to remove some
    higherentropies = { pos:S for pos, S in entropies.items() if S > 0.0 }
    plotData = sorted(tuple(higherentropies.items()), key = lambda x: x[0])
    plotDataDF = pd.DataFrame(plotData, columns = ['Pos', 'Entropy']) # Add 1 to position, because python is 0 based
    plotDataDF = plotDataDF.set_index('Pos')
    barplot(plotDataDF, outputDir = outputDir, title = '{}.entropy'.format(title), figsize = (12,12), T = True, ylab = 'Entropy', xlab = 'Position', textSize = 8 )













