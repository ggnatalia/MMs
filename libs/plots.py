#!/usr/bin/python3

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import math
import scipy.stats as sp
import plotly.express as px
#import psutil kaleido?? #https://plotly.com/python/static-image-export/
#import plotly.io as pio

########################### PLOT FUNCTIONS
#def plot_heatmap(data, outputDir, title, vmin, vmax, center,  legendtitle = 'z', text = None, symmetric = None, figsize = (20, 20)):
    #""" Plot a heatmap """
    #palette = sns.diverging_palette(240, 10, n = 20, sep = 10, center = 'light') # n =21 [-1:1, 0.1]=> 21 breaks
    #mask = np.ones(data.shape, dtype = bool)
    #if symmetric: # plot triangular matrix
    #    mask[np.tril_indices_from(mask, k = 0)] = False
    #else:
    #    mask = False
    #plt.figure(figsize = figsize) # select figure size
    #with sns.axes_style("white"):
    #    if text: # add text values per cell
    #        p = sns.heatmap(data, annot = True, fmt = '.2f' , cmap = palette, square = True, vmin = float(vmin), vmax = float(vmax), center = float(center), mask = mask,  cbar_kws = {'label': legendtitle}) 
    #    else:
    #        p = sns.heatmap(data, annot = False, fmt = '.2f' , cmap = palette, square = True, vmin = float(vmin), vmax = float(vmax), center = float(center), mask = mask,  cbar_kws = {'label': legendtitle}) 
    #p.tick_params(right = False, bottom = True, labelright = False, labelbottom = True, top = False, labeltop = False, left = True, labelleft = True)
    #p.set_title(title, fontsize = 18)
    #plt.legend(fontsize = 'large', title_fontsize = '16')
    #plt.tight_layout()
    #plt.show()
    #plt.savefig('{}/{}.png'.format(outputDir, title))
    
    



def plot_heatmap(data, outputDir, title, zmin, zmax, legendtitle = 'z', symmetric = None):
    """ Plot a heatmap """
    #palette = sns.diverging_palette(240, 10, n = 20, sep = 10, center = 'light') # n =21 [-1:1, 0.1]=> 21 breaks
    mask = np.ones(data.shape, dtype = bool)
    if zmin == 0 and zmax == 1:
         color_continuous_scale = ['white', 'blue']
    else:
        color_continuous_scale = ['red', 'white' , 'green']
    if symmetric: # plot triangular matrix
        mask[np.tril_indices_from(mask, k = 0)] = False
    else:
        mask = False
    fig = px.imshow(data, x = data.index, y = data.columns, color_continuous_scale = color_continuous_scale, zmin = zmin, zmax = zmax , labels = dict(x = xlab, y = ylab, color = legendtitle))
    fig.update_xaxes(side = "bottom")
    fig.update_layout(title = {'text': title, 'xanchor': 'center', 'yanchor': 'top', 'y' : 1, 'x' : 0.5})
    #pio.write_image(fig, '{}/{}.svg'.format(outputDir, title))
    fig.write_html('{}/{}.html'.format(outputDir, title))
    fig.show()
    
    
    
    
    
    
    
    
 
    
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


#def barplotpc(df, title, outputDir, figsize = (20, 20), ylab = 'Axis y', xlab = 'Axis x', textSize = 8):
    #stacked_data = df.apply(lambda x: x*100/sum(x), axis=1)
    # plot
    #stacked_data.plot(kind="bar", stacked=True).legend(loc='center left', bbox_to_anchor=(1, 0.5))
    #plt.title(title)
    #plt.xlabel(xlab)
    #plt.ylabel(ylab)
    #plt.tight_layout()
    #plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.4)
    #plt.show()
    #plt.savefig('{}/{}.png'.format(outputDir, title))
    #plt.close()  
      
   
      

def barplotpc(df, title, outputDir, ylab = 'Axis y', xlab = 'Axis x'):
    stacked_data = df.apply(lambda x: x*100/sum(x), axis=1)
    # Sort by value
    df.sort_values(by = 'counts', axis = 1, ascending = False, inplace = True)
    # plot
    fig = px.bar(x = stacked_data.columns, y = stacked_data.iloc[0])
    fig.update_xaxes(showgrid = True, ticks = "outside")
    fig.update_layout(title = {'text': title, 'xanchor': 'center', 'yanchor': 'top', 'y' : 1, 'x' : 0.5}, xaxis_title = xlab, yaxis_title = ylab)
    #fig.write_image('{}/{}.png'.format(outputDir, title))
    #pio.write_image(fig, '{}/{}.svg'.format(outputDir, title))
    fig.write_html('{}/{}.html'.format(outputDir, title))
    fig.show()






