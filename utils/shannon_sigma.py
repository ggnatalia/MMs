#! /usr/bin/python3

import numpy as np
import random
import matplotlib.pyplot as plt
import math
import statsmodels.formula.api as smf
import argparse
import pandas as pd
from collections import Counter
#from sklearn.preprocessing import PolynomialFeatures
#from sklearn.linear_model import LinearRegression

"""
------------------------------------------------------------------------------------------------------------------------------------
Plot relation between shannon diversity Index (H) and standard deviation from log normal
Run this script with default parameters.
usage: shannon_sigma_FINAL.py [-h] [-sd SDS] [-N NSPECIESSET]
                              [--totalCounts TOTALCOUNTS] [--figsize FIGSIZE]

Process input files

optional arguments:
  -h, --help            show this help message and exit
  -sd SDS, --sds SDS    Standard deviation start, end, step
  -N NSPECIESSET, --nSpeciesset NSPECIESSET
                        Number of different species
  --totalCounts TOTALCOUNTS
                        Number of counts of each sample approx
  --figsize FIGSIZE     Figure size
------------------------------------------------------------------------------------------------------------------------------------
"""


def parse_arguments():
    """Parse the command line arguments and return an object containing them"""
    # Required
    general = argparse.ArgumentParser(description = 'Process input files')
    general.add_argument('-sd','--sdP', type = float, nargs = '+', default = [0.1, 3, 0.1], help = 'Standard deviation Parameters: start, end, step')
    general.add_argument('-N','--nSpeciesset', type = int, nargs = '+', default = [20, 50, 100, 500, 1000, 5000, 10000], help = 'Number of different species')
    general.add_argument('--reads', default = 10000, type = int, help = 'Number of counts of each sample approximately')
    general.add_argument('--figsize', default = (18,10), type = tuple, help = 'Figure size')
    args = general.parse_args()
    return(args)

def lognormal(mu, sd, n):
    """Calculate distribution for n elements (e. g.: sequences)"""
    # log-normal distribution with specified mean, standard deviation, and array shape.
    # Note that the mean and standard deviation are not the values for the distribution itself,
    # but of the underlying normal distribution it is derived from.
    # When sigma is low: all the abundances are similar. When it is high, there is a sequence that prevails
    abunlognormal = np.random.lognormal(mu, sd, n)
    return(abunlognormal)

def shannonindex(sample):
    """Calculate shannon index of the sample"""
    #p = species/sum(abunlognormal)
    total_abun = sum([x[1] for x in sample])
    H = -sum([(abun/total_abun) * np.log(abun/total_abun) for seq, abun in sample])
    return(H)
        
def plotAbundances(sample):
    """ Plot the abundances of each species """
    index = [i for i in range(1,len([x[0] for x in sample])+1)]
    label = ['seq_{}'.format(i) for i in range(1,len([x[0] for x in sample])+1)]
    plt.bar(index, [x[1] for x in sample])
    plt.xlabel('Sequences', fontsize=10)
    plt.ylabel('Abundances', fontsize=10)
    plt.xticks(index, label, fontsize=10, rotation=30)
    plt.title('Log normal distribution')
    plt.show()

def main(args):
    #totalCounts = 10000
    #nSpeciesset = [20, 50, 100, 500, 1000, 5000]
    #sds = np.arange(0.1,  3 + 0.1, 0.1) # sigmas of the  normal distribution (logarithm of abundances)
    reads = args.reads
    sds = np.arange(args.sdP[0], args.sdP[1] + args.sdP[0], args.sdP[2])
    nSpeciesset = args.nSpeciesset
    figsize = args.figsize
    data_to_plot_data = {}
    data_to_plot_formula = {}
    for nSpecies in nSpeciesset:
    # Knowing reads per sample: totalCounts (10000) -> calculate media of reads for each sequence in a normal distribution:
        # E.g. 20 seqs: 10000/20 = 500 | do the np.log of the result
        #mu = np.log(reads/nSpecies)
        mu = 0
    # Names of nSpecies different
        seqs = ['seq_{}'.format(i) for i in range(1, nSpecies+1)]
        counts = {s:0 for s in seqs}
        # Shannon Index: Typical values are generally between 1.5 and 3.5 in most ecological studies, and the index is rarely greater than 4.
        H_values = []
        for sd in sds: 
          H_medio = []
          i=0
          while i < 1000: # replicates, as is a random process: higher sd, higher variations in data subsampled from lognormal. Need to make several samplings
            ln_subset = lognormal(mu, sd, nSpecies)
            ln_normalize = ln_subset/sum(ln_subset)
            subset_reads = np.random.choice(seqs, reads, p = ln_normalize)
            freq = Counter(subset_reads)
            for k,v in freq.items():
                counts[k] = v
            H_medio.append(shannonindex(list(zip(counts.keys(), counts.values()))))
            i = i+1
          H_values.append(np.mean(H_medio))
          # Formula: H = lnS -(1/2)*sd² 
        H_valuesFormula = -0.5*sds**2 + np.log(nSpecies) # This formula is based on sd of the normal distribution associated to the ln(variable)
        data_to_plot_data[nSpecies] = H_values
        data_to_plot_formula[nSpecies] = H_valuesFormula
    #Make the plot
    with open('HvsSd.stats.txt','w') as statsSumm:
        fig = plt.figure(figsize = figsize)
        for num in range(1,len(nSpeciesset)+1):
            x = sds
            y = data_to_plot_data[nSpeciesset[num-1]]
            #z = np.poly1d(np.polyfit(x,y,2)) # same result that smf.ols with the quadratic term
            dataset = pd.DataFrame({'x': x, 'y': y})
            model = smf.ols('y ~ x + I(x**2)', data = dataset)
            result = model.fit()
            statsSumm.write('{}\n'.format(result.summary()))
            z2 = result.params[2]*(x**2) + result.params[1]*x + result.params[0]
            ax = fig.add_subplot(2, 3, num)
            ax.scatter(x, y, color = 'black', s = 2)
            ax.plot(x, data_to_plot_formula[nSpeciesset[num-1]], color = 'blue')
            ax.set_title('S ={}'.format(nSpeciesset[num-1]))
            ax.set_ylabel('H', fontsize = 12)
            ax.set_xlabel('sd', fontsize = 12)
            ax.set_ylim(0)
            plt.plot(x, z2, color='red')
        plt.tight_layout()
        #plt.show()
        plt.savefig('HvsSd.png')
        plt.close()
    
# Additional info related to lognormal distribution, shannon Index
#https://reader.elsevier.com/reader/sd/pii/0022098171900190?token=BF1E2462D20D5F4337E2A715C928C2EC1A039678D93CB3AFD12925CB51216C9DF3463F5CBE4A32B234CDD88418BFCD21
#https://reader.elsevier.com/reader/sd/pii/0040580971900207?token=BBCAEFB363F5ADC280A07DDEDF8B48F89591BD07FD4F94E3638D8329D2D2E84DAC1C819DCDA9461809ACEC2185DA9AB1
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4643760/
#https://www.jstor.org/stable/2529621?seq=1#page_scan_tab_contents
#https://stackoverflow.com/questions/50949089/translating-log-normal-and-log-normal-truncated-simulations-in-risk-to-python
#https://www.sciencedirect.com/topics/engineering/lognormal-distribution
#https://es.wikipedia.org/wiki/Distribuci%C3%B3n_log-normal


#https://www.statsmodels.org/dev/generated/statsmodels.regression.linear_model.OLS.html
#https://es.coursera.org/lecture/regression-modeling-practice/python-lesson-3-polynomial-regression-E3pjy
################################################################################################################
  
if __name__ == '__main__':
    main(parse_arguments())






