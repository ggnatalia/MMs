#!/usr/bin/python3
import random
import pandas as pd
import numpy as np
from scipy.stats import norm
from scipy.stats import nbinom
import scipy.stats as sp


########################### MATH RELATED FUNCTIONS
def calculate_distance(s0, s1):
    """ Given 2 aligned sequences, count the number fo differences in the sequences and divide by the maximun length of both.
        No take into account the '.' or the '-'.
        Count as difference '-/A' but no './A' """
    # Use maximun length of both sequences to be more restrictive with the differences I finally observe
    lengthmax = len(max(s0.replace('.','').replace('-',''), s1.replace('.','').replace('-','')))
    #print(str(lengthmax))
    #print(str(sum([1 for nt1, nt2 in zip(string1, string2) if nt1 != nt2 and nt1 != '.' and nt2 != '.'])))
    return(sum([nt0 != nt1 for nt0, nt1 in zip(s0, s1) if  nt0 != '.' and nt1 != '.'])/lengthmax)


def calculate_distance_set(a, b):
    """ s0 y s1 are dictionaries saving nt positions """
    s0 = a.copy()
    s1 = b.copy()
    dot_index = s1['.']|s0['.'] #to get any position with a '.'. We just can compare align positions: '.' vs. A is NOT an align position
    # Remove index of '.' positions in all nucleotides
    diff = set()
    for nt in ['A','T','C','G','-']:
        s0[nt] = s0[nt]-dot_index
        s1[nt] = s1[nt]-dot_index
        diff.add( len(s0[nt] - s1[nt] )) #a = set([1, 2, 3]), b = set([2, 3, 4]), a.symmetric_difference(b): {1, 4} no puedo hacer esto, porque en el fondo estar?a contando dos veces la misma diferencia
    #print(diff)
    lengthmax = max(len(s0['A']|s0['T']|s0['C']|s0['G']) , len(s1['A']|s1['T']|s1['C']|s1['G']))
    #print(lengthmax)
    # Count diff except in dot positions * choose '.' vs '.' or any '.'
    return(sum(diff)/lengthmax)


def estimate_abundances(Nfeatures, total = 100): 
    """ Create a random vector with percentages """
    abun = random.sample(range(1, total + 1), Nfeatures)
    # NOTE: https://stackoverflow.com/questions/8064629/random-numbers-that-add-to-100-matlab/8068956#8068956  ### FUTURE: No best method, after renormalizing, numbers are not uniformally distributed
    normalizedAbun = [round((i*total/sum(abun))) for i in abun]
    if sum(normalizedAbun) == total:
        return(normalizedAbun)
    else: # Si no suma el total o se pasa, quitarle/add la diferencia al primer elemento.
        diff = total-sum(normalizedAbun) #que pasa si se queda negativo el elemento [0] al hacerlo
        #print(str(diff))
        #print(str(normalizedAbun[0]))
        normalizedAbun[0] = normalizedAbun[0] + diff
        return(normalizedAbun)


def sample_from_multivariate_normal(CorrMatrix, nSamples, mu = 0):
    """ Create a multivariate normal distribution, mu = 0. Correlation matrix CorrMatrix. Sampling nSamples """
    # Number of variables, should be equal to the number of variables present in the correlation matrix. CorrMatrix.shape[0]
    multiNormal = np.random.multivariate_normal([mu]*CorrMatrix.shape[0], CorrMatrix, nSamples)
    return(multiNormal)

def probabilities_from_normal(array):
    """ Estimate the area-probability of a/given matrix of Zscores """
    if isinstance(array, float):
        print('Working with single Zscores. Not a matrix')
        probNormal = norm.cdf(array,  loc = 0, scale = 1)
    else:
        probNormal = np.zeros(array.shape , dtype = float)
        for row in range(array.shape[0]): # each sample
            if isinstance(array[row], float):
                print('Working with array of 1-dim Zscores. Not a matrix')
                probNormal[row] = norm.cdf(array[row],  loc = 0, scale = 1)
            else:
                for col in range(len(array[row])):
                    probNormal[row][col] = norm.cdf(array[row][col],  loc = 0, scale = 1)
    return(probNormal)

def qzinegbin(p, size, pstr0, prob = None, munb = None, nVariables = None):
    """   Percent point function of a Zero Inflated Negative Binomial Distribution   """ # Same nomenclature that R function.
    # 1. Requirements:
    # Size, munb, prob and pstr0 must have the same length than p.
    # Given that each species is a variable and its mean would be different from others, the best way to do it is pass a vector of means, and vector of size
    # This function does NOT work with single values of size, prob, munb, and pstr0. Need lists
    nSpecies = p.shape[1]
    nSamples = p.shape[0]
    if isinstance(size, float) or isinstance(size, int):
        #print('Need a list with one value per variable. Providing the argument nVariables  the value will be repeated nVariables times')
        size = [size] * nVariables # Repeat the same size nVariables time
    if isinstance(prob, float) or isinstance(prob, int):
        #print('Need a list with one value per variable. Providing the argument nVariables the value will be repeated nVariables times')
        prob = [prob] * nVariables
    if isinstance(munb, float) or isinstance(munb, int):
        #print('Need a list with one value per variable. Providing the argument nVariables the value will be repeated nVariables times')
        munb = [munb] * nVariables
    if isinstance(pstr0, float) or isinstance(pstr0, int):
        #print('Need a list with one value per variable. Providing the argument nVariables the value will be repeated nVariables times')
        pstr0 = [pstr0] * nVariables                  
    # 2. Repeated munb, size, prob and pstr0 by each value of the same variable
    if len(munb):
        prob = [s/(s + m) for s, m in zip(size, munb)]
    # Number of values
    LLL = max(len(p.flatten()), len(prob), len(pstr0), len(size), len(munb))
    p = np_rep_len(p.flatten(), LLL)
    if len(pstr0) != LLL:
        pstr0 = np_rep_len(pstr0, LLL)
    if len(prob) != LLL:
        prob = np_rep_len(prob, LLL)
    if len(size) != LLL:
        size = np_rep_len(size, LLL)
    if len(munb) != LLL:
        munb = np_rep_len(munb, LLL)
    # 3. Now everything have the proper length (same values -> same distribution for each variable). 
    # 3.1 Create empty list (should be list to mix 'NA' strings with float, in numpy.array is valid just one datatype)
    ans = list(np.repeat(float('nan'), LLL))
    prob0 = [p**s for p, s in zip(prob, size)]
    deflat_limit=[]
    for i in range(len(prob0)):
        if (1 - prob0[i]) == 0: #1- prob0[i] = 0 only when prob0[i] = 1
            deflat_limit.append(float('-inf'))
            #elif prob0[i] < 0: 
            #    deflat_limit[i] = float('inf')
        elif ((1 - prob0[i]) == 0 and prob0[i] == 0):
            deflat_limit.append(float('nan'))
        else:
            deflat_limit.append(-prob0[i] / (1 - prob0[i]))
    for i in range(len(ans)):
        if p[i] <= pstr0[i]:
            ans[i] = 0
    ind4 =  [(pstr0[i] < p[i]) and (deflat_limit[i] <= pstr0[i]) for i in range(len(p))]
    q = [(p[i] - pstr0[i]) / (1 - pstr0[i]) for i in range(len(p)) if ind4[i]]
    n = [size[i] for i in range(len(size)) if ind4[i]]
    pr = [prob[i] for i in range(len(prob)) if ind4[i]]
    j = 0
    for i in range(len(ind4)):
        if ind4[i]:
            ans[i] = nbinom.ppf(q = q[j], n = n[j], p = pr[j], loc = 0) # This function is not exactly equal to R function. R function return 0 and warnings in some cases, an python return nan
            j = j +1      
    for i in range(len(ans)):
        if pstr0[i] < deflat_limit[i]:
            ans[i] = float('nan')
        if 1 < pstr0[i]:
            ans[i] = float('nan')
        if p[i] < 0:
            ans[i] = float('nan')
        if 1 < p[i]:
            ans[i] = float('nan')
    return(np.array(ans).reshape(nSamples, nSpecies))

def np_rep_len(x, length_out):
    """ Python version of R rep_len (!= of rep) """ # From Internet
    # Flat a numpy matrix into list. Repeat matrix of values (list of list in numpy array) by column (taking 0 element from all the rows, then the 1st element, and so on)
    # Imagine 3 samples, 10 species [S1[sp0,sp1,sp2],S2[sp0,sp1,sp2],S3[sp0,sp1,sp2]]
    # Flat list should have 30 values. np.tile will do 30 * 30 values
    # To just flat the list, it would be necessary 30 * (0+1), (len(x) = 30) hace la division de integers (//)
    return(np.tile(x, length_out // len(x) + 1)[:length_out])

def lognormal(mu, sd, n):
    """Sampling n elements from a lognormal distribution with mu, sd (according to the normal distribution it is derived from)"""
    # NOTE:
    #"""log-normal distribution with specified mean, standard deviation, and array shape.
    #Note that the mean and standard deviation are not the values for the distribution itself,
    #but of the underlying normal distribution it is derived from."""
    # When sigma is low: all the abundances are similar. When it is high, there is a sequence that prevails
    #mu, sd = lognormal2normal(muLOG, sdLOG)
    abunlognormal = np.random.lognormal(mu, sd, n) 
    return(abunlognormal)

def relative_abundance(n, total):
    """ Relative abundance """
    if n == 0:
        return 0
    else:
        return (float(n)/total)

def shannonIndexCalc(data):
    """ Given a list of abundance values , returns the Shannon real Index """
    total = sum(data)
    return(-sum(relative_abundance(n, total) * np.log(float(n)/total) for n in data if n != 0.0))


