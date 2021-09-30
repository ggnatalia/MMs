#!/usr/bin/env python

#import sys
#sys.path.append('/media/natalia/Linux/natalia/opt/makemocks/')

from libs.libs import *
from libs.maths import *
import re
import random
import numpy as np



class Sequence():
    
    #def nt2dict(self):
    #    positions = {'A':set(), 'T':set(), 'G':set(), 'C':set(), '-': set(), '.': set()}
    #    for i, nt in enumerate(self.seq):
    #        positions[nt].add(i)
    #    return(positions)
    
    def nt2array(self):
        dim = len(self.seq)
        equiv = {'A':1, 'T':2, 'G':3, 'C':4, '-':5, '.':0}
        seq2array = np.zeros((dim), dtype = np.uint8)
        for i, nt in enumerate(self.seq):
            seq2array[i] = equiv[nt]
        return(seq2array) 
    
    
    
    def __init__(self, header, seq, tax = None, abun = 0):
        self.header = header
        self.seq = seq
        self.tax = tax
        self.abun = abun
    
    def __repr__(self): 
        return('{}\n{}\n'.format(self.header,self.seq[0:1000]))    
    
    def copy(self):
        SequenceCopy = Sequence(self.header, self.seq, self.tax, self.abun)
        return(SequenceCopy)
    
    def trimregion(self, start, end):
        """ Trim a sequence according to start, end positions """
        start, end = int(start), int(end)
        if self.seq[start:end].replace('.','').replace('-',''): # seq has something in that region, not only '............' or '----------------'
            seq = self.seq[start:end]
            return(Sequence(self.header, seq, self.tax))
        return(Sequence(None, None, None))
    
    def deGap(self): 
        """ Remove '.' and '-' from a sequence """
        seq = self.seq.replace('.','').replace('-','')
        return(Sequence(self.header, seq, self.tax))
        
    def generatemutantASVs(self, Nstrains = None, Nmean = 2, Nposmax = 45, start = 0, end = 50000, by_region = None, by_pos = None, include_original = True): # Nmax: numero de cepas maximas por specie, Nposmax: max n? de posiciones que pueden ser mutadas #WORK
        """ For a given Sequence object, return a set of Sequence objects with one object per 'fake' Nstrains """
        if not Nstrains: #Numero de cepas exactas
            Nstrains = np.random.random_integers(0, (Nmean-1)*2) #  to be exact with ASVmean: (Nmean-1)*2??
        originalSeqName = self.header
        #  Create a set with different Seq objects, the Seq original object, and the strains if it is not the same
        Npos = random.randint(1, Nposmax)
        #print('Nposmax ' + str(Nposmax) + 'Npos ' + str(Npos))
        clusterSeqs = set()
        if include_original:
            clusterSeqs.add(Sequence(originalSeqName, self.trimregion(start, end).seq, self.tax)) ##### NEW LINE TO ADD THE REAL STRAIN TO THE MOCK
            #clusterSeqs.add(Sequence(originalSeqName, self.seq, self.tax)) ##### NEW LINE TO ADD THE REAL STRAIN TO THE MOCK
            for i in range(0, Nstrains):   
                newSequence = Sequence('{}.asv_{}'.format(originalSeqName, i+1),  mutate(string = self.seq, N = Npos, start = start, end = end, regions = by_region, positions = by_pos, header = '{}.asv_{}'.format(originalSeqName, i+1)), self.tax)
                clusterSeqs.add(newSequence) # Add new Seq objects
        else:
            for i in range(0, Nstrains):   
                newSequence = Sequence('{}-{}'.format(originalSeqName, i+1),  mutate(string = self.seq, N = Npos, start = start, end = end, regions = by_region, positions = by_pos, header = '{}-{}'.format(originalSeqName, i+1)), self.tax)
                clusterSeqs.add(newSequence) # Add new Seq objects
        return(clusterSeqs)         # set with different Seq objects from a single Seq object

    @staticmethod
    def makedistances(s):
        """ Calculate distance between two Sequence object in a tuple """
        return(s[0].header, s[1].header, calculate_distance(s[0].seq, s[1].seq))

    @staticmethod
    def assign_taxonomy(header, silva_taxa, taxlevel = None): 
        """ Assign taxonomy from silva database """
        #silva_taxa = loadTaxa(refTax, taxlevel = taxlevel) # dict with the taxonomy of all the sequences from silva
       ####################################### BASED ON SILVA NAMES
        if re.findall('\.asv\_[0-9]*$', header): #Taxonomy of my mutant ASVs
            header = re.sub('\.asv\_[0-9]*$','', header)
        if re.findall('\-[0-9]*$', header): #Taxonomy of my mutant taxa ASVs
            header = re.sub('\-[0-9]*$','', header)
        taxonomy = silva_taxa[header]
        return(taxonomy)

