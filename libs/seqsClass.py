#!/usr/bin/env python3

from libs.libs import *
import re
import random

class Sequence():
    
    def __init__(self, header, seq, tax = None, abun = 0):
        self.header = header
        self.seq = seq
        self.tax = tax
        self.abun = abun
    
    def __repr__(self): 
        return('{}\n{}\n'.format(self.header,self.seq[0:1000]))    
    
    def copy(self):
        """ Copy the object to avoid modifying itself """
        SequenceCopy = Sequence(self.header, self.seq, self.tax, self.abun)
        return(SequenceCopy)
    
    def nt2dict(self):
        """ Convert 'AAG-A' => {'A':{1,2,5}, '.':set(), '-':set(4), 'C':set(), 'T':set(), 'G':set(3)} dict with a set of positions by nt """
        positions = {'A':set(), 'T':set(), 'G':set(), 'C':set(), '-': set(), '.': set()}
        for i, nt in enumerate(self.seq):
            positions[nt].add(i)
        return(positions)
    
    def trimregion(self, start, end):
        """ Trim a sequence according to start, end positions """
        start, end = int(start), int(end)
        if self.seq[start:end].replace('.','').replace('-',''): # seq has something in that region, not only '....' or '---'
            seq = self.seq[start:end]
            return(Sequence(self.header, seq, self.tax))
    
    def deGap(self): 
        """ Remove '.' and '-' from a sequence """
        seq = self.seq.replace('.','').replace('-','')
        return(Sequence(self.header, seq, self.tax))
        
    def generatemutantASVs(self, Nstrains = None, Nmean = 2, Nposmax = 45, start = 0, end = 50000, include_original = True): # Nmax: numero de cepas maximas por specie, Nposmax: max n? de posiciones que pueden ser mutadas #WORK
        """ For a given Sequence object, return a set of Sequence objects with one object per 'fake' Nstrains """
        if not Nstrains: # number of exact ASVs per sequence
            Nstrains = np.random.random_integers(1, Nmean*2)
        originalSeqName = self.header
        #  Create a set with different Seq objects, the Seq original object, and the strains if it is not the same
        Npos = random.randint(1, Nposmax)
        clusterSeqs = set()
        for i in range(0, Nstrains):
            if include_original:
                new_header = '{}.asv_{}'.format(originalSeqName, i+1)
            else:
                new_header = '{}-{}'.format(originalSeqName, i+1)
            newSequence = Sequence(new_header, mutate(string = self.seq, N = Npos, start = start, end = end), self.tax)
            clusterSeqs.add(newSequence) 
        if include_original:
            clusterSeqs.add(Sequence(originalSeqName, self.seq, self.tax))
        return(clusterSeqs)  # set with different Seq objects from a single Seq object


    @staticmethod
    def assign_taxonomy(header, silva_taxa, taxlevel = None): 
        """ Assign taxonomy from silva database """
       ## Including my fake ASVs names
        if re.findall('\.asv\_[0-9]*$', header): # Taxonomy of my mutant ASVs
            header = re.sub('\.asv\_[0-9]*$','', header)
        if re.findall('\-[0-9]*$', header): # Taxonomy of my mutant taxa
            header = re.sub('\-[0-9]*$','', header)
        taxonomy = silva_taxa[header]
        return(taxonomy)
