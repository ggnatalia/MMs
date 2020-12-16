#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import random
import scipy.stats as sp
import math
import os
from collections import defaultdict
import daiquiri
import logging

from libs.plots import *
from libs.maths import *
from libs.libs import *
from libs.seqsClass import Sequence

class Sample():
    
    
    def __init__(self, prefix, sampleName, Seqs, alignment, figsize):
        
        self.prefix = prefix
        self.sample_name = sampleName
        self.Seqs = Seqs
        self.pwd = os.getcwd()
        self.alignment = alignment # list [start, end, region]
        self.figsize = figsize
        
        
        
    def __repr__(self): 
        return(self.sample_name)
    
    def run_sample(self):
        """ Prepare each sample individually """
        write_logfile('info', 'PROCESSING SAMPLE {}'.format(self.sample_name), 'Creating files to run InSilicoSeq')
        self.generate_inSilicoSeq_input_files()
        #self.run_inSilicoSeq(prefix, reads, errormodel = errormodel, cpus = cpus ) #########33 A LO MEJOR LO SACO DE AQUI Y LO HAGO EN GENERAL
        write_logfile('info', 'PROCESSING SAMPLE {}'.format(self.sample_name), 'Plotting abundances')
        self.plot_abundances()
        #write_logfile('info', 'PROCESSING SAMPLE {}'.format(self.sample_name), 'Plotting entropies')
        #self.sample_entropy()
        write_logfile('info', 'PROCESSING SAMPLE {}'.format(self.sample_name), 'Writing align file for the required region')       
        self.write_align_by_region()
        return('{}.{}.sequences16S.fasta'.format(self.prefix, self.sample_name), '{}.{}.abundances'.format(self.prefix, self.sample_name), self.sample_name) 
    def write_align_by_region(self):
        """Write the alignment of the chosen region- for checking purposes"""
        if self.alignment:
            start, end, region  = int(self.alignment[0]), int(self.alignment[1]), str(self.alignment[2])
            sampleAlignRegion = '{}/checkDB/{}.{}.{}.align'.format(self.pwd, self.prefix, self.sample_name, region)
            with open(sampleAlignRegion, 'w') as outputAlign:
                for s in self.Seqs:
                    if s.seq[start:end].replace('.','').replace('-',''): # seq has something in that region
                        outputAlign.write('>{}\n{}\n'.format(s.header, s.seq[start:end]))
            #Always do this:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! QUITAR EN UN FUTURO, por ahora dejar
            if region != '16S':
                sampleAlignRegion = '{}/checkDB/{}.{}.{}.align'.format(self.pwd, self.prefix, self.sample_name, region)
                with open(sampleAlignRegion, 'w') as outputAlign2:
                    for s in self.Seqs:
                        if s.seq[0:5000].replace('.','').replace('-',''): # seq has something in that region
                            outputAlign2.write('>{}\n{}\n'.format(s.header, s.seq[start:end]))
        else:
            pass
            
    def generate_inSilicoSeq_input_files(self):
        """Create what should be the true input files to inSilicoSeqs: fasta, abundances."""
        sampleFasta = '{}/samples/{}.{}.sequences16S.fasta'.format(self.pwd, self.prefix, self.sample_name)
        sampleAbun = '{}/samples/{}.{}.abundances'.format(self.pwd, self.prefix, self.sample_name)
        #abundSeqs = sorted(tuple(self.abundSeqs.items()), key = lambda x: x[1], reverse=True)
        ## Write output files
        with open(sampleFasta, 'w') as seqs, open(sampleAbun, 'w') as abundances: # ADD NNNNNNNNNNNNNNNNNNNNN to simulate part of the reads that are no 16S
            for s in self.Seqs:
                seqs.write('>{}\n{}{}{}\n'.format(s.header, 'N'*500, s.deGap().seq,'N'*500))
                abundances.write('{}\t{}\n'.format(s.header, s.abun/sum([s.abun for s in self.Seqs]))) #remove round(,2)
 
    def plot_abundances(self):
        """ Plot abundance distribution """
        prefix = '{}.{}'.format(self.prefix, self.sample_name)
        outputDir = '{}/samples'.format(self.pwd)
        total = sum([x.abun for x in self.Seqs])
        plotData = sorted({x.header:(x.abun/total)*100 for x in self.Seqs}.items(), key = lambda x : x[1], reverse = True)
        write_logfile('info', 'PROCESSING SAMPLE', 'This is the real shannon index diversity in this sample: {}'.format(str(shannonIndexCalc([v[1] for v in plotData])))) 
        plotDataDF = pd.DataFrame(plotData, columns = ['seqName', 'Percentage'])
        plotDataDF = plotDataDF.set_index('seqName')
        barplot(plotDataDF, outputDir = outputDir, title = '{}.abundances'.format(prefix), T=True, ylab = 'Percentage reads', xlab = 'Sequences', textSize=8, figsize = self.figsize)
    # DE MOMENTO INUTIL 
    def sample_entropy(self):
        """ Calculate entropy of one sample. Take into account coverage """
        prefix = '{}.{}'.format(self.prefix, self.sample_name)
        outputDir = '{}/samples'.format(self.pwd)
        #seqsArray = self.align_to_array()
        #nrow, ncol = np.shape(seqsArray)
        length = len(set([len(s.seq)for s in self.Seqs]))
        # All the squences in silva have length: 50000
        if length == 1: 
            #I take into account if a seq repeats n times.
            #Calculate entropyes ignoring dots
            #entropies = {pos:sp.entropy(list(map(float, dict(sp.itemfreq(np.repeat(seqsArray[:,pos], abundances))).values()))) for pos in range(0,ncol)}
            #entropies = {pos:sp.entropy(list(map(float, dict(sp.itemfreq(np.delete(self.expand_column(seqsArray[:,pos], abundances), np.where(self.expand_column(seqsArray[:,pos], abundances) == '.')))).values()))) for pos in range(ncol)}
            # Expanding all the alignment previously
            seqsArray = align_to_array(self.Seqs, expand = True)
            #entropies = {pos:sp.entropy(list(map(float, dict(sp.itemfreq(np.delete(seqsArray[:,pos], np.where(seqsArray[:,pos] == '.')))).values()))) for pos in range(ncol)}
            plot_entropy(seqsArray, title = prefix, outputDir = outputDir)
        else:
            write_logfile('error', 'PROCESSING sample_entropy', 'Your sequences have different lengths, align them before please')
            
            
    @classmethod
    def init_from_df(cls, prefix, sample, Seqs, alignment, figsize):
        """ Create a Sample object from pandas df 1D and alignSeqs"""
        #sample = row of the df, name = 'S_0'
        sampleName = sample.name
        Seqscopy = [s.copy() for s in Seqs]
        seqsNull = set()
        for s in Seqscopy:
            s.abun = sample.loc[s.header] # Add abun to each object of the Seqs object. #remove round(,2)
            if s.abun == 0:
                seqsNull.add(s.header)
        SeqsSample = {s for s in Seqscopy if s.header not in seqsNull}
        return(Sample(prefix, sampleName, SeqsSample, alignment, figsize))



    @staticmethod
    def run_inSilicoSeq(prefix, seqFile, abunFile, sampleName, reads = 10000, errormodel = 'perfect', cpus = 12, InSilicoparams = (150,200)):  
        """Run inSilicoSeq"""
        # Run inSilicoSeq in samples directory from command line to generate the fasta samples
        command = ['iss', 'generate', '-g', seqFile, '--abundance_file', abunFile, '-o', '{}.{}.{}-{}-{}r-{}i.InSilicoSeq'.format(prefix, sampleName, reads, errormodel, InSilicoparams[0], InSilicoparams[1]), '--n_reads', reads, '--cpus', cpus, '--InSilicoparams', str(InSilicoparams[0]), str(InSilicoparams[1])] # Reads are the total count both pairs, to have the reads required by the user in each file
        if errormodel == 'MiSeq' or errormodel == 'HiSeq' or errormodel == 'NovaSeq':
            command.append('--model')
            command.append(errormodel)
        else:
            command.append('--mode')
            command.append(errormodel) 
        write_logfile('info', 'SAMPLE READS GENERATION', 'Read generation: {}'.format(' '.join(map(str, command))))
        runCommand(command, open('{}.{}.inSilicoSeq.logfile'.format(prefix, sampleName), 'w'))
    
    
             
    
    
    
    
