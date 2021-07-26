#!/usr/bin/env python
import pandas as pd
#import matplotlib.pyplot as plt
import numpy as np
import random
import scipy.stats as sp
import math
import os
import sys
from collections import defaultdict
import daiquiri
import logging

from libs.plots import *
from libs.maths import *
from libs.libs import *
from libs.seqsClass import Sequence


class Sample():
    
    #print('{}/../extlibs/'.format('/'.join(os.path.abspath(__file__).split('/')[:-1])))
    #sys.path.append('{}/../extlibs/'.format('/'.join(os.path.abspath(__file__).split('/')[:-1])))
    def __init__(self, prefix, sampleName, Seqs, alignment, reads, Sim):
        
        self.prefix = prefix
        self.sample_name = sampleName
        self.Seqs = Seqs
        self.pwd = os.getcwd()
        self.alignment = alignment # list [start, end, region]
        self.reads = reads
        self.Sim = Sim
        
    def __repr__(self): 
        return(self.sample_name)
    
    def run_sample(self):
        """ Prepare each sample individually """
        write_logfile('info', 'PROCESSING SAMPLE {}'.format(self.sample_name), 'Creating files to run the sequencing simulator')
        self.generate_Simulator_input_files()
        #self.run_inSilicoSeq(prefix, reads, errormodel = errormodel, cpus = cpus ) #########33 A LO MEJOR LO SACO DE AQUI Y LO HAGO EN GENERAL
        write_logfile('info', 'PROCESSING SAMPLE {}'.format(self.sample_name), 'Plotting abundances')
        self.plot_abundances()
        #write_logfile('info', 'PROCESSING SAMPLE {}'.format(self.sample_name), 'Plotting entropies')
        write_logfile('info', 'PROCESSING SAMPLE {}'.format(self.sample_name), 'Writing align file for the required region')       
        self.write_align_by_region()
        if self.Sim == 'InSilicoSeqs':
            return('{}.{}.sequences16S.fasta'.format(self.prefix, self.sample_name), '{}.{}.abundances'.format(self.prefix, self.sample_name), self.sample_name)
        elif self.Sim == 'NanoSim':
            return('{}.{}.NSgenomes'.format(self.prefix, self.sample_name), '{}.{}.NSabundances'.format(self.prefix, self.sample_name), '{}.{}.NSdl'.format(self.prefix, self.sample_name), self.sample_name)
    
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
            
    def generate_Simulator_input_files(self):
        """Create what should be the true input files to inSilicoSeqs: fasta, abundances."""
        sampleFasta = '{}/samples/{}.{}.sequences16S.fasta'.format(self.pwd, self.prefix, self.sample_name)
        if self.Sim == 'InSilicoSeqs':
            sampleAbun = '{}/samples/{}.{}.abundances'.format(self.pwd, self.prefix, self.sample_name)
            #abundSeqs = sorted(tuple(self.abundSeqs.items()), key = lambda x: x[1], reverse=True)
        ## Write output files
            if not os.path.isfile(sampleFasta) or not os.path.isfile(sampleAbun):
                with open(sampleFasta, 'w') as seqs, open(sampleAbun, 'w') as abundances: # ADD NNNNNNNNNNNNNNNNNNNNN to simulate part of the reads that are no 16S
                    for s in self.Seqs:
                        seqs.write('>{}\n{}{}{}\n'.format(s.header, 'N'*500, s.deGap().seq,'N'*500))
                        abundances.write('{}\t{}\n'.format(s.header, s.abun/sum([s.abun for s in self.Seqs]))) #remove round(,2)
        elif self.Sim == 'NanoSim':
            sampleGL = '{}/samples/{}.{}.NSgenomes'.format(self.pwd, self.prefix, self.sample_name)
            sampleAbun = '{}/samples/{}.{}.NSabundances'.format(self.pwd, self.prefix, self.sample_name)
            sampleDL = '{}/samples/{}.{}.NSdl'.format(self.pwd, self.prefix, self.sample_name)
            if not os.path.isfile(sampleFasta) or not os.path.isfile(sampleAbun) or not os.path.isfile(sampleGL) or not os.path.isfile(sampleDL):
                with open(sampleFasta, 'w') as seqs, open(sampleAbun, 'w') as abundances, open(sampleDL, 'w') as dlNS, open(sampleGL, 'w') as genomes: # ADD NNNNNNNNNNNNNNNNNNNNN to simulate part of the reads that are no 16S
                    abundances.write('Size\t{}\n'.format(self.reads))
                    for s in self.Seqs:
                        seqs.write('>{}\n{}{}{}\n'.format(s.header, 'N'*500, s.deGap().seq,'N'*500))
                        abundances.write('{}\t{}\n'.format(s.header, s.abun/sum([s.abun for s in self.Seqs]*100))) #remove round(,2)
                        dlNS.write('{}\t{}\tlinear\n'.format(s.header, s.header))      
                        genomes.write('{}\t{}\n'.format(s.header, sampleFasta))
    
    def plot_abundances(self):
        """ Plot abundance distribution """
        prefix = '{}.{}'.format(self.prefix, self.sample_name)
        outputDir = '{}/samples'.format(self.pwd)
        total = sum([x.abun for x in self.Seqs])
        plotData = sorted({x.header:(x.abun/total)*100 for x in self.Seqs}.items(), key = lambda x : x[1], reverse = True)
        write_logfile('info', 'PROCESSING SAMPLE', 'This is the real shannon index diversity in this sample: {}'.format(str(shannonIndexCalc([v[1] for v in plotData])))) 
        plotDataDF = pd.DataFrame(plotData, columns = ['seqName', 'Percentage'])
        plotDataDF = plotDataDF.set_index('seqName')
        barplot(plotDataDF.T, outputDir = outputDir, title = '{}.abundances'.format(prefix), ylab = 'Percentage reads', xlab = 'Sequences')
            
            
    @classmethod
    def init_from_df(cls, prefix, sample, Seqs, alignment, reads, Sim):
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
        return(Sample(prefix, sampleName, SeqsSample, alignment, reads, Sim))



    @staticmethod
    def run_inSilicoSeq(prefix, seqFile, abunFile, sampleName, reads = 10000, cpus = 12, ISSerrormodel = 'perfect', ISSparams = (150,200)):  
        """Run inSilicoSeq"""
        # Run inSilicoSeq in samples directory from command line to generate the fasta samples
        command = ['{}/../bin/iss'.format('/'.join(os.path.abspath(__file__).split('/')[:-1])), 'generate', '-g', seqFile, '--abundance_file', abunFile, '-o', '{}.{}.{}-{}-{}r-{}i.InSilicoSeq'.format(prefix, sampleName, reads, ISSerrormodel, ISSparams[0], ISSparams[1]), '--n_reads', reads, '--cpus', cpus, '--inSilicoparams', str(ISSparams[0]), str(ISSparams[1])] # Reads are the total count both pairs, to have the reads required by the user in each file
        if ISSerrormodel == 'MiSeq' or ISSerrormodel == 'HiSeq' or ISSerrormodel == 'NovaSeq':
            command.append('--model')
            command.append(ISSerrormodel)
        else:
            command.append('--mode')
            command.append(ISSerrormodel) 
        write_logfile('info', 'SAMPLE READS GENERATION', 'Read generation: {}'.format(' '.join(map(str, command))))
        runCommand(command, open('{}.{}.inSilicoSeq.logfile'.format(prefix, sampleName), 'w'))
    
    
    @staticmethod
    def run_NanoSim(prefix, seqFile, abunFile, dlFile, sampleName, NSerrormodel = 'perfect', reads = 10000, cpus = 12, NSparams = (1500,50)):
        params = ['-max', NSparams[0], '-min', NSparams[1], '-c', '{}/../../NanoSim/pre-trained_models/metagenome_ERR3152364_Even/training'.format('/'.join(os.path.abspath(__file__).split('/')[:-1]))]
        command = ['python', '{}/../../NanoSim/src/simulator.py'.format('/'.join(os.path.abspath(__file__).split('/')[:-1])), 'metagenome', '-b', 'guppy', '--fastq', '-gl',  seqFile, '-a', abunFile, '-dl', dlFile, '-t', cpus]
        command = command + params
        if NSerrormodel == 'perfect':
            command.append('--perfect')
        else:
            NSerrormodel = 'metagenome'
        command.append('-o')
        command.append('{}.{}.{}-{}-{}max-{}min.NanoSim'.format(prefix, sampleName, reads, NSerrormodel, NSparams[0], NSparams[1]))
        write_logfile('info', 'SAMPLE READS GENERATION', 'Read generation: {}'.format(' '.join(map(str, command))))
        runCommand(command, open('{}.{}.NanoSim.logfile'.format(prefix, sampleName), 'w'))

