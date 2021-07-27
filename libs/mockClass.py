#!/usr/bin/env python
from libs.libs import *#Library with general functions
from libs.maths import *
from libs.plots import *
from sklearn.datasets import make_sparse_spd_matrix
import numpy as np
import math
import pandas as pd
import os
import random
import shutil
from libs.sampleClass import Sample  # verify
from libs.seqsClass import Sequence
from collections import defaultdict, Counter
import logging
import daiquiri
from multiprocessing import Pool

class Mock():
    #CLASS VARIABLE
    multiprocessing_globals_samples = []
    SimInput = []
    
    @classmethod
    def clear_multiprocessing_globals(cls): # Clean class variable instead of overwrite it
        cls.multiprocessing_globals_samples = []
        cls.SimInput = []
        
    def __init__(self, prefix, Seqs, df, samples, alignment = [0, 50000, '16S'], cpus = 12, reads = 10000, Sim = 'InSilicoSeqs', ISSerrormodel = 'perfect', ISSparams = (150,200), NSerrormodel = 'perfect', NSparams = (500, 50)):
        self.prefix = prefix
        self.samples = samples
        self.df = df # row: samples, col: species
        self.gSeqs = Seqs
        self.reads = reads
        self.cpus = cpus
        self.alignment = alignment
        self.Sim = Sim
        self.ISSerrormodel = ISSerrormodel
        self.ISSparams = ISSparams
        self.NSerrormodel = NSerrormodel
        self.NSparams = NSparams
    
    def __repr__(self): 
        return(self.df)
    
    @classmethod    
    def init_and_run_InSilico(cls, prefix, ISSsequence_files, ISSabundance_files, alignment = [0, 50000, '16S'], cpus = 12, reads = 20000, Sim = 'InSilicoSeqs', ISSerrormodel = 'perfect', ISSparams = (150, 200)):
        mock = Mock(prefix, set(), pd.DataFrame(index=[], columns=[]), [], alignment = alignment, cpus = cpus, reads = reads, Sim = Sim, ISSerrormodel = ISSerrormodel,  ISSparams = ISSparams)
        for i, (fasta, abun) in enumerate(zip(ISSsequence_files, ISSabundance_files)):
            sampleName = 'S_' + str(i)
            mock.SimInput.append((fasta, abun, sampleName))
        return(mock.reads_generation())
    
    @classmethod    
    def init_and_run_NanoSim(cls, prefix, NSsequence_files, NSabundance_files, NSdl_files, alignment = [0, 50000, '16S'], cpus = 12, reads = 10000, Sim = 'NanoSim', NSerrormodel = 'perfect', NSparams = (500, 50)):
        mock = Mock(prefix, set(), pd.DataFrame(index=[], columns=[]), [], alignment = alignment, cpus = cpus, reads = reads, Sim = Sim, NSerrormodel = NSerrormodel, NSparams = NSparams)
        for i, (fasta, abun, dl) in enumerate(zip(NSsequence_files, NSabundance_files, NSdl_files)):
            sampleName = 'S_' + str(i)
            mock.SimInput.append((fasta, abun, dl, sampleName))
        return(mock.reads_generation())
    
    @classmethod
    def init_from_previousmock(cls, prefix, sequence_file, abun_table, alignment = [0, 50000, '16S'], cpus = 12, reads = 20000, Sim = 'InSilicoSeqs', ISSerrormodel = 'perfect',  ISSparams = (150,200), NSerrormodel = 'perfect',  NSparams = (500,50)):
        write_logfile('info', 'LOADING A PREVIOUS MOCK', 'Loading sequence and abundance files')
        #Check if the files exists:
        gSeqs = set()
        df = read_table(abun_table)
        #print(df)
        with open(sequence_file, 'r') as fasta:
            print(sequence_file)
            while True:
                header = fasta.readline().lstrip('>').rstrip('\n')
                if not header:
                    break
                seq = fasta.readline().rstrip('\n')
                if not 'N' in seq:
                    gSeqs.add(Sequence(header, seq, None, 0))
        mock = Mock(prefix, gSeqs, df, [], alignment = alignment, cpus = cpus, reads = reads, Sim = Sim, ISSerrormodel = ISSerrormodel, ISSparams = ISSparams, NSerrormodel = NSerrormodel, NSparams = NSparams)
        mock.multiprocessing_globals_samples = [Sample.init_from_df(prefix = prefix, sample = df.loc[sample], Seqs = gSeqs, alignment = alignment, reads = reads, Sim = Sim) for sample in list(df.index.values)] # List of sample objects
        if cpus == 1:
            write_logfile('info', 'Creating input files for squencing simulator', 'Start 1 cpu')
            mock.SimInput = list(map(Sample.generate_Simulator_input_files(), mock.multiprocessing_globals_samples)) #[('seq', 'abun', sample_name),()]
        else:
            write_logfile('info', 'Creating input files for squencing simulator', 'Start multiprocessing')
            with Pool(mock.cpus) as pool:
                mock.SimInput = list(pool.map(Sample.generate_Simulator_input_files, mock.multiprocessing_globals_samples))
        write_logfile('info', 'CREATE MOCK', 'Running Simulator')
        #print(mock.SimInput)
        mock.reads_generation()
        mock.clear_multiprocessing_globals()
        return(mock)




    @classmethod
    def init_from_enviro(cls, Enviro, prefix, Nsamples, shannonIndex, alignment = [0, 50000, '16S'], cpus = 12, alpha = 0.9,  smallest_coef = 0.1, largest_coef = 0.9, reads = 20000, pstr0 = 0.2, size = 1, Sim = 'InSilicoSeqs', ISSerrormodel = 'perfect',  ISSparams = (150,200), NSerrormodel = 'perfect',  NSparams = (500,50)):
        """ Creating a mock from scratch. A set of Samples objects"""
        write_logfile('info', 'CREATE MOCK', 'Estimating abundances')
        if Nsamples == 1:
            abunTable = np.array(list(map(round, cls.make_global_samples_distribution(shannonIndex = shannonIndex, nASVs = len(Enviro.Seqs), reads = reads).T )))#lognormal distribution
            df =  pd.DataFrame(abunTable.reshape(Nsamples, len(Enviro.Seqs)), index = ['S_' + str(i) for i in range(Nsamples)] , columns = [s.header for s in Enviro.Seqs])
        else:
            CorrMatrix = cls.make_correlation_matrix(asvs = [s.header for s in Enviro.Seqs], prefix = prefix, norm_diag = 1, alpha = alpha, smallest_coef = smallest_coef, largest_coef = largest_coef)
            abunTable = cls.ZINBD(nASVs = len(Enviro.Seqs), nSamples = Nsamples, CorrMatrix = CorrMatrix, shannonIndex = shannonIndex, reads = reads, pstr0 = pstr0, size = size)
            df =  pd.DataFrame(abunTable.reshape(Nsamples, len(Enviro.Seqs)), index = ['S_' + str(i) for i in range(Nsamples)] , columns = [s.header for s in Enviro.Seqs]) 
            write_logfile('info', 'CREATE MOCK', 'Plotting abundance correlation')
            plot_heatmap(df.corr(method ='spearman'), outputDir = os.getcwd(), title = '{}.correlationMatrix_fromAbunTable'.format(prefix), zmin = -1, zmax = 1, legendtitle = 'Correlation', symmetric = True)
        # To each sequence, add they global abundance in all the samples
        for s in Enviro.Seqs:
            s.abun = df[s.header].sum()
        mock = Mock(prefix, Enviro.Seqs, df, [], alignment = alignment, cpus = cpus, reads = reads, Sim = Sim, ISSerrormodel = ISSerrormodel, ISSparams = ISSparams, NSerrormodel = NSerrormodel, NSparams = NSparams)
        return(mock.run_mock()) #samples is empty, in run samples, it will be filled
        
    def run_mock(self):
        """ Work sample by sample: create input files for InSilicoSeqs, plots and run InSilicoSeqs """
        self.samples = [Sample.init_from_df(prefix = self.prefix, sample = self.df.loc[sample], Seqs = self.gSeqs, alignment = self.alignment, reads = self.reads, Sim = self.Sim) for sample in list(self.df.index.values)] # List of sample objects
        write_logfile('info', 'CREATE MOCK', 'Writing output abundance tables')
        write_table(self.df, title = 'checkDB/{}.raw.abundances_original'.format(self.prefix), rows = list(self.df.index) , columns = list(self.df.columns), dataframe = None)
        abunTablePercent = make_percent(self.df, outputDir = 'checkDB/', write = True, fileName = '{}.abundances_original'.format(self.prefix), T = False)
        write_logfile('info', 'CREATE MOCK', 'This is the real shannon index diversity in all the mock: {}'.format(str(shannonIndexCalc(self.df.sum()))))
        barplot(self.df, outputDir = os.getcwd(), title = '{}.sampleDistribution'.format(self.prefix),  ylab = 'Abundances', xlab = 'Species', subtitle = True) # Plot asvs distribution by sample: lognormal
        if len(self.samples) > 1:
            barplot(self.df.T, outputDir = os.getcwd(), title = '{}.asvsDistribution'.format(self.prefix), ylab = 'Abundances', xlab = 'Samples', subtitle = False)  # Plot distribution of one asvs in the different samples: ZINBD
        write_logfile('info', 'CREATE MOCK', 'Preparing samples')
        self.multiprocessing_globals_samples = self.samples
        if self.cpus == 1:
            write_logfile('info', 'CONVERT SEQS', 'Start 1 cpu')
            self.SimInput = list(map(Sample.run_sample, self.multiprocessing_globals_samples)) #[('seq', 'abun', sample_name),()]
        else:
            write_logfile('info', 'CONVERT SEQS', 'Start multiprocessing')
            with Pool(self.cpus) as pool:
                self.SimInput = list(pool.map(Sample.run_sample, self.multiprocessing_globals_samples))
        write_logfile('info', 'CREATE MOCK', 'Running Simulator')
        self.reads_generation()
        self.clear_multiprocessing_globals()
        return(self)
               
    def reads_generation(self): 
        """ Call Sample.run_inSilicoSeq/Sample.run_NanoSim for each sample or sequence/abundance files """
        if self.Sim == 'InSilicoSeqs':
            for s in self.SimInput:
                Sample.run_inSilicoSeq(prefix = 'samples/{}'.format(self.prefix), seqFile = 'samples/{}'.format(s[0]), abunFile = 'samples/{}'.format(s[1]), sampleName = s[2], ISSerrormodel = self.ISSerrormodel, reads = self.reads, cpus = self.cpus, ISSparams = self.ISSparams)
        elif self.Sim == 'NanoSim':
            for s in self.SimInput:
                Sample.run_NanoSim(prefix = 'samples/{}'.format(self.prefix), seqFile = 'samples/{}'.format(s[0]), abunFile = 'samples/{}'.format(s[1]), dlFile = 'samples/{}'.format(s[2]), sampleName = s[3], NSerrormodel = self.NSerrormodel, reads = self.reads, cpus = self.cpus, NSparams = self.NSparams)
        write_logfile('info', 'MOCK REPEAT', 'Writing output files')
        self.merge_samples(samplesDir = 'samples')

    def merge_samples(self, samplesDir = 'samples'): ## TESTED!
        """ Merge samples from InSilicoSeqs and rename them (in this case it would not be necessary, but for coherence with the previous procedure) """
        
        
        if self.Sim == 'InSilicoSeqs':
            outputFastq = '{}/{}.allsamples.{}-{}-{}r-{}i.fastq'.format(samplesDir, self.prefix, self.reads, self.ISSerrormodel, self.ISSparams[0], self.ISSparams[1])
            outputFasta = '{}/{}.allsamples.{}-{}-{}r-{}i.fasta'.format(samplesDir, self.prefix, self.reads, self.ISSerrormodel, self.ISSparams[0], self.ISSparams[1])
            outputGroups = '{}/{}.allsamples.{}-{}-{}r-{}i.groups'.format(samplesDir, self.prefix, self.reads, self.ISSerrormodel, self.ISSparams[0], self.ISSparams[1])
            outputEquivalentNames = '{}/{}.allsamples.{}-{}-{}r-{}i.equivalence'.format(samplesDir, self.prefix, self.reads, self.ISSerrormodel, self.ISSparams[0], self.ISSparams[1])
        elif self.Sim == 'NanoSim':
            outputFastq = '{}/{}.allsamples.{}-{}-{}max-{}min.fastq'.format(samplesDir, self.prefix, self.reads, self.NSerrormodel, self.NSparams[0], self.NSparams[1])
            outputFasta = '{}/{}.allsamples.{}-{}-{}max-{}min.fasta'.format(samplesDir, self.prefix, self.reads, self.NSerrormodel, self.NSparams[0], self.NSparams[1])
            outputGroups = '{}/{}.allsamples.{}-{}-{}max-{}min.groups'.format(samplesDir, self.prefix, self.reads, self.NSerrormodel, self.NSparams[0], self.NSparams[1])
            outputEquivalentNames = '{}/{}.allsamples.{}-{}-{}max-{}min.equivalence'.format(samplesDir, self.prefix, self.reads, self.NSerrormodel, self.NSparams[0], self.NSparams[1])
        
        with open(outputFastq,'w') as outfastq, open(outputFasta, 'w') as outfasta, open(outputGroups,'w') as groups, open(outputEquivalentNames, 'w') as equivalence:
            parsedSeqs = []
            equivalences = {}
            for f in sorted(os.listdir(samplesDir)):
                files_header = ''
                if self.Sim == 'InSilicoSeqs':
                    files_header = '.{}-{}-{}r-{}i.InSilicoSeq'.format(self.reads, self.ISSerrormodel, self.ISSparams[0], self.ISSparams[1])
                elif self.Sim == 'NanoSim':
                    files_header = '.{}-{}-{}max-{}min.NanoSim_sample0_aligned'.format(self.reads, self.NSerrormodel, self.NSparams[0], self.NSparams[1])
                if files_header in f and f.endswith('.fastq'):
                    if '.'.join(f.split('.')[0:2]) == self.prefix:
                        sample = f.split('.')[2]
                        if self.Sim == 'InSilicoSeqs':
                            pair = f.split('.')[4].split('_')[1]
                        elif self.Sim == 'NanoSim':
                            pair = 'no'
                        print(f)
                        with open('{}/{}'.format(samplesDir,f),'r') as infastq:
                            name = 1
                            while True:
                                h = infastq.readline().rstrip('\n').lstrip('@')
                                if not h:
                                    break
                                s = infastq.readline().rstrip('\n')
                                c = infastq.readline().rstrip('\n')
                                qual = infastq.readline().rstrip('\n')
                                new_seqName = 'id__{}__seq_{}_pair_{}'.format(sample, name, pair) 
                                name = name + 1
                                if s.replace('N', ''):
                                    ambigpos = [pos for pos in range(len(s)) if s[pos] == 'N']
                                    s = s.replace('N', '') #remove all N for the sequences that exist
                                    qual = list(qual)
                                    for p in ambigpos:
                                        qual[p] = ''
                                    parsedSeqs.append('@{}\n{}\n{}\n{}\n'.format(new_seqName, s, c, ''.join(qual))) #list of complete seqs, seqName included,
                                    equivalences[new_seqName] = h
            permutedSeqs = np.random.permutation(list(parsedSeqs)) #list of complete seqs, seqName included, but shuffle
            for complete_seq in permutedSeqs:
                outfastq.write(complete_seq)
                seqName, seq, comment, qual = complete_seq.rstrip('\n').split('\n')
                outfasta.write('{}\n{}\n'.format(seqName.replace('@','>'), seq))
                seqName = seqName.replace('@','')
                group = seqName.split('__')
                groups.write('{}\t{}\n'.format(seqName, group[1]))
                equivalence.write('{}\t{}\t{}\n'.format(seqName, equivalences[seqName], '_'.join([equivalences[seqName].split('_')[0].lstrip('@'), equivalences[seqName].split('_')[1]])))
            print('Output files:\n{}\n{}\n{}\n{}\n'.format(outputFastq, outfasta, outputGroups, outputEquivalentNames))


    @staticmethod
    def make_correlation_matrix(asvs, prefix, norm_diag = 1, alpha = 0.9, smallest_coef = 0.1, largest_coef = 0.9):
        """ Create a correlation matrix: symmetric, definite positive (diagonal >0) sparse (many 0)"""
        # alpha: The probability that a coefficient is zero (see notes). Larger values enforce more sparsity.
        # norm_diag: Whether to normalize the output matrix to make the leading diagonal elements all 1
        # smallest_coef: The value of the smallest coefficient
        # largest_coef: The value of the largest coefficient
        # rows = ['Sp' + str(i) for i in range(nSpecies)]
        # columns = ['Sp' + str(i) for i in range(nSpecies)]
        corrMatrix = make_sparse_spd_matrix(len(asvs), alpha = alpha, norm_diag = norm_diag, smallest_coef = smallest_coef, largest_coef = largest_coef)
        #corrDf = pd.DataFrame(corrMatrix, index = rows, columns = columns)
        CorrMatrixDF = write_table(corrMatrix, outputDir = os.getcwd(), title = '{}.correlationMatrix'.format(prefix), rows = asvs , columns = asvs, dataframe = True)
        #plot_heatmap(CorrMatrixDF, outputDir = os.getcwd(), vmin = -1, vmax = 1, center = 0, title = '{}.correlationMatrix'.format(prefix), legendtitle = 'Correlation', text = None, symmetric = True)
        return(corrMatrix)
    
    @staticmethod
    def make_global_samples_distribution(shannonIndex, nASVs, reads = 20000):
        """Calculate the abundances distribution based on a shannon Index (H). Formula approach"""  
        if (-2*(shannonIndex - np.log(nASVs))) < 0:
            write_logfile('warning', 'ABUNDANCE ESTIMATION SHANNON', -2*(shannonIndex - np.log(nASVs)))
            exit(-5)
        else:
            sd = math.sqrt(-2*(shannonIndex - np.log(nASVs)))
            mu = 0
            abundistribution = lognormal(mu, sd, nASVs) # NOT USE ROUND, IT PROVOKE THE APARITION OF 0 values, WHEN THIS 0 VALUES ARE IN MUNB OF THE QZINEGBIN, IT RETURNS 1e+100 value in the abundance table
            # Convert it to probabilities and expand to number of required reads
            abundistribution_normalized = abundistribution/sum(abundistribution)
            names = ['sp' + str(i)  for i in range(nASVs)]
            counts = np.random.choice(names, size = reads, replace = True, p = abundistribution_normalized)
            #Counter(counts) - equal to table in R
            #!!!!!!!!!!!!!!!!! cuidadin con los 0, hay que pensarlo bien!!!! si darles una abundancia o media mínimas o directamente quitarlas, si las quito, entonces baja el numero de ASVs pedidas
            #finalcounts = np.append(np.array(list(Counter(counts).values())), np.zeros(nASVs - len(Counter(counts))))
            finalcounts = np.append(np.array(list(Counter(counts).values())), np.full(nASVs - len(Counter(counts)), 0.0001))
            return(finalcounts)
        
    @staticmethod
    def ZINBD(nASVs, nSamples, CorrMatrix, shannonIndex, reads = 20000, pstr0 = 0.2, size = 1):
        """ Sampling from a 'multivariate' ZINBD """
        multinormal = sample_from_multivariate_normal(CorrMatrix, nSamples, mu = 0)
        #multiNormal.shape (3,10) rows-samples/cols-species
        # Calculate the value of the areas under individual normal distributions
        probNormal = probabilities_from_normal(multinormal)
        # Calculate the value that gives that area in the ZINBD
        abundanceTable = qzinegbin(probNormal, size = size, munb = Mock.make_global_samples_distribution(shannonIndex, nASVs, reads), pstr0 = pstr0, nVariables = nASVs)
        return(abundanceTable)
    
