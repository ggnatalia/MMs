#!/usr/bin/env python3


# Basic
import argparse 
import os
import sys
import re 
import random
import numpy as np
import pandas as pd
import shutil
from glob import glob

# My modules
from libs.libs import * 
from libs.maths import *
from libs.plots import *
from libs.enviroClass import Enviro
from libs.seqsClass import Sequence
from libs.mockClass import Mock
from libs.sampleClass import Sample

# Logging
import logging
import daiquiri

def parse_arguments():
    """Parse the command line arguments and return an object containing them"""
    
    # Required
    general = argparse.ArgumentParser(description = 'Process input files')
    
    general.add_argument( '-m', '--mockName', type = str, required = True, help = 'Name of the mock e.g. mock1.' )
    general.add_argument( '-o', '--output', type = str, required = True, help = 'Name of the output directory e.g. aquatic.' )
    general.add_argument( '-N', '--nSamples', type = int, default = 5, help = 'Number of samples to generate. Default 5 samples per mock.' )
    general.add_argument( '-nASVs','--nASVs', type = int, help = 'Number of ASVs.' ) 

    # Options to design your mock community
    general.add_argument( '-H','--shannon', type = float, default = 3, help = 'Shannon diversity Index. Default 3.' )
    general.add_argument( '-r','--rank', type = str, default = 'genus', help = 'Rank to subset taxa: phlyum, order, class, family, genus.' )
    general.add_argument( '-ASVsmean','--ASVsmean', type = int, default = 5, help = 'Mean of the number of mutants per reference.' )
    general.add_argument( '-env', '--enviro', type = str, default = False, help = 'Select one of the predefined environments to simulate the mock. See the file SpeciesperEnviro.tsv for options.' )

    # Options for customizing more your mock community
    general.add_argument( '--taxa', type = str, nargs = '+', default = [], help = 'List of taxa to generate the mock.' )
    general.add_argument( '--seqs', type = str, nargs = '+', default = [], help = 'List of sequences\'s header.' ) 
    general.add_argument( '--minseqs', type = int, default = 5, help = 'Minimun number of sequences to randomly extract from DB. Default 5.' ) 
    general.add_argument( '--taxaAbund', type = int, nargs = '+', default = [], help = 'List of taxa abundances in percentages.' )
    general.add_argument( '--input-file', type = str, help = 'Provide a previous fasta/align file.' )
    
    # Reference files
    general.add_argument( '-ref','--ref', type = str, help = 'SILVA alignment reference formatted for using mothur e.g. silva.nr_v138.align.' )
    general.add_argument( '-refTax','--refTax', type = str, help = 'SILVA taxonomy reference formatted for using mothur e.g. silva.nr_v138.tax.' )    
    general.add_argument( '-refEnv', '--refEnviro', type = str, help = 'Table with most frequent species per environment e. g. SpeciesperEnviro.tsv.' )
    # SILVA
    general.add_argument( '-s','--start', type = int, default = 1,  help = 'Start position in the SILVA alignment reference. Default 1.' ) 
    general.add_argument( '-e','--end', type = int, default = 50000, help = 'End position in the SILVA alignment reference. Default 50000.' ) 
    general.add_argument( '--region', type = str, default = '16S', help = 'Region of the 16S rRNA gene.' )

    # Extra features:
    general.add_argument( '--cutoff', type = float, default = 0.03, help = 'Filter ASVs depending on sequence similarity. Without this flag, sequences can be identical in the studied region.' )
    general.add_argument( '--cpus', type = int, default = 12, help = 'Number of threads.' )   
    general.add_argument( '--force-overwrite', action = 'store_true', help = 'Force overwrite if the output directory already exists.' )
    
    # InSilicoSeqs parameters: add insert size & read length
    general.add_argument( '--error-model', type = str, default = 'perfect', help = 'Mode to simulate reads in InSilicoSeqs.' )
    general.add_argument( '--read-length', type = int, nargs = '+', default = 150, help = 'Read length of the simulated reads. Only available for basic and perfect modes.' )
    general.add_argument( '--insert-size', type = int, nargs = '+', default = 200, help = 'Insert size of the simulated reads. Only available for basic and perfect modes.' )
    general.add_argument( '--sequences-files', type = str, nargs = '+', default = [], help = 'Sequences files for InSilicoSeqs: <projectName><mockName>.<sampleName>.sequences16S.fasta.' )
    general.add_argument( '--abundance-files', type = str, nargs = '+', default = [], help = 'Abundances files for InSilicoSeqs: <projectName><mockName>.<sampleName>.abundances.' )
    general.add_argument( '--repeat-inSilicoSeqs-autocomplete', action = 'store_true', help = 'Repeat reads simulation using previous files.' )
    
    # Other customizable parameters
    general.add_argument( '--reads', type = int, default = 50000, help = 'Number of reads approximately in total, counting both pairs.' )
    general.add_argument( '--alpha', type = float, default = 0.9, help = 'Correlation Matrix: Probability that a coefficient is zero. Larger values enforce more sparsity.' )
    general.add_argument( '--pstr0', type = float, default = 0.2, help = 'ZINBD: Probability of structure 0.' )
    general.add_argument( '--size', type = int, default = 1, help = 'ZINBD: Size - dispersion of ZINBD.' )
    # Plots
    general.add_argument( '--figsize', type = tuple, default = (20, 20), help = 'Size of plots.' )
    
    args = general.parse_args()
    return(args)

def main(args): 
    """ Set variables, set the way: environment, random, list of taxa """
    
    # SET GENERAL VARIABLES
    MMs_home = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
    DB = MMs_home + '/DB'
    
    path = os.getcwd()
    projectPrefix = args.output
    projectPath = os.getcwd() + '/' + args.output
    mockName = args.mockName
    mockPrefix = '{}.{}'.format(args.output, mockName)
    mockPath = projectPath + '/' + mockName
    
    # SET LOGGING
    if args.force_overwrite:
        os.remove('{}.log'.format(mockPrefix))
    daiquiri.setup(level = logging.INFO, outputs = (daiquiri.output.STDERR, daiquiri.output.File('{}.log'.format(mockPrefix))))

    write_logfile( 'info', 'GENERAL: Launch ', os.path.abspath(__file__) + ' ' + ' '.join(sys.argv[1:]) )
    write_logfile( 'info', 'GENERAL: Location ', os.getcwd() )
    
    
    # SET MORE VARIABLES
    if args.ref:
        ref  =  args.ref
    else:
        ref = DB + '/silva.nr_v138.align'        
    if args.refTax:
        refTax = args.refTax
    else:
        refTax = DB + '/silva.nr_v138.tax'
    if args.refEnviro:
        refEnviro = args.refEnviro
    else:
        refEnviro = DB + '/SpeciesperEnviro.tsv'

    if args.start == 0:
        sys.exit('M&Ms uses 1-based start position. To refer to the first position, fix args.start = 1, we\'ll satisfy python indexes for you.')
    else:
        start = args.start - 1
    end = args.end
    region = args.region
    alignment = [start, end, region]
    rank = rank2number(args.rank)
    
    shannon = args.shannon
    S = args.nSamples
    ASVsmean = args.ASVsmean
    
    enviro = args.enviro
    taxa = args.taxa
    seqs = args.seqs
    minseqs = args.minseqs
    taxaAbund = args.taxaAbund
    inputfile = args.input_file

    # Extra features:
    reads = args.reads
    pstr0 = args.pstr0
    size = args.size
    alpha = args.alpha
    figsize = args.figsize
    cutoff = args.cutoff
    cpus = args.cpus
    force_overwrite = args.force_overwrite
    errormodel = args.error_model
    inSilicoparams = tuple(args.read_length, args.insert_size)
    sequences_files  = args.sequences_files
    abundance_files = args.abundance_files
      
    if args.repeat_inSilicoSeqs_autocomplete:
        sequences_files = [f.split('/')[-1] for f in sorted(glob('{}/samples/{}.S_*.sequences16S.fasta'.format(mockPath, mockPrefix)))]
        sequences_files.sort(key =lambda e : e.split('_')) # sort by sample name
        abundance_files = [f.split('/')[-1] for f in sorted(glob('{}/samples/{}.S_*.abundances'.format(mockPath, mockPrefix)))]
        abundance_files.sort(key =lambda e : e.split('_')) # sort by sample name
    
    # OPTION 1: REPEAT READ GENERATION: sequences and abundances files are provided. 
    if len(sequences_files) == len(abundance_files) and len(sequences_files) > 0: 
        write_logfile('info', 'MOCK REPEAT', 'M&Ms will only repeat the reads generation.')
        if os.path.isdir(mockPath):
            os.chdir(mockPath)
        else:
            write_logfile('error', 'REPEAT READS SIMULATION', '{} doesn\'t exist. Remember to launch M&Ms outside your output directory.'.format(mockPath))
            exit(2) # Exit -2: dir not found
        if not os.path.isdir(mockPath + '/samples/'):
            write_logfile('error', 'REPEAT READS SIMULATION', '{} doesn\'t exist. Remember to launch M&Ms outside your output directory. Within your output directory, should be a mock folder with a checkDB & samples directories'.format(mockPath + '/samples/'))
            exit(2) # Exit -2: dir not found
        Mock.init_and_run_inSilico( mockPrefix, sequences_files = sequences_files, abundance_files = abundance_files, errormodel = errormodel, alignment = alignment, reads = reads, inSilicoparams = inSilicoparams, figsize = figsize, cpus = cpus )
    else: 
        # Create the project directory
        try:
            os.mkdir(projectPath)
        except OSError as e:
            if e.errno != 17:  #[Errno 17] File exists: 'output'
                raise
            elif not force_overwrite and inputfile:
                write_logfile('info', 'OUTPUT DIRECTORY', 'M&Ms\'ll use your file {}'.format(inputfile))
                pass
            elif not force_overwrite and not inputfile:
                write_logfile('warning', 'OUTPUT DIRECTORY', 'The directory {} already exists. Please, remove it or choose other name for the output directory'.format(projectPath))
                exit(17)
            else:
                write_logfile('warning', 'OUTPUT DIRECTORY', 'You activated the flag --force-overwrite. The directory {} will be removed and a new one will be created'.format(projectPath))
                shutil.rmtree(projectPath, ignore_errors = True)
                os.mkdir(projectPath)
    
        os.chdir(projectPath)
        if inputfile:
            env = Enviro.init_from_file(prefix = projectPrefix, path = projectPath, inputfile = inputfile, refTax = refTax)
            alignment = [] # If fasta, there aren't align sequences
        else: # No align, do mutant ASVs
            nASVs = int(args.nASVs)
            if enviro:
                rank = 5
                env = Enviro.init_from_enviro(nASVs = nASVs, prefix = projectPrefix, rank = rank, enviro = enviro, refEnviro = refEnviro, refTax = refTax, ref = ref, nTaxa = 10000)
            elif taxa: 
                Env = Enviro.init_from_taxa(nASVs = nASVs, prefix = projectPrefix, rank = rank, taxa = taxa, taxa_abundances = taxaAbund, refTax =  refTax, ref = ref)
            else: # seqs
                env = Enviro.init_from_seqs(nASVs = nASVs, prefix = projectPrefix, rank = rank, seqs = seqs,  minseqs = minseqs, refTax =  refTax, ref = ref)
            
            env.makeASVs(region = region, start = start, end = end, ASVsmean = ASVsmean, cutoff = cutoff, cpus = cpus, figsize = figsize)
            write_logfile('info', 'ENVIRONMENT', 'Writing output files')
            env.write_output()
            write_logfile('info', 'ENVIRONMENT', 'Plotting taxa')
            env.plot_taxonomy(rank = rank)
        
        # 2_CREATE THE MOCK
        write_logfile('info', 'MOCK GENERATION', 'Sequences have been chosen! Alea iacta est.')
        try:
            os.mkdir(mockPath)
        except OSError as e: #[Errno 17] File exists: 'mockDir'
            if e.errno != 17:
                raise
            else:    
                write_logfile('warning', 'MOCK', 'The directory {}/{} already exists. Please, remove it or choose other name for the output directory'.format(mockPath, mockPrefix))
                write_logfile('info', 'MOCK', 'To not repeat sequences selection, provide a fasta/align using the flag --input-file and create a new mock using those sequences')
                exit(-17)
        
        os.chdir(mockPath)
        os.mkdir('{}/samples'.format(mockPath))
        os.mkdir('{}/checkDB'.format(mockPath))
        
        mock = Mock.init_from_enviro(enviro = env, prefix = mockPrefix, nSamples = S, shannon = shannon, alignment = alignment, alpha = alpha, smallest_coef = smallest_coef, largest_coef = largest_coef, reads = reads, pstr0 = pstr0, size = size, figsize = figsize, inSilicoparams = inSilicoparams)
        os.chdir(path)
        command = ['python3', '{}/utils/shannonIndex_sweep.py'.format(MMs_home), '-o', projectPrefix, '-m', mockName, '--align', '{}.align'.format(projectPrefix)]
        print(' '.join(command))
        runCommand(command)
        
        # Taxonomy subset check 
        if enviro:
            command = ['Rscript', '{}/check_scripts/checkTaxonomy.R'.format(os.path.dirname(os.path.abspath(__file__))), projectPrefix, projectPath, mockPrefix, mockPath, enviro, rank, 15 ]
            #runCommand(command)
################################################################################################################
    
if __name__ == '__main__':
    main(parse_arguments())
