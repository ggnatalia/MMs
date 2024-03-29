#!/usr/bin/env python3


# Basic
import argparse 
import os
import sys
import re # regular expressions
import random
import numpy as np
import pandas as pd
import shutil
from glob import glob
# Other
#from collections import defaultdict

# My modules
from libs.libs import * # verify
from libs.maths import *
from libs.plots import *
from libs.enviroClass import Enviro
from libs.seqsClass import Sequence
from libs.mockClass import Mock
from libs.sampleClass import Sample

# Plot
#import matplotlib.pyplot as plt
#import seaborn as sns


# Logging
import logging
import daiquiri

def parse_arguments():
    """Parse the command line arguments and return an object containing them"""
    # Required
    general = argparse.ArgumentParser(description = 'Process input files')
    general.add_argument( '-m', '--mockName', type = str, required = True, help = 'Mock name. Ex: mock1')
    general.add_argument( '-o', '--output', type = str, required = True, help = 'Output directory. Preferably, name of the environment. Ex: Aquatic')
    general.add_argument( '-s','--start', default = 1, type = int, help = 'SILVA alignment reference positions-START. Default 1. 1-based') 
    general.add_argument( '-e','--end', default = 50000, type = int, help = 'SILVA alignment reference positions-END. Default 50000. 1-based')
    general.add_argument( '--by_region', default = None, nargs = '+', help = 'File with defined regions to introduce point mutations or list of V1-V9 regions')
    general.add_argument( '--by_position', default = None, help = 'File with defined positions to introduce point mutations') 

    general.add_argument( '--region', default = '16S', help = 'Name of the studied region')    
    general.add_argument( '-H','--shannonIndex', required=True, type = float, help = 'ShannonIndex') 
    general.add_argument( '-N', '--nSamples', required=True, type = int, help = 'Number of samples')
    general.add_argument( '-r', '--reads', required = True, type = int, help = 'Number of reads')

    # Options for making your mock randomly
    general.add_argument( '--rank', type = str, default = 'phylum', help = 'Rank to subset taxa: phlyum, order, class, family, genus')
    general.add_argument( '-ASVsmean','--ASVsmean',type = int, default = 6, help = 'Mean of mutant ASV per silva sequence')
    general.add_argument( '-nASVs','--nASVs', type = int, required = False, help = 'Number of ASVs') 
    general.add_argument( '-env', '--enviro', type = str, default = False, help = 'Let the user simulate an environmental mock. Look refEnv for options')

    # Options for customizing more your mock community
    general.add_argument( '-tx', '--taxa', type = str, nargs = '+', default = [], help = 'List of taxa')
    general.add_argument( '-seqs', '--seqs', type = str, nargs = '+', default = [], help = 'List of sequences\' names') 
    general.add_argument( '--minseqs', type = int, default = 100, help = 'Minimun number of sequences to extract from DB') # Subsettting random from silva means that you can subset 1 seq, produce, two strains and you'll want 5, no possibility to reach that number. Repeat strain generation
    general.add_argument( '-txAbund', '--taxaAbund', type = int, nargs = '+', default = [], help = 'List of abundances of taxa')
    general.add_argument( '--inputfile', type = str, help = 'The user can provide an align file without creating it from scratch.')
    
    # Reference files: 
    general.add_argument( '-ref','--ref', type = str, default = '', help = 'SILVA alignment reference')
    general.add_argument( '-refTax','--refTax', type = str, default = '', help = 'SILVA alignment TAX reference')    
    general.add_argument('-refEnv', '--refEnviro', type = str , default = '', help = 'Environment reference')

    # Extra features:
    general.add_argument('--threshold', type = float, default = 0.97, help = 'Make & filter strains using distances. Without this flag, sequences can be identical in the studied region')
    general.add_argument('--cpus', default = 12, type = int, help = 'Number of threads')   
    general.add_argument('--force-overwrite', action = 'store_true', help = 'Force overwrite if the output directory already exists')
        
    
    # NanoSim for long reads
    general.add_argument( '--Sim', default = 'InSilicoSeqs', choices = ['InSilicoSeqs', 'NanoSim'], help = 'Choose read simulator: InSilicoSeqs or NanoSim')
    
    # InSilicoSeqs parameters: add insert size & read length?
    general.add_argument( '--ISSerrormodel', default = 'perfect', type = str, help = 'Mode to generate InSilicoSeqs')
    general.add_argument( '--ISSparams', default = [150, 200], type = int, nargs = '+', help = '(Read length, insert size')
    general.add_argument( '--ISSsequences_files', type = str, nargs = '+', default = [], help = 'List of files to obtain reads: sequences <projectName><mockName>.<sampleName>.sequences16S.fasta. Same order that paired abundance files')
    general.add_argument( '--ISSabundance_files', type = str, nargs = '+', default = [], help = 'List of files to obtain reads: abundances<projectName><mockName>.<sampleName>.abundances. Same order that paired sequence files')
    general.add_argument( '--repeat_ISS_autocomplete', action = 'store_true', help = 'If True use mock and samples directly from the directory, without writing one by one the files')
    
    # NanoSim
    general.add_argument( '--NSerrormodel', default = 'perfect', choices = ['perfect', 'metagenome'], help = 'Mode to generate NanoSim')
    general.add_argument( '--NSparams', default = [500, 50], type = int, nargs = '+', help = '(Maximum read length, Minimum read length')
    general.add_argument( '--repeat_NS_autocomplete', action = 'store_true', help = 'If True use mock and samples directly from the directory, without writing one by one the files')

    
    # Other customizable parameters
    general.add_argument('--alpha', default = 0.9, type = float, help = 'Correlation Matrix: Probability that a coefficient is zero. Larger values enforce more sparsity.')
    general.add_argument('--pstr0', default = 0.2, type = float, help = 'ZINBD: Probability of structure 0')
    general.add_argument('--size', default = 1, type = int, help = 'ZINBD: Size - dispersion of ZINBD')
    general.add_argument('--ambiguities', default = 500, type = int, help = 'Number of Ns at the beginning or the end of the sequence')
    
    # For testing
    general.add_argument( '--just_taxa_selection', action = 'store_true', help = 'If True: do the selection of the sequences and stop.')
    general.add_argument( '--repeat_previous_mock', action = 'store_true', help = 'If True: With a previous mock, repeat reads simulation using another sequencing simulator.')

    args = general.parse_args()
    return(args)


def main(args): 
    
    # SET GENERAL VARIABLES
    path = os.getcwd()
    projectPrefix = args.output
    projectPath = os.getcwd() + '/' + args.output
    mockPrefix = '{}.{}'.format(args.output, args.mockName)
    mockPath = projectPath + '/' + args.mockName
    mockName = args.mockName
    # SET LOGGING
    if args.force_overwrite:
        if os.path.isfile('{}.log'.format(mockPrefix)):
            os.remove('{}.log'.format(mockPrefix))
        else:
            write_logfile('info', 'GENERAL', 'Not found a previous logfile, maybe you\'re using --force-overwrite but it\'s the first time you run it?')
    daiquiri.setup(level = logging.INFO, outputs = (daiquiri.output.STDERR, daiquiri.output.File('{}.log'.format(mockPrefix))))

    write_logfile('info', 'GENERAL: script launch ', os.path.abspath(__file__) + ' ' + ' '.join(sys.argv[1:]))
    write_logfile('info', 'GENERAL: os.getcwd() ', os.getcwd())
    
    
    # SET MORE VARIABLES
    start = int(args.start - 1) # user uses 1-based, convert to 0 based: a='12345' to select all the sequence: a[1-1:5]='12345'
    end = int(args.end)
    if isinstance(args.by_region, str):
        by_region = read_mutation_regions(args.by_region)
    elif isinstance(args.by_region, list):
        by_region = load_mutations_regions(args.by_region)
    else:
        by_region = args.by_region
    
    if args.by_position:
        by_position = load_positions(args.by_position)
    else:
        by_position = args.by_position
    
    
    region = args.region
    alignment = [start, end, region]
    
    rank = rank2number(args.rank)
    
    shannon = args.shannonIndex
    S = int(args.nSamples)
    ASVsmean = int(args.ASVsmean)
    enviro = args.enviro
    
    taxa = args.taxa
    seqs = args.seqs
    minseqs = args.minseqs
    taxaAbund = args.taxaAbund
    inputfile = args.inputfile

    Sim = args.Sim
    makemocks_home = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
    
    if not args.ref or not args.refTax:
        makemocks_home = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
        DB = makemocks_home + '/DB/'
        ref  = DB + '/' + 'silva.nr_v138.align'
        refTax = DB + '/' + 'silva.nr_v138.tax'
        refEnviro = DB + '/' + 'SpeciesperEnviro.tsv'
    else:    
        ref  =  args.ref
        refTax = args.refTax
        refEnviro = args.refEnviro
    
    # Check if DB are present
    if not os.path.isfile(refTax) and not os.path.isfile(ref):
        write_logfile('warning', 'DB:', 'DB {} has not been found. Please be sure to run /path/to/MMs/bin/make_databases.py to configure databases'.format(refTax.replace('.tax', '')))
        exit(-1)
    
    # Extra features:
    threshold = 1-args.threshold
    cpus = args.cpus
    force_overwrite = args.force_overwrite
    # InSilicoSeqs options
    ISSerrormodel = args.ISSerrormodel
    ISSparams = tuple(args.ISSparams)
    ISSsequences_files  = args.ISSsequences_files
    ISSabundance_files = args.ISSabundance_files
    # NanoSim options
    NSparams = tuple(args.NSparams)
    NSerrormodel = args.NSerrormodel
    
    if args.repeat_ISS_autocomplete:
        ISSsequences_files = [f.split('/')[-1] for f in sorted(glob('{}/samples/{}.S_*.sequences16S.fasta'.format(mockPath, mockPrefix)))]
        ISSsequences_files.sort(key =lambda e : e.split('_')) # sort by sample name
        ISSabundance_files = [f.split('/')[-1] for f in sorted(glob('{}/samples/{}.S_*.abundances'.format(mockPath, mockPrefix)))]
        ISSabundance_files.sort(key =lambda e : e.split('_')) # sort by sample name
    
    if args.repeat_NS_autocomplete:
        NSsequences_files = [f.split('/')[-1] for f in sorted(glob('{}/samples/{}.S_*.NSgenomes'.format(mockPath, mockPrefix)))]
        NSsequences_files.sort(key =lambda e : e.split('_')) # sort by sample name
        NSabundance_files = [f.split('/')[-1] for f in sorted(glob('{}/samples/{}.S_*.NSabundances'.format(mockPath, mockPrefix)))]
        NSabundance_files.sort(key =lambda e : e.split('_')) # sort by sample name
        NSdl_files = [f.split('/')[-1] for f in sorted(glob('{}/samples/{}.S_*.NSdl'.format(mockPath, mockPrefix)))]
        NSdl_files.sort(key =lambda e : e.split('_')) # sort by sample name
    
    
    
    # Other customizable parameters
    reads = args.reads
    pstr0 = args.pstr0
    size = args.size
    alpha = args.alpha
    ambiguidities = args.ambiguities
    # OPTION 1: REPEAT ONLY THE READ GENERATION: sequences and abundances files are provided, run Insilico, and prepare files to smiTE. 
    
    if (Sim == 'InSilicoSeqs' and args.repeat_ISS_autocomplete) or (ISSsequences_files and ISSabundance_files and args.repeat_ISS_autocomplete): 
        if len(ISSsequences_files) == len(ISSabundance_files) and len(ISSsequences_files) > 0: 
            write_logfile('info', 'MOCK REPEAT', 'Assuming you have a previous mock and you ONLY want to repeat the reads generation')
            if os.path.isdir(mockPath):
                os.chdir(mockPath)
            else:
                write_logfile('error', 'MOCK REPEAT', '{} does not exist. Please be sure you are running this script outside your output directory.'.format(mockPath))
                exit(-2) # Exit -2: dir not found
            if not os.path.isdir(mockPath + '/samples/'):
                write_logfile('error', 'MOCK REPEAT', '{} does not exist. Please be sure you are running this script outside your output directory.'.format(mockPath + '/samples/'))
                exit(-2) # Exit -2: dir not found
            Mock.init_and_run_InSilico(  mockPrefix, ISSsequences_files, ISSabundance_files, alignment = alignment, cpus = cpus, reads = reads, Sim = Sim, ISSerrormodel = ISSerrormodel, ISSparams = ISSparams )
    elif Sim == 'NanoSim' and args.repeat_NS_autocomplete:
        if len(NSsequences_files) == len(NSabundance_files) and len(NSsequences_files) > 0: 
            write_logfile('info', 'MOCK REPEAT', 'Assuming you have a previous mock and you ONLY want to repeat the reads generation')
            if os.path.isdir(mockPath):
                os.chdir(mockPath)
            else:
                write_logfile('error', 'MOCK REPEAT', '{} does not exist. Please be sure you are running this script outside your output directory.'.format(mockPath))
                exit(-2) # Exit -2: dir not found
            if not os.path.isdir(mockPath + '/samples/'):
                write_logfile('error', 'MOCK REPEAT', '{} does not exist. Please be sure you are running this script outside your output directory.'.format(mockPath + '/samples/'))
                exit(-2) # Exit -2: dir not found
            Mock.init_and_run_NanoSim(  mockPrefix, NSsequence_files = NSsequences_files, NSabundance_files = NSabundance_files, NSdl_files = NSdl_files, alignment = alignment, cpus = cpus, reads = reads, Sim =  Sim, NSerrormodel = NSerrormodel, NSparams = NSparams)
    elif args.repeat_previous_mock:
        write_logfile('info', 'MOCK REPEAT', 'Assuming you have a previous mock and you ONLY want to repeat the reads generation but with a different simulator')
        if os.path.isdir(mockPath):
            os.chdir(mockPath)
            if not os.path.isfile('../{}.fasta'.format(projectPrefix)) or not os.path.isfile('checkDB/{}.raw.abundances_original.tsv'.format(mockPrefix)):
                write_logfile('error', 'MOCK REPEAT', '../{}.fasta and checkDB/{}.raw.abundances_original.tsv do not exist. Please be sure you are running this script outside your output directory and that you have these files from a previous mock community.'.format(projectPrefix, mockPrefix))
            else:
                Mock.init_from_previousmock(mockPrefix, sequence_file = '../{}.fasta'.format(projectPrefix), abun_table = 'checkDB/{}.raw.abundances_original.tsv'.format(mockPrefix), alignment = alignment, cpus = cpus, reads = reads, Sim = Sim, ISSerrormodel = ISSerrormodel,  ISSparams = ISSparams, NSerrormodel = NSerrormodel,  NSparams = NSparams, ambiguidities = ambiguidities)
        else:
            write_logfile('error', 'MOCK REPEAT', '{} does not exist. Please be sure you are running this script outside your output directory.'.format(mockPath))
            exit(-2) # Exit -2: dir not found
            if not os.path.isdir(mockPath + '/samples/'):
                write_logfile('error', 'MOCK REPEAT', '{} does not exist. Please be sure you are running this script outside your output directory.'.format(mockPath + '/samples/'))
            
    else: 
        
        # Create the project directory:
        try:
            os.mkdir(projectPath)
        except OSError as e: #[Errno 17] File exists: 'output'
            if e.errno != 17:
                raise
            elif not force_overwrite and  inputfile:
                write_logfile('warning', 'OUTPUT DIRECTORY GENERATION', 'The user wants to use a previous file: {}'.format(inputfile))
                pass # User wants to use a previous align
            elif not force_overwrite and not inputfile:
                write_logfile('warning', 'OUTPUT DIRECTORY GENERATION', 'The directory {} already exists. Please, remove it or choose other name for the output directory'.format(projectPath))
                exit(-17)
            else:
                write_logfile('warning', 'OUTPUT DIRECTORY GENERATION', 'You ran --force-overwrite. The previous directory {} will be deleted and a new one will be created'.format(projectPath))
                shutil.rmtree(projectPath, ignore_errors = True)
                os.mkdir(projectPath)
    
        os.chdir(projectPath)
        if inputfile:
            Env = Enviro.init_from_file(prefix = projectPrefix, path = projectPath, inputfile = inputfile, cpus = cpus)
            alignment = [] # No tengo las secuencias alineadas
        else: # No align, do mutant ASVs
            if not args.nASVs:
                write_logfile('warning', 'nASVs parameter', 'Please, as you\'re not using a previous input file, indicate the number of ASVs to be simulated')
                exit(-1)
            nASVs = int(args.nASVs)
            if enviro:
                Env = Enviro.init_from_enviro(nASVs = nASVs, prefix = projectPrefix, enviro = enviro, refEnviro = refEnviro, refTax = refTax, ref = ref, nTaxa = 10000, rank = 5, cpus = cpus, by_region = by_region)
                rank = 5
            elif taxa:
                if os.path.isfile('{}/{}'.format(path,taxa[0])):
                    taxa, taxaAbund = parseTaxas('{}/{}'.format(path,taxa[0]))
                print('\t'.join(taxa))
                print('\t'.join(map(str,taxaAbund)))
                Env = Enviro.init_from_taxa(nASVs = nASVs, prefix = projectPrefix, rank = rank, taxa = taxa, taxa_abundances = taxaAbund, refTax =  refTax, ref = ref, cpus = cpus)
            else: # seqs
                Env = Enviro.init_from_seqs(prefix = projectPrefix, rank = rank, seqs = seqs, nASVs = nASVs, minseqs = minseqs, refTax =  refTax, ref = ref, cpus = cpus)
            write_logfile('info', 'ENVIRONMENT', 'Simulating microdiversity')
            Env.makeASVs(region, start, end, ASVsmean, threshold , cpus, by_region = by_region, by_pos = by_position)   # Only with sequences that are not in the align file. Assume align file provided by the user is ok!
            write_logfile('info', 'ENVIRONMENT', 'Writing output files')
            Env.write_output()
            write_logfile('info', 'ENVIRONMENT', 'Plotting taxa')
            Env.plot_taxonomy(rank = rank)
        
        if not args.just_taxa_selection:
            # 2_CREATE THE MOCK!!
            write_logfile('info', 'MOCK GENERATION', os.getcwd())
            try:
                os.mkdir(mockPath)
            except OSError as e: #[Errno 17] File exists: 'mockDir'
                if e.errno != 17:
                    raise
                else:    
                    write_logfile('warning', 'MOCK GENERATION', 'The directory {} already exists. Please, remove it or choose other name for the output directory'.format(mockPath))
                    write_logfile('info', 'MOCK GENERATION', 'To avoid repeating the environment generation, provide a align with the align flag and repeat only the mock generation with othe mock name')
                    exit(-17)
            
            os.chdir(mockPath)
            os.mkdir('{}/samples'.format(mockPath))
            os.mkdir('{}/checkDB'.format(mockPath))
            #Define mock abundances
            mock = Mock.init_from_enviro(Env, mockPrefix, S, shannon, alignment = alignment, cpus = cpus, alpha = alpha,  smallest_coef = 0.1, largest_coef = 0.9, reads = reads, pstr0 = pstr0, size = size, Sim = Sim, ISSerrormodel =ISSerrormodel, ISSparams = ISSparams, NSerrormodel = NSerrormodel, NSparams = NSparams, ambiguidities = ambiguidities)
            os.chdir(path)
            # ShannonIndex correspondence (mothur is REQUIRED : CHECK)
            write_logfile('info', 'SHANNON INDEX', 'To know the correspondence of Shannon Index in the different taxonomic ranks, mothur is required')
            command = ['{}/utils/shannonIndex_sweep.py'.format(makemocks_home), '-o', projectPrefix, '-m', mockName, '--align', '{}.align'.format(projectPrefix)] # in path in conda. Within utils in the github repo
            runCommand(command)
        
################################################################################################################
    
if __name__ == '__main__':
    main(parse_arguments())





