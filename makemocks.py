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
    general.add_argument( '-H','--shannonIndex', type = float, default = 3, help = 'Shannon diversity Index. Default 3.' )
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
    general.add_argument( '--repeat-InSilicoSeqs-autocomplete', action = 'store_true', help = 'Repeat reads simulation using previous files.' )
    
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
    
    MMs_home = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
    DB = MMs_home + '/DB'
    
    # SET GENERAL VARIABLES
    path = os.getcwd()
    projectPrefix = args.output
    projectPath = os.getcwd() + '/' + args.output
    mockPrefix = '{}.{}'.format(args.output, args.mockName)
    mockPath = projectPath + '/' + args.mockName
    mockName = args.mockName
    # SET LOGGING
    if args.force_overwrite:
        os.remove('{}.log'.format(mockPrefix))
    daiquiri.setup(level = logging.INFO, outputs = (daiquiri.output.STDERR, daiquiri.output.File('{}.log'.format(mockPrefix))))

    write_logfile('info', 'GENERAL: script launch ', os.path.abspath(__file__) + ' ' + ' '.join(sys.argv[1:]))
    write_logfile('info', 'GENERAL: os.getcwd() ', os.getcwd())
    
    
    # SET MORE VARIABLES
    start = int(args.start - 1) # user uses 1-based, convert to 0 based: a='12345' to select all the sequence: a[1-1:5]='12345'
    end = int(args.end)
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
        
    

    # Extra features:
    cutoff = args.cutoff
    cpus = args.cpus
    force_overwrite = args.force_overwrite
    errormodel = args.error_model
    InSilicoparams = tuple(args.read_length, args.insert_size)
    sequences_files  = args.sequences_files
    abundance_files = args.abundance_files
    
    if args.repeat_InSilicoSeqs_autocomplete:
        sequences_files = [f.split('/')[-1] for f in sorted(glob('{}/samples/{}.S_*.sequences16S.fasta'.format(mockPath, mockPrefix)))]
        sequences_files.sort(key =lambda e : e.split('_')) # sort by sample name
        abundance_files = [f.split('/')[-1] for f in sorted(glob('{}/samples/{}.S_*.abundances'.format(mockPath, mockPrefix)))]
        abundance_files.sort(key =lambda e : e.split('_')) # sort by sample name
    
    
    
    # Other customizable parameters
    reads = args.reads
    pstr0 = args.pstr0
    size = args.size
    alpha = args.alpha
    figsize = args.figsize
    # OPTION 1: REPEAT ONLY THE READ GENERATION: sequences and abundances files are provided, run Insilico, and prepare files to smiTE. 
    
    if len(sequences_files) == len(abundance_files) and len(sequences_files) > 0: 
        write_logfile('info', 'MOCK REPEAT', 'Assuming you have a previous mock and you ONLY want to repeat the reads generation')
        if os.path.isdir(mockPath):
            os.chdir(mockPath)
        else:
            write_logfile('error', 'MOCK REPEAT', '{} does not exist. Please be sure you are running this script outside your output directory.'.format(mockPath))
            exit(-2) # Exit -2: dir not found
        if not os.path.isdir(mockPath + '/samples/'):
            write_logfile('error', 'MOCK REPEAT', '{} does not exist. Please be sure you are running this script outside your output directory.'.format(mockPath + '/samples/'))
            exit(-2) # Exit -2: dir not found
        #Mock.reads_generation(prefix = mockPrefix, sequences_files = sequences_files, abundance_files = abundance_files, errormodel = errormodel, reads = reads, cpus = cpus)
        Mock.init_and_run_InSilico(  mockPrefix, sequences_files, abundance_files, errormodel = errormodel, alignment = alignment, reads = reads, InSilicoparams = InSilicoparams, figsize = figsize, cpus = cpus)
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
            Env = Enviro.init_from_file(prefix = projectPrefix, path = projectPath, inputfile = inputfile, refTax = refTax)
            alignment = [] # No tengo las secuencias alineadas
        else: # No align, do mutant ASVs
            nASVs = int(args.nASVs)
            if enviro:
                Env = Enviro.init_from_enviro(nASVs = nASVs, prefix = projectPrefix, enviro = enviro, refEnviro = refEnviro, refTax = refTax, ref = ref, nTaxa = 10000, rank = 5)
                rank = 5
            elif taxa: 
                Env = Enviro.init_from_taxa(nASVs = nASVs, perfix = projectPrefix, rank = rank, taxa = taxa, taxa_abundances = taxaAbund, refTax =  refTax, ref = ref)
            else: # seqs
                Env = Enviro.init_from_seqs(prefix = projectPrefix, rank = rank, seqs = seqs, nASVs = nASVs, minseqs = minseqs, refTax =  refTax, ref = ref)
            Env.makeASVs(region, start, end, ASVsmean, cutoff , cpus, figsize)   # Only with sequences that are not in the align file. Assume align file provided by the user is ok!
            write_logfile('info', 'ENVIRONMENT', 'Writing output files')
            Env.write_output()
            write_logfile('info', 'ENVIRONMENT', 'Plotting taxa')
            Env.plot_taxonomy(rank = rank)
            #write_logfile('warning', 'ENVIRONMENT', 'Plotting taxon entropies: USE A LOT OF RAM')
            #Env.plot_taxon_entropy(ref = ref, refTax = refTax, rank = rank)
        
        # 2_CREATE THE MOCK!!
        write_logfile('info', 'MOCK GENERATION', os.getcwd())
        try:
            os.mkdir(mockPath)
        except OSError as e: #[Errno 17] File exists: 'mockDir'
            if e.errno != 17:
                raise
            else:    
                write_logfile('warning', 'MOCK GENERATION', 'The directory {}/{} already exists. Please, remove it or choose other name for the output directory'.format(mockPath, mockPrefix))
                write_logfile('info', 'MOCK GENERATION', 'To avoid repeating the environment generation, provide a align with the align flag and repeat only the mock generation with othe mock name')
                exit(-17)
        
        os.chdir(mockPath)
        os.mkdir('{}/samples'.format(mockPath))
        os.mkdir('{}/checkDB'.format(mockPath))
        #Define mock abundances
        mock = Mock.init_from_enviro(Env, mockPrefix, S, shannon, alignment = alignment, alpha = alpha,  smallest_coef = 0.1, largest_coef = 0.9, reads = reads, pstr0 = pstr0, size = size, figsize = figsize, InSilicoparams = InSilicoparams)
        #mock.run_mock()
        os.chdir(path)
        # ShannonIndex correspondence (mothur is REQUIRED : CHECK)
        write_logfile('info', 'SHANNON INDEX', 'To know the correspondence of Shannon Index in the different taxonomic ranks, mothur is required')
        command = ['python3', '{}/utils/shannonIndex_sweep.py'.format(os.path.dirname(os.path.abspath(__file__))), '-o', projectPrefix, '-m', mockName, '--align', '{}.align'.format(projectPrefix)]
        runCommand(command)
        # Taxonomy subset check 
        if enviro:
            command = ['Rscript', '{}/check_scripts/checkTaxonomy.R'.format(os.path.dirname(os.path.abspath(__file__))), projectPrefix, projectPath, mockPrefix, mockPath, enviro, rank, 15 ]
            runCommand(command)
################################################################################################################
    
if __name__ == '__main__':
    main(parse_arguments())


#refTax = '/home/natalia/Projects/natalia/DB/silva.nr_v138/silva.nr_v138.tax'
#newSeqs = Sequence.set_sequences(fastaFile = 'test/Saline/Saline.align', refTax = '/home/natalia/Projects/natalia/DB/silva.nr_v138/silva.nr_v138.tax', degap = False, cleanHeader = False, splitChar = '\t', conservative = False)


