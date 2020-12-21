#!/usr/bin/env python3

import pandas as pd
import numpy as np
import random
import subprocess
import logging
import daiquiri
################################################################################################### Useful functions, common in different steps

########################### WORK WITH FASTA/ ALIGN FILES

def simplifyString(string, splitChar = '\t', conservative = False):
    """ For a string: change splitChar by '_' or split by it """
    if conservative:# Replace splitChar by '_' and remove blank spaces
        string = string.replace(splitChar, '_').replace(' ','')
    else:# Maintain first element of the splitted string
        if splitChar in string:  
            string = string.split(splitChar)[0]
    return(string)

def fasta2dict(fasta):
        SeqDic = {}
        with open(fasta) as f:
            for sequece in f.read().strip().lstrip('>').split('>'):
                name, seq = sequece.split('\n',1)
                #print(name)
                name = name.split('\t')[0]
                if not 'N' in seq:
                    SeqDic[name] = seq.replace('\n','')
        return(SeqDic)


def fastq2fasta(fastq, fasta):
    """ Convert fastq files to fasta file. Ignoring qualities True """
    with open(fastq, 'r') as fastqfile, open('{}.fasta'.format(fasta), 'w') as fastafile:
        while True:
            header = fastqfile.readline().rstrip('\n').lstrip('@')
            if not header:
                break
            seq = fastqfile.readline().rstrip('\n')
            coment= fastqfile.readline().rstrip('\n')
            qual = fastqfile.readline().rstrip('\n')
            fastafile.write('>{}\n{}\n'.format(header, seq))

########################### LIST

def flattened(l):
    """ Flatten a list """
    return([y for x in l for y in x])


def reversedict(d):
    """ Convert keys into values and vice versa. Values are returned in lists for repeated values """
    # https://stackoverflow.com/questions/1031851/how-do-i-exchange-keys-with-values-in-a-dictionary
    inverse = {}
    for i in d:  
        inverse.setdefault(d[i],[]).append(i)
    return(inverse)



########################### WORK WITH TABLES
def make_percent(data, outputDir = '.', write = True, fileName = None, T = False):
    """ Calculate the percentages %1 by rows of an abudance table """
    if T:
        data = data.T
    pcdata = data[data.columns].apply(lambda x: x/x.sum(), axis=1)
    if T:
        pcdata = pcdata.T
    if write:
        write_table(pcdata, title = '{}/{}'.format(outputDir, fileName), rows = list(pcdata.index) , columns = list(pcdata.columns), dataframe = None)
    return(pcdata)

def write_table(data, title, outputDir = '.', rows = None , columns = None, dataframe = None):
    """ Write a csv file from matrix. Return a pandas Dataframe if asked """
    if not rows:
        rows = ['R_' + str(i) for i in range(data.shape()[0])]
    if not columns:
        columns = ['C_' + str(i) for i in range(data.shape()[1])]
    df = pd.DataFrame(data, index = rows , columns = columns)
    export_csv = df.to_csv ('{}/{}.tsv'.format(outputDir, title), index = True, header=True, sep = '\t')
    export_csv = df.to_csv ('{}/{}.csv'.format(outputDir, title), index = True, header=True, sep = ',')
    if dataframe:
        return(df)
    
def load_table(filepath, rows = None , cols = None, path = '.', sep = ','):
    """ Load a csv file as pandas dataframe """
    if rows  and cols:
        data = pd.read_csv(filepath, sep = sep, index_col = 0)
    else:
        data = pd.read_csv(filepath, sep = sep, header=None)
    return(data)

########################### MUTATIONS
def mutate(string, N, start = None, end = None, randomly = True): 
    """ Generate a mutate string from original one in N positions """
    bases = ['A','T','C','G','-']
    strings_list = np.array(list(string), dtype ='U1')
    if not start: start = 0
    if not end: end = len(string)-1
    # If mutation site includes '.' sample another mutation site. No repeat positions
    # List of possible positions without '.'
    possiblePos = [i for i, nt in enumerate(strings_list) if nt != '.']
    mutation_sites = list(np.random.choice(possiblePos, N, replace = False))
    for pos in mutation_sites: # Avoid cases in which the mutated base is identical to the original
        strings_list[pos] = random.sample(list(filter(lambda x: x != string[pos], bases)), 1)[0]
    
    return (''.join(strings_list.tolist()))

########################### TAXONOMY
def rank2number(rank): 
    """ Convert taxonomic ranks into numerical factors 0 based """
    if rank == 'kingdom':
        return(0)
    elif rank == 'phylum':
        return(1)
    elif rank == 'class':
        return(2)
    elif rank == 'order':
        return(3)
    elif rank == 'family':
        return(4)
    elif rank == 'genus':
        return(5)
    else:
        write_logfile('error', 'LIBS: rank2number' , '{} is not a valid taxonomic rank. Valid taxonomic ranks are: phylum, class, order, family, genus'.format(rank))
        exit(0)


def estimate_mutations(rank, length = 1500): #1500 bp approximately 16S rRNA
    if rank == 1: #phylum
        cutoff =  1-0.75
        write_logfile( 'info', 'FAKE TAXA', '{} cutoff {}'.format(rank, cutoff) )
        return(cutoff * length)
    elif rank == 2: #class
        cutoff =  1-0.785
        write_logfile( 'info', 'FAKE TAXA', '{} cutoff {}'.format(rank, cutoff) )
        return(cutoff * length)
    elif rank == 3: #order
        cutoff =  1-0.82
        write_logfile( 'info', 'FAKE TAXA', '{} cutoff {}'.format(rank, cutoff) )
        return(cutoff * length)
    elif rank == 4: #family
        cutoff =  1-0.865
        write_logfile( 'info', 'FAKE TAXA', '{} cutoff {}'.format(rank, cutoff) )
        return(cutoff * length)
    elif rank == 5: #genus
        cutoff =  1-0.945
        write_logfile( 'info', 'FAKE TAXA', '{} cutoff {}'.format(rank, cutoff) )
        return(cutoff * length)
    else:
        return(100)


def loadTaxa(refTax = '/home/natalia/Projects/natalia/DB/silva.nr_v138/silva.nr_v138.tax', rank = None):
    """ For a seq return taxonomy. Prepare for parsing silvaDB taxonomy"""
    with open(refTax, 'r') as taxReference:
        all_sequences_list = taxReference.readlines()
        if not rank:
            taxas = {line.rstrip('\n').split('\t')[0] : line.rstrip('\n').split('\t')[1] for line in all_sequences_list}
        else:
            taxas = {line.rstrip('\n').split('\t')[0] : line.rstrip('\n').split('\t')[1].split(';')[rank] for line in all_sequences_list}
        return(taxas)
    
def loadEnviro(referenceEnviro = '/media/natalia/Linux/natalia/opt/makemocks/DB/speciesperEnvironment.tsv'):
    """ For the environment table subset an enviro """
    enviro_df = pd.read_csv(referenceEnviro, sep = '\t', comment = '#', index_col = 0) # OK
    return(enviro_df)
    
########################### SYSTEM

def runCommand(command, stdout = None, stderr = None):
    """Run the command and check the success of the subprocess. Return exit if it went wrong"""    
    write_logfile('info', 'Command' , '{}'.format(' '.join(map(str, command))))
    exitcode = subprocess.call(map(str, command), stdout = stdout, stderr = stderr)
    if exitcode != 0:
        write_logfile('error', 'LIBS: runCommand' , 'There must be some problem with "{}".\nIt\'s better to stop and check it'.format(' '.join(map(str, command))))
        exit(-1)

def write_logfile(level, step, message):
    """ Write log file using daiquiri and logging """
    if level == 'debug':
        daiquiri.getLogger(step).debug(message)
    elif level == 'info':
        daiquiri.getLogger(step).info(message)
    elif level == 'warning':
        daiquiri.getLogger(step).warning(message)
    elif level == 'error':
        daiquiri.getLogger(step).error(message)
    else:
        daiquiri.getLogger(step).critical(message)



    
    
    