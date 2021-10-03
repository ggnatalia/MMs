#!/usr/bin/env python

import pandas as pd
import numpy as np
import random
import subprocess
import logging
import daiquiri
import os
import shutil
################################################################################################### Useful functions, common in different steps

########################### WORK WITH FASTA/ ALIGN FILES

def simplifyString(string, splitChar = '\t', conservative = False): # TESTED!
    """ For a string: change splitChar by '_' or split by it """
    if conservative:# Replace splitChar by '_' and remove blank spaces
        string = string.replace(splitChar, '_').replace(' ','')
    else:# Maintain first element of the splitted string
        if splitChar in string:  
            string = string.split(splitChar)[0]
    return(string)

def fasta2dict(fasta): # if fasta is huge, it consumes a lot of RAM memory
    SeqDic = {}
    with open(fasta) as f:
        for sequece in f.read().strip().lstrip('>').split('>'):
            name, seq = sequece.split('\n',1)
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
    #export_csv = df.to_csv ('{}/{}.csv'.format(outputDir, title), index = True, header=True, sep = ',')
    if dataframe:
        return(df)

def read_table(table):
    df = pd.read_csv(table, sep = '\t', header=0, index_col = 0)
    return(df)



def load_table(filepath, rows = None , cols = None, path = '.', sep = ','):
    """ Load a csv file as pandas dataframe """
    if rows  and cols:
        data = pd.read_csv(filepath, sep = sep, index_col = 0)
    else:
        data = pd.read_csv(filepath, sep = sep, header=None)
    return(data)

########################### MUTATIONS
def mutate(string, N, start = None, end = None, regions = None, header = None, positions = None): # WORK HACER QUE FUNCIONE QUITANDO LOS GAPS A LA SECUENCIA 
    """ Generate a mutate string from original one in N positions """
    bases = ['A','T','C','G','-']
    if not start:
        start = 0
    if not end:
        end = len(string) - 1
    string = string[start:end]
    original_string = np.array(list(string), dtype ='U1')
    string_degap = string.replace('-','').replace('.','')
    strings_list = np.array(list(string_degap), dtype ='U1')
    pos_equivalence = {}
    #print(str(len(string)))
    #print(str(len(strings_list)))
    if positions: # List of positions given by the user
        mutation_sites = positions
        print(positions)
        j = 0 #Initial pos without '.' and '-'
        for i, old_pos in enumerate(list(string)):
            pos_equivalence[j] = i
                #print('o' + str(i))
                #print('n' + str(j))
            j += 1
        strings_list = original_string
        #print(pos_equivalence)
    else:
        j = 0 #Initial pos without '.' and '-'
        for i, old_pos in enumerate(list(string)):
            if old_pos == '.' or old_pos == '-':
                continue
            else:
                pos_equivalence[j] = i
                #print('o' + str(i))
                #print('n' + str(j))
            j += 1

        # If mutation site includes '.' sample another mutation site. No repeat positions
        # List of possible positions without '.'
        possiblePos = [i for i,nt in enumerate(strings_list)] #list of indexes
        #print('possiblePos ' + str(len(possiblePos)))
        #print(possiblePos)
        if not isinstance(regions, list): 
            mutation_sites = list(np.random.choice(possiblePos, N, replace=False))
        else:
            by_region = regions.copy()
            #print(by_region)
            # Change probability of positions to be muted. If sum 1, only introduce mutations in that regions, otherwise introduce mutations in any region
            prob_pos = [list(range(r[1], r[2])) for r in by_region] 
            # Remove those positions that are not possible
            prob_pos_good = []
            for j in range(len(prob_pos)):
                if prob_pos[j][0] not in possiblePos:
                    #print('# Remove complete region')
                    possiblePos = possiblePos[:prob_pos[j-1][-1]]
                    #print('byregion ' + len(by_region))
                else:
                    #print('iterate')
                    prob_pos_good.append([]) # add the region
                    for i in prob_pos[j]: # Adjust last region
                        if i not in possiblePos:
                            #print('adjust')
                            #print(str(i))
                            continue
                            #prob_pos_good[j] = list(range(prob_pos[j][0], (i- 1)))
                        else:
                            prob_pos_good[j].append(i)
            #print(len(prob_pos_good))
            #print([len(reg) for reg in prob_pos])
            #print(len(possiblePos))
            #print(prob_pos[-1])
            #print(prob_pos_good[-1])
            #exit(-1)
            # If some region is not present, remove those possitions from possiblePos
    
            prob_values = list(np.repeat(0, len(possiblePos)))
            #prob_values = list(np.repeat(0, len(flattened(prob_pos))))
            for j in range(len(prob_pos_good)):
                for i in prob_pos_good[j][:-1]:
                    #print('probvalues')
                    #print(str(i))
                    #print(len(prob_values))
                    #print(prob_values[i])
                    #print(str(len(possiblePos)))
                    #print(str(prob_pos_good[j]))
                    #print(prob_pos_good[j])
                        #print(str(len(strings_list)))
                    prob_values[i] = by_region[j][3]/(len(prob_pos_good[j][:-1]))  # Probabilidad de que sea par 1/2: 50%, cual es la probabilidad de cada uno de los valores pares: 1/6==0.50/3
            if not round(sum([float(r[3]) for r in by_region])) == 1:
                probdivide = (1 - sum(prob_values))/sum([True if v == 0 else False for v in prob_values]) # Calculate percentage for the rest of regions until sum 1
                prob_values = [probdivide if v == 0 else v for v in prob_values]
                # If we have removed a region, we modify the probabilities per region, so we have to add the probability of that region to the others
            if sum(prob_values) < 1:
                #print('diff')
                diff = 1 - float(sum(prob_values))
                #prob_values = [v + float(diff/sum([len(block)-1 for block in prob_pos_good])) if v in [p for r in prob_pos_good for p in r] else v for v in prob_values]
                prob_values[0] = prob_values[0] + diff
                #print(str(sum(prob_values)))
            #print('mutate N ' + str(N) + ' ' + str(sum([True if f>0 else False for f in prob_values])))    
            if N < sum([True if f>0 else False for f in prob_values]): # More mutations that possible positions
                mutation_sites = list(np.random.choice(possiblePos, N, p = prob_values, replace = False))
            else:
               # print('Warning: More mutation sites than possible options are requested. Please, select wider regions or disable strict mode to introduce point mutations in different regions although with less probability')
                #exit(-1)
                N = sum([True if f>0 else False for f in prob_values])
                mutation_sites = list(np.random.choice(possiblePos, N, p = prob_values, replace = False))
        #print(sorted(mutation_sites))
    try:
        os.mkdir('mutations')
    except OSError as e: #[Errno 17] File exists: 'output'
        if e.errno != 17:
            raise
    with open('mutations/{}.nt_mutations.tsv'.format(header), 'w') as nt_out:
        for pos in sorted(mutation_sites): # Avoid cases in which the mutated base is identical to the original
            original_pos = pos_equivalence[pos] # see what is the position in the sequence with '.' and '-'
            random_nt = random.sample(list(filter(lambda x: x != ''.join(strings_list)[pos], bases)), 1)[0]
            nt_out.write('{}\t{}\t{}\t{}\t{}\n'.format(header, original_pos, pos, original_string[original_pos], random_nt))
            original_string[original_pos] = random_nt
        #print(original_string[original_pos])
    #tring_mutate = ''.join(strings_list.tolist())
    return(''.join(original_string.tolist()))


def load_positions(positions_file):
    with open(positions_file, 'r') as positions_align:
        positions=positions_align.read().rstrip()
        pos = np.where(np.array(list(positions)) == '1')[0].tolist()
    return(pos)

########################### TAXONOMY
def rank2number(rank): 
    """ Convert taxonomic ranks into numerical factors 0 based"""
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
        write_logfile('error', 'LIBS: convertTaxlevel' , '{} is not a valid taxlevel. Valid taxlevel are: kingdom, phylum, order, class, family, genus'.format(rank))
        exit(-3)


def estimate_mutations(rank, length = 1500): #1500 bp approximately 16S rRNA
    if rank == 1: #phylum
        cutoff =  1-0.75
        #write_logfile('info', 'FAKE TAXA', '{} cutoff {}'.format(rank, cutoff))
        return(cutoff * length)
    elif rank == 2: #class
        cutoff =  1-0.785
        #write_logfile('info', 'FAKE TAXA', '{} cutoff {}'.format(rank, cutoff))
        return(cutoff * length)
    elif rank == 3: #order
        cutoff =  1-0.82
        #write_logfile('info', 'FAKE TAXA', '{} cutoff {}'.format(rank, cutoff))
        return(cutoff * length)
    elif rank == 4: #family
        cutoff =  1-0.865
        #write_logfile('info', 'FAKE TAXA', '{} cutoff {}'.format(rank, cutoff))
        return(cutoff * length)
    elif rank == 5: #genus
        cutoff =  1-0.945
        #write_logfile('info', 'FAKE TAXA', '{} cutoff {}'.format(rank, cutoff))
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


def read_mutation_regions(f):
    """ Parse mutation regions and probabilities file: nameRegion\tstart\end\probability """ # Paper: PMC2562909: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2562909/
    #strict = f.readline().rstrip('\n').lstrip('#') # First line indicating if the mutations only appear in the defined regions. If the user wants to just rise the probability of mutations in some regions select Strict/NoStrict
    with open(f, 'r') as ff:
        l = [(line.rstrip().split('\t')[0], int(line.rstrip().split('\t')[1]), int(line.rstrip().split('\t')[2]), float(line.rstrip().split('\t')[3])) for line in ff]
    #l.insert(0, strict)
        return(l) # this will be input of by_region
    
def load_mutations_regions(regions):
    """ Load the default mutation regions just in case the user wants to use them """
    regions_default = {'V1': (69, 99), 'V2': (137, 242), 'V3': (433, 497), 'V4': (576, 682) , 'V5': (822, 879), 'V6': (986, 1043), 'V7': (1117, 1173), 'V8': (1243, 1294), 'V9': (1435, 1465)}
    # Pick the regions that the user wants and distribute probabilities between them
    l = []
    number_regions = len(regions)
    for r in regions:
        l.append((r, regions_default[r][0], regions_default[r][1], 1.0/number_regions)) # the sum of all probabilities must be 1
    return(l)




    
    
