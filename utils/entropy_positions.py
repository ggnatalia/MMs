#!/usr/bin/env python3


# Basic
import argparse 
import os
import sys
import numpy as np
import pandas as pd
import shutil
from multiprocessing import Pool
import itertools
import scipy.stats as sp
import matplotlib.pyplot as plt

def parse_arguments():
    """Parse the command line arguments and return an object containing them"""
    # Required
    general = argparse.ArgumentParser(description = 'Process input files')
    general.add_argument( '-o', '--output', type = str, default = 'variable_positions.txt', required = True, help = 'Output file.')
    general.add_argument( '--raw_positions', nargs = '+', default = '', help = 'List of positions to be muted o file with those positions. 0 based. ')
    general.add_argument( '--length', default = 50000, type = int, help = 'Total length of the sequences, by default 50000 which is from SILVA alignment')
    general.add_argument( '--rank', type = str, default = 'family', help = 'Rank to subset taxa: phylum, order, class, family, genus')
    general.add_argument( '-tx', '--taxa', type = str, nargs = '+', default = [], help = 'List of taxa')
    general.add_argument( '-ref','--ref', type = str, default = '', help = 'SILVA alignment reference')
    general.add_argument( '-refTax','--refTax', type = str, default = '', help = 'SILVA alignment TAX reference')    
    general.add_argument('--cutoff', type = float, default = 0.5, help = 'Cutoff to select positions with and entropy higher than it')
    general.add_argument('--cpus', default = 12, type = int, help = 'Number of threads')   
    args = general.parse_args()
    return(args)






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
        print('{} is not a valid taxlevel. Valid taxlevel are: kingdom, phylum, order, class, family, genus'.format(rank))
        exit(-3)


def loadTaxa(refTax = '/home/natalia/Documentos/TRABAJO/MMs/DB/silva.nr_v138.tax', rank = None):
    """ For a seq return taxonomy. Prepare for parsing silvaDB taxonomy"""
    with open(refTax, 'r') as taxReference:
        all_sequences_list = taxReference.readlines()
        if not rank:
            taxas = {line.rstrip('\n').split('\t')[0] : line.rstrip('\n').split('\t')[1] for line in all_sequences_list}
        else:
            taxas = {line.rstrip('\n').split('\t')[0] : line.rstrip('\n').split('\t')[1].split(';')[rank] for line in all_sequences_list}
        return(taxas)


def validate_Sequence(combo):
    header = combo[0]
    seq = combo[1]
    header = header.lstrip('>').rstrip('\n').split('\t')[0] # assuming SILVA db
    Seq = ()
    if header in sequences: 
        if not 'N' in seq: #checkpoint: some sequences from silva contain N, exclude those sequences
        #add sequence to the set
            Seq = seq.rstrip('\n')
        return(Seq)

        
def set_sequences(fastaFile, cpus = 1):
    """ Make a list of sequence objects (or select some based on the header) from fasta """
    with Pool(cpus) as pool:
        with open(fastaFile, 'r') as fasta:
            Seqs = set(list(filter(None, pool.imap(validate_Sequence, itertools.zip_longest(*[fasta]*2)))))
    return(Seqs)


def extract_positions(taxa, ref = '/home/natalia/Documentos/TRABAJO/MMs/DB/silva.nr_v138.align', refTax = '/home/natalia/Documentos/TRABAJO/MMs/DB/silva.nr_v138.tax', rank = 'family', cutoff = 0.2, cpus = 6):
    """ From a particular taxa, extract more variable positions regarding a cutoff"""
    positions = []
    # First load the taxonomy for all the sequences in the reference
    silva_taxa = loadTaxa(rank = rank2number(rank))
    # Collect the name of all the sequences from a corresponding taxa
    global sequences
    sequences = set()
    for k, v in silva_taxa.items():
        if v in taxa:
            sequences.add(k)
    # Load those sequences to calculate entropy per position
    Seqs = set_sequences(ref, cpus = cpus )
    print(len(Seqs))
    # Convert sequences to array
    seqsArray = np.array(list(tuple(s) for s in Seqs))
    #print(seqsArray)
    # Plot entropies
    nrow, ncol = seqsArray.shape
    # SIMPLY CODE:
    entropies = {}
    for pos in range(ncol):
        entropies[pos] = sp.entropy(list(map(float,{nt:i for nt, i in zip(*np.unique(seqsArray[:,pos], return_counts=True)) if nt != '.'}.values())))
    # entropies = {pos: sp.entropy(list(map(float,{nt:i for nt, i in zip(*np.unique(array[:,pos], return_counts=True)) if nt != '.'}.values()))) for pos in range(ncol)}
    # Impossible to see 50000 positions, need to remove some
    higherentropies = { pos:S for pos, S in entropies.items() if S > cutoff }
    #plotData = sorted(tuple(higherentropies.items()), key = lambda x: x[0])
    #plotDataDF = pd.DataFrame(plotData, columns = ['Pos', 'Entropy']) #  0 based
    #plotDataDF = plotDataDF.set_index('Pos')
    #barplot(plotDataDF, outputDir = outputDir, title = '{}.entropy'.format('_'.join(taxa)), figsize = (12,12), T = True, ylab = 'Entropy', xlab = 'Position', textSize = 8 )
    y_pos = np.arange(len(higherentropies.keys()))
    # Create bars
    plt.bar(y_pos, higherentropies.values())
    # Create names on the x-axis
    plt.xticks(y_pos, higherentropies.keys())
    plt.xticks(rotation = 90)
    plt.ylabel('Entropy')
    plt.xlabel('Positions in SILVA alignment (0 based)')
    plt.savefig('Entropy_{}.png'.format('_'.join(taxa)))
    # Show graphic
    #plt.show()
    return(higherentropies)



def main(args): 
    
    # Create final string
    finalSeq = np.zeros(args.length, dtype=int)
    if args.raw_positions: # read a file with positions or a list of positions
        if os.path.isfile(args.raw_positions[0]) and len(args.raw_positions) == 1:
            with open(args.raw_positions[0], 'r') as rawpositionsfile:
                positions = rawpositionsfile.read().rstrip('\n').split('\n')
        elif isinstance(args.raw_positions, list):
            positions = args.raw_positions
    else: # Look for entropy higher than cutoff
        if not args.ref or not args.refTax:
            makemocks_home = os.path.abspath(os.path.dirname(os.path.realpath(__file__))) + '/../' # this script is within utils
            DB = makemocks_home + '/DB/'
            ref  = DB + '/' + 'silva.nr_v138.align'
            refTax = DB + '/' + 'silva.nr_v138.tax'
        else:    
            ref  =  args.ref
            refTax = args.refTax
    
    # Check if DB are present
        if not os.path.isfile(refTax) and not os.path.isfile(ref):
            print('DB {} has not been found. Please be sure to run /path/to/MMs/bin/make_databases.py to configure databases'.format(refTax.replace('.tax', '')))
            exit(-1)
    
        higherentropies = extract_positions(taxa = args.taxa, ref = ref, refTax = refTax, rank = args.rank, cutoff = args.cutoff, cpus = args.cpus)
        positions = list(higherentropies.keys())
    # Change positions
    print(positions)
    finalSeq[list(map(int,positions))] = 1        
    # OUTPUT
    with open('{}.txt'.format(args.output), 'w') as output:
        # Write the numpy array as an string of 0=False, 1 = True
        output.write(''.join(map(str,list(finalSeq))))



################################################################################################################
    
if __name__ == '__main__':
    main(parse_arguments())



