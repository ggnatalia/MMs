#!/usr/bin/env python3
import argparse 
import os
import shutil
import pandas as pd
import matplotlib.pyplot as plt
import subprocess

def parse_arguments():
    """Parse the command line arguments and return an object containing them"""
    # Required
    general = argparse.ArgumentParser(description = 'Process input files')
    general.add_argument( '-m', '--mockName', type = str, required = True, help = 'Mock name. Ex: mock1')
    general.add_argument( '-o', '--output', type = str, required = True, help = 'Output directory. Preferably, name of the environment. Ex: Aquatic')
    general.add_argument( '--align', type = str, help = 'The user can provide an align file without creating it from scratch.')
    general.add_argument('--start', default = 0, type=int, help = 'Start position')
    general.add_argument('--end', default = 50000, type = int, help = 'End position')
    args = general.parse_args()
    return(args)


def runCommand(command, stdout = None, stderr = None):
    """Run the command and check the success of the subprocess. Return exit if it went wrong"""    
    exitcode = subprocess.call(map(str, command), stdout = stdout, stderr = stderr)
    if exitcode != 0:
        print('There must be some problem with "{}".\nIt\'s better to stop and check it'.format(' '.join(map(str, command))))
        exit(-1)

def load_table(filepath, rows = None , cols = None, path = '.', sep = ','):
    """ Load a csv file as pandas dataframe """
    if rows  and cols:
        data = pd.read_csv(filepath, sep = sep, index_col = 0)
    else:
        data = pd.read_csv(filepath, sep = sep, header=None)
    return(data)        

        
def main(args): 
    path = os.getcwd()
    projectPrefix = args.output
    projectPath = os.getcwd() + '/' + args.output
    mockPrefix = '{}.{}'.format(args.output, args.mockName)
    mockPath = projectPath + '/' + args.mockName
    os.chdir(mockPath)
    # List of sequences. 
    cutoffs = [0.00, 1-0.97, 1-0.945, 1-0.865, 1-0.82, 1-0.785, 1-0.75]
    labels = ['ASVs', 'species', 'genus', 'family', 'order', 'class', 'phylum']
    # Make distances with mothur
    command = ['mothur', '#pcr.seqs(fasta = {}/{})'.format(projectPath, args.align)]
    runCommand(command)
    command = ['mothur','#dist.seqs(fasta = {}/{}.pcr.align)'.format(projectPath, args.align.replace('.align', ''))]
    runCommand(command)
    # Create an abundance table like mothur requires
    df = load_table('{}/checkDB/{}.raw.abundances_original.tsv'.format(mockPath, mockPrefix), rows = True , cols = True, path = '.', sep = '\t')
    count_table = df.T
    count_table['total'] = count_table[list(count_table.columns)].sum(axis=1)
    count_table = count_table[['total']  + [col for col in count_table if col != 'total']]
    seqs2remove = set(count_table[count_table.total == 0].index)
    count_table = count_table[count_table.total != 0]
    count_table.index = count_table.index.rename('Representative_Sequence')
    count_table.to_csv ('{}/{}.count_table'.format(mockPath, mockPrefix), index = True, header=True, sep = '\t') 
    # Remove total 0 sequences from distance file
    dist = load_table('{}/{}.dist'.format(projectPath, projectPrefix), rows = False , cols = False, path = '.', sep = ' ')
    distcheck = dist[~dist[0].isin(seqs2remove) & ~dist[1].isin(seqs2remove)] # pick seqs ~ df.isin(x)  (== not in)
    distcheck.to_csv ('{}/{}.dist'.format(mockPath, mockPrefix), index = False, header=False, sep = ' ')
    
    shannonIndex = {}
    for c, label in zip(cutoffs, labels):
        # Cluster sequences
        command = ['mothur','#cluster(column = {}/{}.dist, count = {}/{}.count_table, cutoff = {})'.format(mockPath, mockPrefix, mockPath, mockPrefix, c)] # Saline.dist
        runCommand(command)
        shutil.move('{}.opti_mcc.list'.format(mockPrefix), '{}.{}.opti_mcc.list'.format(mockPrefix, label))
        command = ['mothur', '#get.sabund(list = {}.{}.opti_mcc.list, count = {}.count_table)'.format(mockPrefix, label, mockPrefix)]
        runCommand(command)
        command = ['mothur', '#summary.single(sabund = {}.{}.opti_mcc.sabund, calc = shannon)'.format(mockPrefix, label)]
        runCommand(command)
        # Save ShannonIndex value
        Hdf = load_table('{}.{}.opti_mcc.summary'.format(mockPrefix, label), rows = True , cols = True, path = '.', sep = '\t')
        shannonIndex[label] = float(Hdf['shannon'])
    
    H = pd.DataFrame.from_dict(shannonIndex, orient='index', columns = ['H'])
    H.to_csv ('{}/{}.shannonIndex.rank.tsv'.format(mockPath, mockPrefix), index = True, header=True, sep = '\t') 
    axes = H.plot.bar(rot=0, legend = False)
    plt.title('{} Shannon diversity index by rank'.format(mockPrefix))
    plt.xlabel('rank')
    plt.ylabel('Shannon Index (H)')
    plt.savefig('{}.shannonIndex.rank.png'.format(mockPrefix))

################################################################################################################
    
if __name__ == '__main__':
    main(parse_arguments())
