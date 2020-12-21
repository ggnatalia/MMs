#!/usr/bin/env python3

from libs.libs import *#Library with general functions
from libs.maths import *
from libs.plots import *
from libs.seqsClass import Sequence
from libs.mockClass import Mock

import pandas as pd
import numpy as np
from collections import defaultdict, Counter
import random
import os
from itertools import combinations, combinations_with_replacement
import itertools
from multiprocessing import Pool

class Enviro():
    
    # Class variables
    multiprocessing_globals = tuple()
    multiprocessing_globals_seqs = tuple()
    multiprocessing_globals_combinations = tuple()
    
    def __init__(self, prefix, enviro, Seqs, nASVs):
        self.prefix = prefix 
        self.enviro = enviro 
        self.Seqs = Seqs 
        self.nASVs = nASVs
    
    @classmethod
    def init_from_enviro(cls, nASVs, prefix, enviro, refEnviro = '/home/natalia/Projects/natalia/opt/makemocks/utils/speciesperEnviro.check.collapse.tsv', refTax =  '/home/natalia/Projects/natalia/DB/silva.nr_v138/silva.nr_v138.tax', ref = '/home/natalia/Projects/natalia/DB/silva.nr_v138/silva.nr_v138.align', nTaxa = 10000, rank = 5, figsize = (20, 20)):
        """ Return a Environment object from a given environment """
        print('Environment: \'{}\''.format(enviro))
        taxa = cls.subset_taxa_from_environment( enviro, refEnviro = refEnviro, nTaxa = nTaxa ) 
        # Collapse abundances of equal taxa
        taxAbun =  dict(Counter(taxa)) # Sequences from most abundant taxa will have been selected more times
        # tuple([(i,j) for i, j in taxAbun.items()])
        headers, neededSeqs = cls.subsetSilvaproportions( taxAbun, refTax = refTax, rank = rank )
        #print('headers ' + str(len(headers)))
        Seqs = cls.set_sequences( fastaFile = ref, refTax = refTax, degap = False, cleanHeader = True, splitChar = '\t', conservative = False, selected = list(headers)) # Original Seqs from Silva
        FakeSeqs = cls.make_fake_taxa(neededSeqs, rank, ref, refTax)
        TotalSeqs = Seqs|FakeSeqs
        #print('Seqs ' + str(len(Seqs)))
        #print('Fake ' + str(len(FakeSeqs)))
        #print('TotalSeqs' + str(len(TotalSeqs)))
        cls.plotTaxonomy2(TotalSeqs, '{}.silvasubset'.format(prefix), rank, figsize)
        return( Enviro(prefix, enviro, TotalSeqs, nASVs) ) 
    
    @classmethod
    def make_fake_taxa(cls, neededSeqs, rank, ref, refTax):
        """ Return a Seqs set with the fake taxa that fall short """
        fakeSeqs = set()
        seqs2fake = {}
        # Add sequences from each taxa that fail because there were not enough sequences
        for abun, seqsHeaders in neededSeqs: # With one seq, it can be generated so many fake species as needed, but for the shake of variety, it's better to distribute the number of species from which the fake will be generated
            if abun <= len(seqsHeaders):
                for s in random.sample(seqsHeaders, abun):
                    seqs2fake[s]= 1 
            else:
                basic = int(abun/len(seqsHeaders))
                for s in random.sample(seqsHeaders, len(seqsHeaders)):
                    seqs2fake[s] = basic 
                privilegeSeq = random.sample(list(seqs2fake.keys()), 1)[0] # Add the rest of the division to one random sequence.
                seqs2fake[privilegeSeq] = seqs2fake[privilegeSeq] + int(abun%len(seqsHeaders)) # Add the rest of the division to one random sequence.
        moreSeqs = cls.set_sequences( fastaFile = ref, refTax = refTax, degap = False, cleanHeader = True, splitChar = '\t', conservative = False, selected = list(seqs2fake.keys())) # Set pf sequences to make more sequences from them
        for s in moreSeqs:
            Nposmax = int(estimate_mutations(rank, length = len(s.seq.replace('.','').replace('-',''))))
            #print('Nposmax ' , str(Nposmax))
            Nstrains = int(seqs2fake[s.header]) 
            #print('mutant seqs I want ' + s.header + ' ' + str(Nstrains))
            newS = s.generatemutantASVs(Nstrains = Nstrains, Nposmax = Nposmax, start = 0, end = 50000, include_original = False)
            #print('newS ' + str(len(newS)))
            fakeSeqs.update(newS)                
        return(fakeSeqs)
             
    
    
    
    @classmethod
    def init_from_taxa(cls, nASVs, prefix, rank, taxa, taxa_abundances = [], refTax =  '/home/natalia/Projects/natalia/DB/silva.nr_v138/silva.nr_v138.tax', ref = '/home/natalia/Projects/natalia/DB/silva.nr_v138/silva.nr_v138.align', figsize = (20, 20)):
        """ Create an environment from a list of user taxa or from random taxa at the selected rank """
        enviro = 'taxas'
        write_logfile('info', 'SUBSET SEQUENCES', 'You have pass a list of taxa {} from which take sequences and make ASVs'.format(taxa))
        if len(taxa_abundances) > 0 :
            taxa_abundances = [i*100 for i in  taxa_abundances]  # To subset 10000 sequences from SILVA
        elif len(taxa_abundances) == 0 :
            # Recreate the list taxa_abundances but with random values
            write_logfile('critical', 'init_from_taxa', 'NOT WORKING when using ESTMATE ABUNDANCES FUNCTION!!!')
            taxa_abundances = estimate_abundances(len(taxa), 100) # Do % of each taxa
        #else: # select random taxa instead of sequences from DB
        #    taxa = Enviro.subset_random_taxa(args.refTax, taxlevel, nTaxa = 0)
        #           nTaxa = len(taxa)
        #    taxaAbund = estimate_abundances(nTaxa, args.nASVs)
        taxAbun = dict(zip(taxa, taxa_abundances))
        #taxAbun = tuple(zip(taxa, taxa_abundances))
        headers, neededSeqs = cls.subsetSilvaproportions( taxAbun, refTax = refTax, rank = rank )
        Seqs = cls.set_sequences( fastaFile = ref, refTax = refTax, degap = False, cleanHeader = True, splitChar = '\t', conservative = False, selected = list(headers) ) 
        FakeSeqs = cls.make_fake_taxa(neededSeqs, rank, ref, refTax)
        TotalSeqs = Seqs|FakeSeqs
        print('TotalSeqs' + str(len(TotalSeqs)))
        ### AHORA MISMO INUTIL
        cls.plotTaxonomy2(TotalSeqs, '{}.silvasubset'.format(prefix), rank, figsize)
        return( Enviro(prefix, enviro, TotalSeqs, nASVs) ) 
    
    @classmethod   
    def init_from_seqs(cls, prefix, rank, seqs, nASVs, minseqs = 1, refTax =  '/home/natalia/Projects/natalia/DB/silva.nr_v138/silva.nr_v138.tax', ref = '/home/natalia/Projects/natalia/DB/silva.nr_v138/silva.nr_v138.align'):
        enviro = 'seqs'
        write_logfile('info', 'SEQUENCES PROVIDED', 'Sequences {} are treated as OTUs, mutant ASVs will be generated until complete {} ASVs'.format(seqs, nASVs))
        if seqs:
            if isinstance(seqs, list):
                if nASVs < len(seqs) and seqs:
                    write_logfile('warning', 'SEQUENCES PROVIDED', 'You want less total ASVs ({}) than provided sequences ({})'.format(nASVs, len(seqs)))
                    print('You want less total ASVs ({}) than provided sequences ({})'.format(nASVs, len(seqs)))
                    print('Reduce your sequences, rise the number of ASVs and if you do not want strains, fix ASVmean = 0')
                    Seqs = cls.set_sequences(fastaFile = ref, refTax = refTax, degap = False, cleanHeader = True, splitChar = '\t', conservative = False, selected = seqs) # Subset the list of sequences from SILVA DB
                else:
                    write_logfile('info', 'RANDOM SEQUENCES', 'You have not pass any align, environment, list of sequences or taxas, random sequences will be subsetted to make ASVs')
                    Seqs = cls.set_sequences(fastaFile = ref, refTax = refTax, degap = False, cleanHeader = True, splitChar = '\t', conservative = False, Nrandom = random.randint(minseqs, nASVs)) # Subset the list of sequences from SILVA DB
            else: #seqs is a fasta file with the sequences:
                seqs = list(fasta2dict(seqs).keys()) # Extract header of the sequences. I do this to take the alignment sequence from silva
                Seqs = cls.set_sequences(fastaFile = ref, refTax = refTax, degap = False, cleanHeader = True, splitChar = '\t', conservative = False, selected = seqs)
        return( Enviro(prefix, enviro, Seqs, nASVs) ) 
    
    @classmethod
    def init_from_file(cls, prefix, path, inputfile, refTax = '/home/natalia/Projects/natalia/DB/silva.nr_v138/silva.nr_v138.tax'):
        """ Open and load an align file """
        if os.path.isfile('{}/{}'.format(path, inputfile)): 
            write_logfile('info', 'PREVIOUS ALIGN', 'The file {}/{} will be used as reference'.format(path, inputfile))
            enviro = 'inputfile'
            Seqs = cls.set_sequences(fastaFile = '{}/{}'.format(path, inputfile), refTax = refTax, degap = False, cleanHeader = False, splitChar = '\t', conservative = False)
        else:
            write_logfile('warning', 'INPUT ALIGN', '{} not found in {}.'.format(inputfile, projectPath))
            exit(-2)
        return( Enviro(prefix, enviro, Seqs, len(Seqs)) )
    
    
    
    
    @classmethod 
    def distances(cls, Seqs, prefix, region, start, end, cutoff, cpus, complete = False): # igual se puede combinar con makeASVs
        """ Calculate distances in an specific region (plot a heatmap) and filter sequences based on them. """
        cutoff = float(cutoff)
        cpus = int(cpus)
        # Trim sequences according to start end position. Return set of Seqs objects
        # https://www.geeksforgeeks.org/copy-python-deep-copy-shallow-copy
        cls.multiprocessing_globals = {s.trimregion(start, end) for s in Seqs}
        equivalence = np.array([s.header for s in cls.multiprocessing_globals]) #"https://docs.python.org/2/library/stdtypes.html#mapping-types-dict: If keys, values and items views are iterated over with no intervening modifications to the dictionary, the order of items will directly correspond"
        if cpus == 1 or len(Seqs) < 100:
            write_logfile('info', 'CONVERT SEQS', 'Start 1 cpu')
            cls.multiprocessing_globals_seqs = tuple(map(Sequence.nt2dict, cls.multiprocessing_globals))
        else:
            write_logfile('info', 'CONVERT SEQS', 'Start multiprocessing')
            with Pool(cpus) as pool:
                cls.multiprocessing_globals_seqs = tuple(pool.map(Sequence.nt2dict, cls.multiprocessing_globals))
        write_logfile('info', 'CONVERT SEQS', 'Finish')
        cls.multiprocessing_globals = tuple()
        cls.multiprocessing_globals_combinations = tuple(combinations_with_replacement(list(range(len(cls.multiprocessing_globals_seqs))),2))
        if cpus == 1 or len(Seqs) < 100:
            write_logfile('info', 'DISTANCE CALCULATIONS', 'Start 1 cpu')
            distances = list(map(cls.calc_dist, cls.multiprocessing_globals_combinations))
        else:
            write_logfile('info', 'DISTANCE CALCULATIONS', 'Start multiprocessing')
            with Pool(cpus) as pool:
                distances = list(pool.map(cls.calc_dist, cls.multiprocessing_globals_combinations))
        write_logfile('info', 'DISTANCE CALCULATIONS', 'Finish')
        
        length2split = list(reversed(range(1, len(equivalence) + 1)))   # split the list into one row per original sequence: for 3 sequences: [(1,1),(1,2),(1,3)],[(2,2), (2,3)],[(3,3)]: split 3, 2, 1
        iterable_distances = iter(distances) 
        dist_trian = [list(itertools.islice(iterable_distances, elem)) for elem in length2split]
        
        fill = list(range(0, len(equivalence)))
        # Prepare 0 to complete the df as a matrix
        write_logfile('info', 'DISTANCE CALCULATIONS', 'Finish')
        arr = np.array(list(map(lambda x, y: x*[0] + y, fill, dist_trian)))
        # Make a df with distance calcs for all the sequences
        df = pd.DataFrame(arr, columns = equivalence, index = equivalence)
        # Logical df with False when value ==-1, to remove this columns
        mask = df.values!=-1
        # Remove these columns and the corresponding rows
        df = df.drop(index = list(np.where(mask==False)[1]), columns = list(np.where(mask==False)[1]))
        seqstoremove = list(np.where(mask==False)[1])
        SeqsFilteredheaders = list(np.delete(equivalence, seqstoremove))
        SeqsFiltered = {s for s in Seqs if s.header in SeqsFilteredheaders}
        if complete:
            write_table(df, title = '{}.distances'.format(prefix), rows = list(df.index) , columns = list(df.columns), dataframe = None)
            return(SeqsFiltered, df)
        else:
            return(SeqsFiltered)
        
    @classmethod
    def clear_multiprocessing_globals(cls): # Clean class variable instead of overwrite it
        cls.multiprocessing_globals_seqs = tuple()
        cls.multiprocessing_globals_combinations = tuple()
        cls.multiprocessing_globals = tuple()
    
    @classmethod
    def calc_dist(cls, args): 
        """ args is a tuple of two index combinations with the index of the sequences between calculate distances """
        d = calculate_distance(cls.multiprocessing_globals_seqs[args[0]], cls.multiprocessing_globals_seqs[args[1]])
        if d != 0:
            return(d)
        elif args[0] == args[1]: # Si es el mismo indice contra si mismo, d = 0 they are identical
            return(0.0)
        else:    
            return(-1.0) # Return -1 to filter by this value after and remove those sequences which are identical

    # Methods for Enviro objects
    def makeASVs(self, region, start, end, ASVsmean, cutoff , cpus, figsize = (20,20)): # ADD ASV/OTU calculation
        """ Create fake ASVs from sequences """
        # Subset Seqs, do ASVs, until complete nASVs
        ASVs = set()
        sampleSeqs = set()
        while len(ASVs) < self.nASVs:  # Repeat until complete number of ASVs
            newS = random.sample(self.Seqs, 1)[0] # Subset one sequence. Return a list with one element
            if newS.header not in [s.header for s in sampleSeqs]: # if that sequence has not been yet taken 
                sampleSeqs.add(newS)
                newSasvs = newS.generatemutantASVs(Nstrains = None, Nmean = ASVsmean, Nposmax = 45, start = start, end = end)
                #Check distances
                ASVsdiff = self.distances(Seqs = newSasvs, prefix = self.prefix, region = region, start = start, end = end, cutoff = cutoff, cpus = cpus)
                #write_logfile('warning', 'SUBSET RANDOM SEQS', 'len(ASVs) {}\tlen(ASVsdiff) {}'.format(len(newSavs), len(ASVsdiff)))
                i=0
                while len(ASVsdiff) < len(newSasvs): # If one sequence have been removed 'cause it is identical to another one, repeat to take ASVs from that sequence, otherwise some taxa can be minusvalorated
                    newSavs = newS.generatemutantASVs(Nstrains = None, Nmean = ASVsmean, Nposmax = 45, start = start, end = end)
                    #Check distances
                    ASVsdiff = self.distances(Seqs = newSasvs, prefix = self.prefix, region = region, start = start, end = end, cutoff = cutoff, cpus = cpus)
                    i = i+1
                    if i > 10:
                        #write_logfile('warning', 'SUBSET RANDOM SEQS2', 'len(ASVsdiff) {}\tlen(newSavs) {}'.format(len(newSavs), len(ASVsdiff)))
                        break
                ASVs.update(ASVsdiff)
            else:
                continue
        ASVsdiff, distdf = self.distances(Seqs = ASVs, prefix = self.prefix, region = region, start = start, end = end, cutoff = cutoff, cpus = cpus, complete = True) 
        #ASVsdiffselected = set(random.sample({s for s in ASVs if s.header not in [so.header for so in SeqsSilvaselected]}, args.nASVs-len(SeqsSilvaselected))) # No matter if original sequences are included or not
        ASVsdiffselected = set(random.sample(ASVsdiff, self.nASVs))
        #self.Seqs = SeqsSilvaselected.union(ASVsdiffselected) # No matter if original sequences are included or not
        self.Seqs = ASVsdiffselected # I REALLY WANT TO MODIFY IT, I want to remove extra Sequences
        write_logfile('info', 'ENVIRONMENT', 'Plotting distances heatmap')
        self.plot_distances(df = distdf, region = region, cutoff = cutoff, text = None, symmetric = True, figsize = (20, 20))
        return(self)# I really want to modify it
    
    
    #### PLOTS & OUTPUTS
    def plot_distances(self, df, region = '16S', cutoff = 0.03, text = None, symmetric = None, figsize = (20, 20)):
        """ Plot distances among selected sequences """
        # Update df with the final sequences that finally have been included:
        shortdf = df.loc[ [idseq.header for idseq in self.Seqs] , [idseq.header for idseq in self.Seqs] ] # select columns and rows
        plot_heatmap(shortdf, outputDir = os.getcwd(), title = '{}.{}.distances'.format(self.prefix, region), vmin=0, vmax=1, center= 0, legendtitle = 'distance', text = text, symmetric = symmetric, figsize = figsize)
    
    def plot_taxonomy(self, rank):
        """ Barplot stacked and percent of the taxonomy of the sequences included in the mock """
        taxa = [s.tax.split(';')[rank] for s in self.Seqs]
        taxAbun =  dict(Counter(taxa))
        # Create a df of taxonomy and counts
        df = pd.DataFrame.from_dict(taxAbun, orient='index', columns=['counts']).transpose()
        barplotpc(df, outputDir = os.getcwd(), title = '{}.pctax'.format(self.prefix), figsize = (20, 20), ylab = 'taxon pc %', xlab = 'Mock', textSize = 8 )
        write_table(df, title = '{}.pctax'.format(self.prefix), outputDir = '.', rows = list(df.index) , columns = list(df.columns), dataframe = None)

    def write_output(self): 
        """ Write an align file, a fasta file and a taxonomy file of the sequences that will take part in the mocks """
        with open('{}.align'.format(self.prefix), 'w') as align, open('{}.fasta'.format(self.prefix), 'w') as fasta, open('{}.taxonomy'.format(self.prefix), 'w') as taxonomy:
                Seqscopy = [s.copy() for s in self.Seqs]
                for s in Seqscopy:
                    align.write('>{}\n{}\n'.format(s.header, s.seq))
                    fasta.write('>{}\n{}\n'.format(s.header, s.deGap().seq))
                    taxonomy.write('{}\t{}\n'.format(s.header, s.tax))
                return(True)

    # DE MOMENTO CANCELADA
    def plot_taxon_entropy(self, ref = '/home/natalia/Projects/natalia/DB/silva.nr_v138/silva.nr_v138.align', refTax =  '/home/natalia/Projects/natalia/DB/silva.nr_v138/silva.nr_v138.tax', rank = 3, title = None):
        """ Calculate and plot the entropy of the all the taxons """
        # Need to look for all the sequences that have the same taxon
        taxons = set() #taxons = {'Micrococcales', 'Alteromodalaes'}
        for s in self.Seqs:
            taxons.add(s.tax.split(';')[rank])
        print(taxons)
        # New strategy: read everything in silva tax, extract sequences from the corresponding taxa and then select sequences of this taxa    
        silvaTaxas = loadTaxa(refTax  = refTax, rank = rank)    
        taxaselectedSequences = { taxa : Sequence.set_sequences(ref, refTax = refTax, rank = rank, selected = [k for k,v in silvaTaxas.items() if v == taxa ]) for taxa in taxons}    
        for taxa in taxaselectedSequences.keys():
            print(taxa)
            length = len(set([len(s.seq)for s in taxaselectedSequences[taxa]]))
            # All the squences in silva have length: 50000
            if length == 1:
                #Now create the array with just the split sequences:
                seqsArray = align_to_array(taxaselectedSequences[taxa], expand = False) # dtype='<U1': string 1 length
                title = '{}.{}'.format(self.prefix, taxa)
                plot_entropy(seqsArray, title, outputDir = self.prefix)

    @staticmethod
    def set_sequences(fastaFile, refTax, degap = False, cleanHeader = True, splitChar = '\t', conservative = False, rank = None, selected = [], Nrandom = 0): # TESTED! ################## MIRAR CASO FASTA MULTILINEA
        """ Make a list of sequence objects (or select some based on the header) from fasta"""
        with open(fastaFile, 'r') as fasta:
            all_sequences_list = fasta.readlines()
            silva_taxa = loadTaxa(refTax, rank = rank)
            Seqs = set()
            if Nrandom != 0:
                selected = random.sample(silva_taxa.keys(), Nrandom)
            for i in range(0, len(all_sequences_list), 2): 
                if cleanHeader:
                    header = simplifyString(all_sequences_list[i].lstrip('>').rstrip('\n'), splitChar, conservative)
                else:
                    header = all_sequences_list[i].lstrip('>').rstrip('\n')
                if not 'N' in all_sequences_list[i + 1]: #checkpoint: some sequences from silva contain N, exclude those sequences
                    if not selected: # select all the sequences from the fasta file
                        #write_logfile('debug', 'Sequence.set_Sequences', selected)
                        Seq = Sequence(header, all_sequences_list[i + 1].rstrip('\n'), Sequence.assign_taxonomy(header, silva_taxa, rank), 0)
                        #write_logfile('debug', 'Sequence.set_Sequences NOT selected ', Seq.header)
                    else: # Select just the sequences of the list
                        #write_logfile('debug', 'Sequence.set_Sequences', selected)
                        if header in selected: #add sequence to the set
                            Seq = Sequence(header, all_sequences_list[i + 1].rstrip('\n'), Sequence.assign_taxonomy(header, silva_taxa, rank), 0)
                            #write_logfile('debug', 'Sequence.set_Sequences SELECTED', Seq.header)
                        else:
                            continue
                    if degap:
                        #write_logfile('debug', 'Sequence.set_Sequences degap', degap)
                        Seq = Seq.deGap()
                    write_logfile('debug', 'Sequence.set_Sequences ADD', Seq.header)
                    Seqs.add(Seq)  
                else:
                    continue
            write_logfile('debug', 'Sequence.set_sequences', 'If there are repeated sequences in the fasta, there will be repeated sequences in the set Sequences, each Sequence object is different although it contains the same elements')
            return(Seqs)
    
    @staticmethod
    def subset_taxa_from_environment(enviro, refEnviro = '/home/natalia/Projects/natalia/opt/makemocks/utils/speciesperEnviro.check.collapse.tsv',  nTaxa = 0): #the table consists of genus
        """ From a particular environment, select different taxas that are supposed to be more abundant in them """
        nTaxa = int(nTaxa)
        # Open and parse enviro reference
        enviro_df = loadEnviro(refEnviro)
        if enviro == 'random': # Subset an environment randomly
            enviro = ''.join(random.sample(list(enviro_df.columns), 1) )
        else:
            enviro = enviro
        # Subset environment from enviro reference
        enviro_subset = enviro_df[enviro]
        if nTaxa == 0:
            nTaxa = random.randint(1, len(enviro_df.index))
        # Subset taxa from enviro. Each taxa has a different possibility of be subsetted.
        #taxas = list(enviro_df.index)
        # 0 times in the environments means: There is none in that environment or it is in low abundance, so it has not been detected.
        # We cannot distinguish it, so we give a low weigth to appear
        enviro_subset_no0 = enviro_subset.replace({0: float(0.005)})
        # Remove 'Unresolved' row, because there isn't in silva DB
        #enviro_subset_no0 = enviro_subset_no0.drop(['Unresolved']) # no more needed
        #enviro_subset_no0 = enviro_subset_no0.drop(['Inconsistent']) #no more needed
        weights = enviro_subset_no0.apply(lambda x: x/float(enviro_subset_no0.sum())) # % 1
        #print(weights)
        # pd.sample renormalize to 1
        # Subset directly from weigths pandas series # ADD random_state to fix a seed? 
        selectedTaxas = list(enviro_subset_no0.sample(n = nTaxa, replace = True, weights = weights).index)
        #selectedTaxas = list(enviro_subset_no0.sample(n = nTaxa, replace = True, weights = enviro_subset_no0).index)
        #print(selectedTaxas)
        return(selectedTaxas)

    @staticmethod
    def subsetSilvaproportions(taxAbun, refTax =  '/home/natalia/Projects/natalia/DB/silva.nr_v138/silva.nr_v138.tax', ref = '/home/natalia/Projects/natalia/DB/silva.nr_v138/silva.nr_v138.align', rank = 3):
        """ From a dict of taxas {'tax':Abun} take, sequences from from silva """
        silvaTaxas = loadTaxa(refTax  = refTax, rank = rank) # {seq:tax}
        silvaSeqs = list(fasta2dict(ref).keys())
        for seq in list(silvaTaxas.keys()):
            if not seq in silvaSeqs:
                silvaTaxas.pop(seq)
        taxonSeqs = reversedict(silvaTaxas) # {taxa:[seq,...]}
        # Some sequences from silva have 'N', we do not want to include them, so exclude these sequences
         # Only sequences without 'N'
        # Check good sequences in taxonSeqs, remove those that does not appear.
        
        subsetSeqs = set()
        moreSeqs = set() # set of species for which more sequences are required {(abun,(seq1,seq2)),(abun,(seq3,seq4,seq5))}
        for tx, abun in taxAbun.items(): # si taxAbun es una tupla: for tx, abun in taxAbun:
            if tx in list(taxonSeqs.keys()):
                if abun > len(taxonSeqs[tx]):
                    write_logfile('info', 'SUBSET RANDOM SEQS', '{} Ask for {}, but only {}. Fake representative sequences will be created to satisfy this number'.format(tx, abun, len(taxonSeqs[tx])))
                    moreSeqs.add((abun - len(taxonSeqs[tx]) , tuple(taxonSeqs[tx])))
                    abun = len(taxonSeqs[tx])
                # If there are not enough sequences, create more at  the given taxonomic level!!!!!!!!! MODIFY IN FUTURE
                seq = random.sample(taxonSeqs[tx], abun)
                #print(seq)
                # If the sequence has been added in a previous loop, it does not matter 'cause it is a set, no duplicates
                subsetSeqs = subsetSeqs.union(set(seq))
            else:
                write_logfile('error', 'SUBSET RANDOM SEQS', '{} not present in the DB {} sequences'.format(tx, abun))
        return(subsetSeqs, moreSeqs)

    @staticmethod
    def subset_random_taxa(refTax, rank, nTaxa = 0):
        """ Subset random taxa from reference Taxonomy DB """
        silvaTaxas = loadTaxa(refTax  = refTax, rank = rank) # {seq:tax}
        taxonSeqs = reversedict(silvaTaxas) #{taxa:[seq,...]}
        if nTaxa == 0:        
            nTaxa = random.randint(1, len(silvaTaxas))
        taxa = random.sample(silvaTaxas.keys(), nTaxa)
        return(taxa)

#### PROVISIONAL
    @staticmethod # Este es ahora para probar si el subset de silva respeta proporciones, luego se borra
    def plotTaxonomy2(Seqs, prefix, rank, figsize = (20, 20)):
        """ Bar plot stacked and percent of the taxonomy of the sequences included in the silva subset """
        taxas = [s.tax.split(';')[rank] for s in Seqs]
        count = {}
        for i in taxas: 
            count[i] = count.get(i, 0) + 1    
        # Create a df of taxonomy and counts
        df = pd.DataFrame.from_dict(count, orient='index', columns=['counts'])
        df = df.transpose()
        barplotpc(df, outputDir = os.getcwd(), title = '{}.pctax'.format(prefix), figsize = figsize, ylab = 'taxon pc %', xlab = 'Mock', textSize = 8 )   
        write_table(df, title = '{}.pctax'.format(prefix), outputDir = '.', rows = list(df.index) , columns = list(df.columns), dataframe = None)















