#!/usr/bin/env python

from libs.libs import *
from libs.maths import *
from libs.plots import *
from libs.seqsClass import Sequence
from libs.mockClass import Mock
from libs.distances import calculate_distances_cython

import pandas as pd
import numpy as np
from collections import defaultdict, Counter
import random
import os
from itertools import combinations, combinations_with_replacement
import itertools
from multiprocessing import Pool
import cython


class Enviro():
    
    # Class variables
    multiprocessing_globals = tuple()
    multiprocessing_globals_seqs = tuple()
    multiprocessing_globals_combinations = tuple()
    
    silva_taxa = dict()
    rank = None
    degap = False
    notax = False # To avoid taxonomic annotation
    selected = []
    
    def __init__(self, prefix, enviro, Seqs, nASVs):
        self.prefix = prefix 
        self.enviro = enviro 
        self.Seqs = Seqs 
        self.nASVs = nASVs
    
    @classmethod
    def init_from_enviro(cls, nASVs, prefix, enviro, refEnviro = '/home/natalia/Projects/natalia/opt/makemocks/utils/speciesperEnviro.check.collapse.tsv', refTax =  '/home/natalia/Projects/natalia/DB/silva.nr_v138/silva.nr_v138.tax', ref = '/home/natalia/Projects/natalia/DB/silva.nr_v138/silva.nr_v138.align', nTaxa = 10000, rank = 5, cpus = 20, by_region = None):
        """ Return a Environment object from a given environment """
        print('Environment: \'{}\''.format(enviro))
        taxa = cls.subset_taxa_from_environment( enviro, refEnviro = refEnviro, nTaxa = nTaxa ) 
        #print(len(taxa))
        # Collapse abundances of equal taxa
        taxAbun =  dict(Counter(taxa)) # Sequences from most abundant taxa will have been selected more times
        #print('subsetSilvaproportions')
        write_logfile('info', 'Enviro', 'Subset sequences')
        headers, neededSeqs = cls.subsetSilvaproportions( taxAbun, refTax = refTax, ref = ref, rank = rank, cpus = cpus )
        #print('headers ' + str(len(headers)))
        #print('neededSeqs ' + str(len(neededSeqs)))
        write_logfile('info', 'Enviro', 'Select sequences')
        Seqs = cls.set_sequences( fastaFile = ref, refTax = refTax , cpus = cpus, rank = rank, selected = list(headers), degap = False) # Original Seqs from Silva
        #print('fakeASVs')
        write_logfile('info', 'Enviro', 'Complete sequences')
        if neededSeqs:
            FakeSeqs = cls.make_fake_taxa(neededSeqs, rank = rank, cpus = cpus, ref = ref, refTax = refTax, by_region = by_region, by_pos = None)
            TotalSeqs = Seqs|FakeSeqs
        else:
            TotalSeqs = Seqs
        return( Enviro(prefix, enviro, TotalSeqs, nASVs) ) 
    
    @classmethod
    def make_fake_taxa(cls, neededSeqs, rank, cpus, ref, refTax, by_region = None, by_pos = None):
        """ Return a Seqs set with the fake taxa that fall short """

        fakeSeqs = set()
        seqs2fake = {} # number of mutants for each sequence that will be generated        
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
        #moreSeqs = cls.set_sequences( fastaFile = ref, refTax = refTax, degap = False, cleanHeader = True, splitChar = '\t', conservative = False, selected = list(seqs2fake.keys())) # Set pf sequences to make more sequences from them
        #print('fake taxa ' + str(sum(seqs2fake.values())))
        #print(seqs2fake)
        moreSeqs = cls.set_sequences( fastaFile = ref, refTax = refTax, cpus = cpus, rank = rank, selected = list(seqs2fake.keys()), degap = False)
        for s in moreSeqs:
            #print(s.header)
            Nposmax = int(estimate_mutations(rank, length = len(s.seq.replace('.','').replace('-',''))))
            #print('Nposmax ' , str(Nposmax))
            Nstrains = int(seqs2fake[s.header]) 
            #print('mutant seqs I want ' + s.header + ' ' + str(Nstrains))
            newS = s.generatemutantASVs(Nstrains = Nstrains, Nposmax = Nposmax, start = 0, end = 50000, include_original = False, by_region = by_region, by_pos = by_pos)
            #print('newS ' + str(len(newS)))
            fakeSeqs.update(newS)                
        return(fakeSeqs)
             
    
    @classmethod
    def init_from_taxa(cls, nASVs, prefix, rank, taxa, taxa_abundances = [], refTax =  '/home/natalia/Projects/natalia/DB/silva.nr_v138/silva.nr_v138.tax', ref = '/home/natalia/Projects/natalia/DB/silva.nr_v138/silva.nr_v138.align', cpus = 20):
        """ Create an environment from a list of user taxa or from random taxa at the selected rank """
        enviro = 'taxas'
        #write_logfile('info', 'SUBSET SEQUENCES', 'You have passed a list of taxa {} from which take sequences and make ASVs'.format(taxa))
        if taxa_abundances:
            taxa_abundances = [i*100 for i in  taxa_abundances]  # To subset 10000 sequences from SILVA
        else:
            # Recreate the list taxa_abundances but with random values
            taxa_abundances = estimate_abundances(len(taxa), 10000) # Do % of each taxa
        #else: # select random taxa instead of sequences from DB
        #    taxa = Enviro.subset_random_taxa(args.refTax, taxlevel, nTaxa = 0)
        #           nTaxa = len(taxa)
        #    taxaAbund = estimate_abundances(nTaxa, args.nASVs)
        taxAbun = dict(zip(taxa, taxa_abundances))
        #taxAbun = tuple(zip(taxa, taxa_abundances))
        write_logfile('info', 'Taxa', 'Subset sequences')
        headers, neededSeqs = cls.subsetSilvaproportions( taxAbun, refTax = refTax, ref = ref, rank = rank )
        write_logfile('info', 'Enviro', 'Select sequences')
        Seqs = cls.set_sequences( fastaFile = ref, refTax = refTax , cpus = cpus, rank = rank, selected = list(headers), degap = False)
        write_logfile('info', 'Taxa', 'Complete sequences')
        if neededSeqs:
            FakeSeqs = cls.make_fake_taxa(neededSeqs, rank = rank, cpus = cpus, ref = ref, refTax = refTax)
            TotalSeqs = Seqs|FakeSeqs
        else:
            TotalSeqs = Seqs
        #print('TotalSeqs' + str(len(TotalSeqs)))
        return( Enviro(prefix, enviro, TotalSeqs, nASVs) ) 
    
    @classmethod   
    def init_from_seqs(cls, prefix, rank, seqs, nASVs, minseqs = 100, refTax =  '/home/natalia/Projects/natalia/DB/silva.nr_v138/silva.nr_v138.tax', ref = '/home/natalia/Projects/natalia/DB/silva.nr_v138/silva.nr_v138.align', cpus = 20):
        enviro = 'seqs'
        #write_logfile('info', 'SEQUENCES PROVIDED', 'Sequences {} are treated as OTUs, mutant ASVs will be generated until complete {} ASVs'.format(seqs, nASVs))
        if seqs:
            if isinstance(seqs, list):
                if nASVs < len(seqs) and seqs:
                    write_logfile('warning', 'SEQUENCES PROVIDED', 'Less total ASVs ({}) than provided sequences ({})'.format(nASVs, len(seqs)))
                    print('Less total ASVs ({}) than provided sequences ({})'.format(nASVs, len(seqs)))
                    print('Reduce the number of sequences, rise the number of ASVs and if not microdiversity is needed: ASVmean = 0')
                    Seqs = cls.set_sequences( fastaFile = ref, refTax = refTax , cpus = cpus, rank = rank, selected = seqs, degap = False)
        else:
            write_logfile('info', 'RANDOM SEQUENCES', 'No align, environment, list of sequences or taxas have been found, random sequences will be subsetted to make ASVs')
            Seqs = cls.set_sequences( fastaFile = ref, refTax = refTax , cpus = cpus, rank = rank, Nrandom = random.randint(minseqs, nASVs), degap = False)
        return( Enviro(prefix, enviro, Seqs, nASVs) ) 
    
    @classmethod
    def init_from_file(cls, prefix, path, inputfile, cpus = 20):
        """ Open and load an align file """
        if os.path.isfile('{}/{}'.format(path, inputfile)): 
            write_logfile('info', 'PREVIOUS ALIGN', 'The file {}/{} will be used as reference'.format(path, inputfile))
            enviro = 'inputfile'
            Seqs = cls.set_sequences( fastaFile = '{}/{}'.format(path, inputfile), cpus = cpus, rank = None, refTax = '', selected = [], degap = False)
        else:
            write_logfile('warning', 'INPUT ALIGN', '{} not found in {}.'.format(inputfile, path))
            exit(-2)
        return( Enviro(prefix, enviro, Seqs, len(Seqs)) )
    
    @classmethod 
    def distances(cls, Seqs, prefix, region, start, end, cpus, complete = False, trimmed = False): # igual se puede combinar con makeASVs
        """ Calculate distances in an specific region (plot a heatmap) and filter sequences based on them. """
        cpus = int(cpus)
        nSeqs = len(Seqs)
        # Trim sequences according to start end position. Return set of Seqs objects
        # https://www.geeksforgeeks.org/copy-python-deep-copy-shallow-copy
        if trimmed:
            cls.multiprocessing_globals = iter([s for s in Seqs])
        else:
            cls.multiprocessing_globals = iter([s.trimregion(start, end) for s in Seqs])
        equivalence = np.array([s.header for s in Seqs]) #"https://docs.python.org/2/library/stdtypes.html#mapping-types-dict: If keys, values and items views are iterated over with no intervening modifications to the dictionary, the order of items will directly correspond"
        if cpus == 1 or nSeqs < 100:
            write_logfile('info', 'CONVERT SEQS', 'Start 1 cpu')
            cls.multiprocessing_globals_seqs = list(map(Sequence.nt2array, cls.multiprocessing_globals)) # Una vez use cls.multiprocessing_globals_seqs se vacía automáticamente el iterable
        else:
            write_logfile('info', 'CONVERT SEQS', 'Start multiprocessing')
            with Pool(cpus) as pool:
                cls.multiprocessing_globals_seqs = list(pool.map(Sequence.nt2array, cls.multiprocessing_globals, chunksize = 50)) #iterable, imap consume iterable y devuelve iterable
        write_logfile('info', 'CONVERT SEQS', 'Finish')
        cls.multiprocessing_globals_combinations = combinations_with_replacement(list(range(nSeqs)),2)
        if cpus == 1 or nSeqs < 100:
            write_logfile('info', 'DISTANCE CALCULATIONS', 'Start 1 cpu')
            arr = np.zeros(shape=(nSeqs, nSeqs), dtype = np.float32)
            distances = map(cls.calc_dist, cls.multiprocessing_globals_combinations)
            for (i, j, d) in map(cls.calc_dist, cls.multiprocessing_globals_combinations):
                arr[i, j] = d 
        else:
            write_logfile('info', 'DISTANCE CALCULATIONS', 'Start multiprocessing')
            arr = np.zeros(shape=(nSeqs, nSeqs), dtype = np.float32)
            with Pool(cpus) as pool:
                #print(list(distances))
                for (i, j, d) in pool.imap(cls.calc_dist, cls.multiprocessing_globals_combinations, chunksize = 50):
                    arr[i, j] = d 
        write_logfile('info', 'DISTANCE CALCULATIONS', 'Finish')
        #arr = np.array(list(map(lambda x, y: x*[0] + y, fill, dist_trian)))
        #print(arr)
        # Make a df with distance calcs for all the sequences
        df = pd.DataFrame(arr, columns = equivalence, index = equivalence)
        # Logical df with False when value ==-1, to remove this columns
        #print('mask')
        mask = df.values!=-1
        # Remove these columns and the corresponding rows
        df = df.drop(index = equivalence[list(np.where(mask==False)[1])], columns = equivalence[list(np.where(mask==False)[1])])
        seqstoremove = list(np.where(mask==False)[1])
        SeqsFilteredheaders = list(np.delete(equivalence, seqstoremove))
        SeqsFiltered = {s for s in Seqs if s.header in SeqsFilteredheaders}
        cls.clear_multiprocessing_globals()
        if complete:
            print('writing')
            write_table(df, title = '{}.distances'.format(prefix), rows = list(df.index) , columns = list(df.columns), dataframe = None)
            return(SeqsFiltered, df)
        else:
            return(SeqsFiltered)
        
    @classmethod
    def clear_multiprocessing_globals(cls): # Clean class variable instead of overwrite it
        cls.multiprocessing_globals_seqs = ''
        cls.multiprocessing_globals_combinations = ''
        cls.multiprocessing_globals = ''
    
    @classmethod
    def calc_dist(cls, args): 
        """ args is a tuple of two index combinations with the index of the sequences between calculate distances """
        d = calculate_distances_cython(cls.multiprocessing_globals_seqs[args[0]], cls.multiprocessing_globals_seqs[args[1]])
        if d != 0:
            return(args[0], args[1], d)
            #return(d)
        elif args[0] == args[1]: # Si es el mismo indice contra si mismo, d = 0 they are identical
            return(args[0], args[1], 0.0)
            #return(0.0)
        else:    
            return(args[0], args[1], -1.0) # Return -1 to filter by this value after and remove those sequences which are identical
            #return(-1.0)
    # Methods for Enviro objects
       # Methods for Enviro objects
    def makeASVs(self, region, start, end, ASVsmean, cutoff , cpus, by_region = None, by_pos = None ): # ADD ASV/OTU calculation
        """ Create fake ASVs from sequences """
        # Subset Seqs, do ASVs, until complete nASVs
        ASVs = set()
        sampleSeqs = set()
        j = 0
        if start and end:
            trimmed = True
        else:
            trimmed = False
        while len(ASVs) < self.nASVs:  # Repeat until complete number of ASVs
            #print('Seqs' + str(len(self.Seqs)))
            newS = random.sample(self.Seqs, 1)[0] # Subset one sequence. Return a list with one element
            if newS.trimregion(start, end).seq:
                if newS.header not in [s.header for s in sampleSeqs]: # if that sequence has not been yet taken 
                    sampleSeqs.add(newS)
                    Nposmax = int(round(cutoff * len(newS.trimregion(start, end).deGap().seq)))
                    #print(newS.seq)
                    #print(len(newS.seq))
                    newSasvs = newS.generatemutantASVs(Nstrains = None, Nmean = ASVsmean, Nposmax = Nposmax, start = start, end = end, by_region = by_region, by_pos = by_pos)
                    #Check distances
                    ASVsdiff = self.distances(Seqs = newSasvs, prefix = self.prefix, region = region, start = start, end = end, cpus = cpus, trimmed = trimmed)
                    #write_logfile('warning', 'SUBSET RANDOM SEQS', 'len(ASVs) {}\tlen(ASVsdiff) {}'.format(len(newSavs), len(ASVsdiff)))
                    i=0
                    while len(ASVsdiff) < len(newSasvs): # If one sequence have been removed 'cause it is identical to another one, repeat to take ASVs from that sequence, otherwise some taxa can be minusvalorated
                        newSavs = newS.generatemutantASVs(Nstrains = None, Nmean = ASVsmean, Nposmax = 45, start = start, end = end, by_region = by_region, by_pos = by_pos)
                        #Check distances
                        ASVsdiff = self.distances(Seqs = newSasvs, prefix = self.prefix, region = region, start = start, end = end, cpus = cpus, trimmed = trimmed)
                        i = i+1
                        if i > 10:
                            #write_logfile('warning', 'SUBSET RANDOM SEQS2', 'len(ASVsdiff) {}\tlen(newSavs) {}'.format(len(newSavs), len(ASVsdiff)))
                            break
                    ASVs.update(ASVsdiff)
                else:
                    continue
            else:
                j = j + 1
            maxj = 0.1*len(ASVs)
            if j > maxj:
                write_logfile('warning', 'ENVIRONMENT', 'More than the 10% of the total of requested sequences that we have subset are empty in the region included between your start/end positions, so we have decided to stop iterating. Please, remember that these positions are referred to the SILVA alignment and not the length of the 16S')
                exit()
        ASVsdiff, distdf = self.distances(Seqs = ASVs, prefix = self.prefix, region = region, start = start, end = end, cpus = cpus, complete = True, trimmed = trimmed) 
        #ASVsdiffselected = set(random.sample({s for s in ASVs if s.header not in [so.header for so in SeqsSilvaselected]}, args.nASVs-len(SeqsSilvaselected))) # No matter if original sequences are included or not
        #print('ASVsdiff ' + str(len(ASVsdiff)))
        #print('nASVs ' + str(self.nASVs))
        ASVsdiffselected = set(random.sample(ASVsdiff, self.nASVs))
        #self.Seqs = SeqsSilvaselected.union(ASVsdiffselected) # No matter if original sequences are included or not
        self.Seqs = ASVsdiffselected # I REALLY WANT TO MODIFY IT, I want to remove extra Sequences
        write_logfile('info', 'ENVIRONMENT', 'Plotting distances heatmap')
        self.plot_distances(df = distdf, region = region, text = None, symmetric = True)
        return(self)# I really want to modify it 
    
    #### PLOTS & OUTPUTS
    def plot_distances(self, df, region = '16S',  text = None, symmetric = None):
        """ Plot distances among selected sequences """
        # Update df with the final sequences that finally have been included:
        Seqs_headers = [s.header for s in self.Seqs]
        shortdf = df.loc[ [idseq for idseq in df.index if idseq in Seqs_headers] , [idseq for idseq in df.columns if idseq in Seqs_headers] ] # select columns and rows
        plot_heatmap(shortdf.T, outputDir = os.getcwd(), title = '{}.{}.distances'.format(self.prefix, region), zmin = 0, zmax = 1, legendtitle = 'distance', symmetric = symmetric)
    
    def plot_taxonomy(self, rank):
        """ Barplot stacked and percent of the taxonomy of the sequences included in the mock """
        taxa = [s.tax.split(';')[rank] for s in self.Seqs]
        #taxa = [s.tax for s in self.Seqs]
        taxAbun =  dict(Counter(taxa))
        # Create a df of taxonomy and counts
        df = pd.DataFrame.from_dict(taxAbun, orient='index', columns=['counts']).transpose()
        #barplotpc(df, outputDir = os.getcwd(), title = '{}.pctax'.format(self.prefix),  ylab = 'taxon pc %', xlab = 'Mock', textSize = 8 )
        # Sort by value
        df.sort_values(by = 'counts', axis = 1, ascending = False, inplace = True)
        barplotpc(df, outputDir = os.getcwd(), title = '{}.pctax'.format(self.prefix), ylab = 'taxon pc %', xlab = 'Taxonomy at the rank level')
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
    
    @classmethod
    def set_sequences( cls, fastaFile, refTax, rank, cpus = 1, Nrandom = 0, selected = [], degap = False):
        """ Make a list of sequence objects (or select some based on the header) from fasta """
        cls.rank = rank
        cls.degap = degap
        if not refTax:
            cls.notax = True
        else:
            cls.silva_taxa = loadTaxa(refTax = refTax)
            if Nrandom != 0:
                cls.selected = random.sample(cls.silva_taxa.keys(), Nrandom)
            else:
                cls.selected = selected
        with Pool(cpus) as pool:
            with open(fastaFile, 'r') as fasta:
                Seqs = set(list(filter(None, pool.imap(cls.validate_Sequence, itertools.zip_longest(*[fasta]*2)))))
        return(Seqs)
    
    @classmethod
    def validate_Sequence(cls, combo):
        header = combo[0]
        seq = combo[1]
        #if cleanHeader: # if it is other DB distict from SILVA
        #    header = simplifyString(header.lstrip('>').rstrip('\n'))
        #else:
        #    header = header.lstrip('>').rstrip('\n')
        header = simplifyString(header.lstrip('>').rstrip('\n')) # assuming SILVA db
        if not 'N' in seq: #checkpoint: some sequences from silva contain N, exclude those sequences
            Seq = None
            if not cls.selected: # select all the sequences from the fasta file
                #write_logfile('debug', 'Sequence.set_Sequences', cls.selected)
                if not cls.notax:
                    Seq = Sequence(header, seq.rstrip('\n'), Sequence.assign_taxonomy(header, cls.silva_taxa), 0)
                else:
                    Seq = Sequence(header, seq.rstrip('\n'), None, 0) # No taxonomic annotation
                #write_logfile('debug', 'Sequence.set_Sequences NOT selected ', Seq.header)
            else:
                #write_logfile('debug', 'Sequence.set_Sequences', cls.selected)
                if header in cls.selected: #add sequence to the set
                    Seq = Sequence(header, seq.rstrip('\n'), Sequence.assign_taxonomy(header, cls.silva_taxa), 0)
                #write_logfile('debug', 'Sequence.set_Sequences SELECTED', Seq.header)
            if Seq:
                if cls.degap:
                    Seq = Seq.deGap()
                    write_logfile('debug', 'Sequence.set_Sequences ADD', Seq.header)
                return(Seq)
    
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


    @classmethod
    def checkN(cls, combo):
        header = combo[0]
        seq = combo[1]
        header = simplifyString(header.lstrip('>').rstrip('\n')) # assuming SILVA db
        if not 'N' in seq: #checkpoint: some sequences from silva contain N, exclude those sequences
            return(header)






    @classmethod
    def subsetSilvaproportions(cls, taxAbun, refTax =  '/home/natalia/Projects/natalia/DB/silva.nr_v138/silva.nr_v138.tax', ref = '/home/natalia/Projects/natalia/DB/silva.nr_v138/silva.nr_v138.align', rank = 3, cpus = 12):
        """ From a dict of taxas {'tax':Abun} take, sequences from silva """
        rank = rank
        silva_taxa = loadTaxa(refTax = refTax, rank = rank)
        #silvaSeqs = list(fasta2dict(ref).keys())
        with Pool(cpus) as pool:
            with open(ref, 'r') as fasta:
                silvaSeqs = set(list(pool.map(cls.checkN, itertools.zip_longest(*[fasta]*2))))
        for seq in list(silva_taxa.keys()):
            if not seq in silvaSeqs:
                silva_taxa.pop(seq)
        taxonSeqs = reversedict(silva_taxa) # {taxa:[seq,...]}
        # Some sequences from silva have 'N', we do not want to include them, so exclude these sequences
         # Only sequences without 'N'
        # Check good sequences in taxonSeqs, remove those that does not appear.
        
        subsetSeqs = set()
        moreSeqs = set() # set of species for which more sequences are required {(abun,(seq1,seq2)),(abun,(seq3,seq4,seq5))}
        for tx, abun in taxAbun.items(): # si taxAbun es una tupla: for tx, abun in taxAbun:
            if tx in list(taxonSeqs.keys()):
                if abun > len(taxonSeqs[tx]):
                    #write_logfile('info', 'SUBSET RANDOM SEQS', '{} Ask for {}, but only {}. Mutant representative sequences will be created'.format(tx, abun, len(taxonSeqs[tx])))
                    moreSeqs.add((abun - len(taxonSeqs[tx]) , tuple(taxonSeqs[tx])))
                    abun = len(taxonSeqs[tx])
                # If there are not enough sequences, create more at  the given taxonomic level!!!!!!!!! MODIFY IN FUTURE
                seq = random.sample(taxonSeqs[tx], abun)
                #print(seq)
                # If the sequence has been added in a previous loop, it does not matter 'cause it is a set, no duplicates
                subsetSeqs = subsetSeqs.union(set(seq))
            else:
                write_logfile('error', 'SUBSET RANDOM SEQS', '{} not present'.format(tx))
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














