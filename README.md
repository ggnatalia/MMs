# MMs
Software for generating Metagenomic samples of realistic Microbial Mock communities

usage: makemocks.py [-h] -m MOCKNAME -o OUTPUT -N NSAMPLES -nASVs NASVS [-H SHANNON] [-r RANK] [-env ENVIRO] [--cpus CPUS] [--force-overwrite] 

Process input files

optional arguments:
  -h, --help            show this help message and exit    
  
  -m, --mockName: Name of the mock e.g. mock1.    
  
  -o, --output:     Name of the output directory e.g. aquatic.   
  
  -N, --nSamples: Number of samples to generate. Default 5 samples per mock.   
  
  -nASVs, --nASVs:   Number of ASVs.   
  
  -H, --shannon:   Shannon diversity Index. Default 3.    
  
  -r, --rank   Rank to subset taxa: phlyum, order, class, family, genus.    
  
  -ASVsmean, --ASVsmean: Mean of the number of mutants per reference.   
  
  -env, --enviro: Select one of the predefined environments to simulate the mock. See the file SpeciesperEnviro.tsv for options.    
  
  --taxa : List of taxa to generate the mock.    
  
  --seqs :  List of sequences's header.    
  
  --minseqs : Minimun number of sequences to randomly extract from DB. Default 5.    
  
  --taxaAbund: List of taxa abundances in percentages.    
  
  --input-file: Provide a previous fasta/align file.    
  
  -ref, --ref: SILVA alignment reference formatted for using mothur e.g. silva.nr_v138.align.   
  
  -refTax, --refTax: SILVA taxonomy reference formatted for using mothur e.g. silva.nr_v138.tax.    
  
  -refEnv, --refEnviro: Table with most frequent species per environment e. g. SpeciesperEnviro.tsv.   
  
  -s, --start : Start position in the SILVA alignment reference. Default 1.    
  
  -e, --end : End position in the SILVA alignment reference. Default 50000.  
  
  --region : Region of the 16S rRNA gene.      
  
  --cutoff: Filter ASVs depending on sequence similarity. Without this flag, sequences can be identical in the studied region.    
  
  --cpus : Number of threads.    
  
  --force-overwrite : Force overwrite if the output directory already exists.  
  
  --error-model : Mode to simulate reads in InSilicoSeqs.    
  
  --read-length : Read length of the simulated reads. Only available for basic and perfect modes.   
  
  --insert-size : Insert size of the simulated reads. Only available for basic and perfect modes.   
  
  --sequences-files : Sequences files for InSilicoSeqs: <projectName><mockName>.<sampleName>.sequences16S.fasta.    
  
  --abundance-files : Abundances files for InSilicoSeqs: <projectName><mockName>.<sampleName>.abundances.    
  
  --repeat-inSilicoSeqs-autocomplete : Repeat reads simulation using previous files.   
  
  --reads : Number of reads approximately in total, counting both pairs.    
  
  --alpha : Correlation Matrix: Probability that a coefficient is zero. Larger values enforce more sparsity.    
  
  --pstr0 : ZINBD: Probability of structure 0.    
  
  --size : ZINBD: Size - dispersion of ZINBD.    
  
  --figsize : Size of plots.    

