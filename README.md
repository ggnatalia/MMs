# M&Ms
Software for generating Metagenomic samples of realistic Microbial Mock communities

## 1. What is M&Ms?  

M&Ms is a user-friendly open-source bioinformatic tool to produce realistic amplicon datasets from reference sequences, based on pragmatic ecological parameters.   
It creates sequence libraries for 'in silico' microbial communities with user-controlled richness, evenness, microdiversity, and source environment. It also provides additional figures and files with extensive details on how each synthetic community is composed.
M&Ms allows the user to generate simple to complex read datasets based on real parameters that can be used in developing bioinformatic software or in benchmarking current tools.  

Simulation with M&Ms has four stages (Fig.1):  
![](https://github.com/ggnatalia/MMs/blob/images/Figure_1.png)
*Basic workflow of M&Ms.*  

  1.- Selection of the community members  
  2.- Microdiversity simulation  
  3.- Microbial abundance distribution assignment    
  4.- Sequencing data simulation using InSilicoSeqs (Gourlé *et al*., 2019) to produce realistic Illumina reads.    

## 2. Installation

M&Ms is intended to be run in a x86_64 Linux OS (tested in Ubuntu and CentOS). The easiest way to install it is by using conda (https://anaconda.org/ggnatalia/mms).

  1.-Create the enviroment:
    
    conda create -y --name MMs python=3.7

  This will create a new conda environment named MMs, which must then be activated.

  2.-Activate it:
    
    conda activate MMs

  3.-Install M&Ms package and its dependencies:
    
    conda install -c anaconda -c bioconda -c conda-forge -c ggnatalia mms

  Once M&Ms have finished, deactivate the environment:
    
    conda deactivate

  To remove the M&Ms enviroment:
    
    conda remove --name MMs --all

When using conda, all the scripts from the M&Ms distribution will be available on $PATH.


## 3. Downloading the databases or using a previous database

M&Ms uses two databases: the mothur formatted SILVA database (Schloss *et al*., 2009; Yilmaz *et al*., 2014) and a genus per environment matrix (Tamames et al., 2016). The script make_databases.py can be run to download from source the databases or link the users' databases to M&Ms. 
Databases will be in the /path/to/MMs/DB
    
    # Usage
    # Downloading SILVA databases:
    /path/to/MMs/bin/make_databases.py
    # Using a previous downloaded SILVA database
    /path/to/MMs/bin/make_databases.py -ref /path/to/silva.nr_v138/silva.nr_v138.align -refTax /path/to/silva.nr_v138/silva.nr_v138.tax


 ## 4. Execution and running scripts

### Scripts location
The main script composing M&Ms can be found in the /path/to/MMs/makemocks.py. The scripts composing M&Ms classes and modules are in the /path/to/MMs/libs
Other utility scripts can be found in the /path/to/MMs/utils directory. 

### Execution
The command for running M&Ms has the following syntax:

    makemocks.py -m <mockname> -o <output> -N <samplesnumber> -nASVs <ASVsnumber> -H <Shannonindex> -r <reads> <options> 

##### Mandatory parameters
    
    -m, --mockName: Name of the mock e.g. mock1    
    -o, --output: Name of the output directory e.g. aquatic    
    -N, --nSamples: Number of samples to generate        
    -H, --shannon: Shannon diversity Index   
    -r, --reads: Number of reads
    
M&Ms accepts one of the following five independent types of inputs depending on the step from which starts M&Ms:

- To generate a mock microbial communitity from the taxa selection (step 1):

  First, the user must set the desired number of sequences with:
   
      -nASVs, --nASVs: Number of different sequences
    
  Then, choose one of the following options:

   a) Use a list of sequences names from SILVA: 
   
       -seqs, --seqs: List of sequences' names from SILVA db
   
   b) Use a list of taxa (and optionally their abundances):
   
      -tx, --taxa: List of taxa
      #Optionally, the user can provide a second argument with the abundances of each taxa:
      -txAbund, --taxaAbund: List of abundances of taxa
    
   c) Select an environment   

      -env, --enviro: Let the user simulate an environmental mock by choosing species from a specific environment

- To recreate microbial abundances and sequencing simulation (step 3):   
   d) Use a FASTA file  
    
      --inputfile: The user can provide an align file or a FASTA file without creating it from scratch.

- To repeat just sequencing simulation (step 4):   
   e) Repeat just the sequencing simulation.  
   
      # The user has two options, provide the names of the FASTA files and the abundances' files with the following options:
      --sequences_files: List of files to obtain reads: sequences <projectName><mockName>.<sampleName>.sequences16S.fasta. Same order that paired abundance files
      --abundance_files: List of files to obtain reads: abundances<projectName><mockName>.<sampleName>.abundances. Same order that paired sequence files
      # Or do it automatically with:
      --repeat_InSilicoSeqs_autocomplete: If True use mock community and samples directly from the directory, without writing one by one the files
   
##### Options

 - Reference files (if the user has not previously set the databases): 
      
       -ref, --ref: SILVA alignment reference
       -refTax, --refTax: SILVA alignment TAX reference    
       -refEnv, --refEnviro: Taxa per environment reference
 
 - Select the taxonomic level (not needed if the user recreates microbial communities from specific environments):
        
       --rank: Rank to subset taxa: phlyum, order, class, family, genus. Default phylum
 
 - Select the microdiversity per species:
 
       -ASVsmean,--ASVsmean: Mean of mutant sequences (microdiversity) per reference sequence. Default 5.
       --cutoff: Cutoff to identify species and its microdiversity (distance calculations). Default 0.03

 - InSilicoSeqs parameters:    

       --errormodel: Mode to generate InSilicoSeqs (see InSilicoSeqs documentation: https://insilicoseq.readthedocs.io/en/latest/index.html): basic, perfect,    MiSeq, HiSeq, NovaSeq
       --InSilicoparams: List of two elements in the following order: Read length, insert size. Default = [150, 200]. 
  
  - Region of the 16S to consider: 
    
        -s, --start: SILVA alignment reference positions-START. Default 1. 1-based
        -e, --end: SILVA alignment reference positions-END. Default 50000. 1-based
        --region: Name of the studied region. Default '16S'    
   
   - Making a mock randomly 
    
         --minseqs: Minimun number of sequences to extract from DB. Default 5
    
   - Performance
 
         --cpus: Number of threads
         --force-overwrite: Force overwrite if the output directory already exists
         --just_taxa_selection: If True: do the selection of the sequences and stop
    
   - Other customizable parameters
    
         --alpha: Correlation Matrix: Probability for a coefficient to be zero. Larger values enforce more sparsity. Default 0.9
         --pstr0: ZINBD: Probability of structure 0. Default 0.2
         --size: 'ZINBD: Size - dispersion of ZINBD. Default 1
