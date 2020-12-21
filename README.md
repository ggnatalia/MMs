# MMs
Software for generating Metagenomic samples of realistic Microbial Mock communities

usage: makemocks.py [-h] -m MOCKNAME -o OUTPUT [-N NSAMPLES] [-nASVs NASVS] [-H SHANNON] [-r RANK] [-ASVsmean ASVSMEAN] [-env ENVIRO] [--taxa TAXA [TAXA ...]] [--seqs SEQS [SEQS ...]] [--minseqs MINSEQS] [--taxaAbund TAXAABUND [TAXAABUND ...]] [--input-file INPUT_FILE]
                    [-ref REF] [-refTax REFTAX] [-refEnv REFENVIRO] [-s START] [-e END] [--region REGION] [--cutoff CUTOFF] [--cpus CPUS] [--force-overwrite] [--error-model ERROR_MODEL] [--read-length READ_LENGTH [READ_LENGTH ...]]
                    [--insert-size INSERT_SIZE [INSERT_SIZE ...]] [--sequences-files SEQUENCES_FILES [SEQUENCES_FILES ...]] [--abundance-files ABUNDANCE_FILES [ABUNDANCE_FILES ...]] [--repeat-inSilicoSeqs-autocomplete] [--reads READS] [--alpha ALPHA] [--pstr0 PSTR0]
                    [--size SIZE] [--figsize FIGSIZE]

Process input files

optional arguments:
  -h, --help            show this help message and exit
  -m MOCKNAME, --mockName MOCKNAME
                        Name of the mock e.g. mock1.
  -o OUTPUT, --output OUTPUT
                        Name of the output directory e.g. aquatic.
  -N NSAMPLES, --nSamples NSAMPLES
                        Number of samples to generate. Default 5 samples per mock.
  -nASVs NASVS, --nASVs NASVS
                        Number of ASVs.
  -H SHANNON, --shannon SHANNON
                        Shannon diversity Index. Default 3.
  -r RANK, --rank RANK  Rank to subset taxa: phlyum, order, class, family, genus.
  -ASVsmean ASVSMEAN, --ASVsmean ASVSMEAN
                        Mean of the number of mutants per reference.
  -env ENVIRO, --enviro ENVIRO
                        Select one of the predefined environments to simulate the mock. See the file SpeciesperEnviro.tsv for options.
  --taxa TAXA [TAXA ...]
                        List of taxa to generate the mock.
  --seqs SEQS [SEQS ...]
                        List of sequences's header.
  --minseqs MINSEQS     Minimun number of sequences to randomly extract from DB. Default 5.
  --taxaAbund TAXAABUND [TAXAABUND ...]
                        List of taxa abundances in percentages.
  --input-file INPUT_FILE
                        Provide a previous fasta/align file.
  -ref REF, --ref REF   SILVA alignment reference formatted for using mothur e.g. silva.nr_v138.align.
  -refTax REFTAX, --refTax REFTAX
                        SILVA taxonomy reference formatted for using mothur e.g. silva.nr_v138.tax.
  -refEnv REFENVIRO, --refEnviro REFENVIRO
                        Table with most frequent species per environment e. g. SpeciesperEnviro.tsv.
  -s START, --start START
                        Start position in the SILVA alignment reference. Default 1.
  -e END, --end END     End position in the SILVA alignment reference. Default 50000.
  --region REGION       Region of the 16S rRNA gene.
  --cutoff CUTOFF       Filter ASVs depending on sequence similarity. Without this flag, sequences can be identical in the studied region.
  --cpus CPUS           Number of threads.
  --force-overwrite     Force overwrite if the output directory already exists.
  --error-model ERROR_MODEL
                        Mode to simulate reads in InSilicoSeqs.
  --read-length READ_LENGTH [READ_LENGTH ...]
                        Read length of the simulated reads. Only available for basic and perfect modes.
  --insert-size INSERT_SIZE [INSERT_SIZE ...]
                        Insert size of the simulated reads. Only available for basic and perfect modes.
  --sequences-files SEQUENCES_FILES [SEQUENCES_FILES ...]
                        Sequences files for InSilicoSeqs: <projectName><mockName>.<sampleName>.sequences16S.fasta.
  --abundance-files ABUNDANCE_FILES [ABUNDANCE_FILES ...]
                        Abundances files for InSilicoSeqs: <projectName><mockName>.<sampleName>.abundances.
  --repeat-inSilicoSeqs-autocomplete
                        Repeat reads simulation using previous files.
  --reads READS         Number of reads approximately in total, counting both pairs.
  --alpha ALPHA         Correlation Matrix: Probability that a coefficient is zero. Larger values enforce more sparsity.
  --pstr0 PSTR0         ZINBD: Probability of structure 0.
  --size SIZE           ZINBD: Size - dispersion of ZINBD.
  --figsize FIGSIZE     Size of plots.

