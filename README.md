# MotiveSearch

## Introduction

## Workflow
1. Identify differentially expressed genes. (e.g. RNAseq)
2. Aquire genomic sequence of host organism and sequences of promoters.
3. Specify parameters in config.txt. An example is given in the "Input files" section.
4. Determine most frequent kmers in the promoters of the differentially expressed genes.
5. Determine baseline frequency of those kmers in the entire genome.
6. Cluster kmers that are found significantly more frequent in the used promoters than in the genome.


## Usage
The script offers two 

python script.py -m PrepareGenome -i genome.fasta

python script.py -m AnalysePromoters -i counttable.txt -t counttable
## Input files