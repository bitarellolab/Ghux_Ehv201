#!/usr/bin/bash

# Host data 

## top link: https://protists.ensembl.org/Emiliania_huxleyi/Info/Index
## download host whole genome fasta
mkdir -p ~/scratch/CCMP1516
cd ~/scratch/CCMP1516
wget http://ftp.ensemblgenomes.org/pub/protists/release-60/fasta/emiliania_huxleyi/dna/Emiliania_huxleyi.Emiliana_huxleyi_CCMP1516_main_genome_assembly_v1.0.dna.toplevel.fa.gz

## download host cdna fasta
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/protists/release-60/fasta/emiliania_huxleyi/cdna/Emiliania_huxleyi.Emiliana_huxleyi_CCMP1516_main_genome_assembly_v1.0.cdna.all.fa.gz

## download host gff3
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/protists/release-60/gff3/emiliania_huxleyi/Emiliania_huxleyi.Emiliana_huxleyi_CCMP1516_main_genome_assembly_v1.0.60.gff3.gz

## download ncRNA
wget http://ftp.ensemblgenomes.org/pub/protists/release-60/fasta/emiliania_huxleyi/ncrna/Emiliania_huxleyi.Emiliana_huxleyi_CCMP1516_main_genome_assembly_v1.0.ncrna.fa.gz

## virus data

## install ncbi_datasets-cli: https://github.com/ncbi/datasets
conda create -n ncbi_datasets

conda activate ncbi_datasets

conda install -c conda-forge ncbi-datasets-cli

datasets download genome accession GCA_010974645.1 --include gff3,rna,cds,protein,genome,seq-report

unzip ncbi_dataset.zip -d ~/scratch/EhV201

rm ncbi_dataset.zip

conda deactivate
