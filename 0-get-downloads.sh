#!/usr/bin/bash

# Host data 

## top link: https://protists.ensembl.org/Emiliania_huxleyi/Info/Index
## download host whole genome fasta

base_dir="/home/bbitarello" # change as needed
cd ${base_dir}

#download ensembl protists v60 (when project began) and v62 (when project ended)
for vers in 60 62;
do
mkdir -p ${base_dir}/scratch/CCMP1516/v${vers}/
cd ${base_dir}/scratch/CCMP1516/v${vers}

wget "http://ftp.ensemblgenomes.org/pub/protists/release-${vers}/fasta/emiliania_huxleyi/dna/Emiliania_huxleyi.Emiliana_huxleyi_CCMP1516_main_genome_assembly_v1.0.dna.toplevel.fa.gz"

## download host cdna fasta
wget "https://ftp.ensemblgenomes.ebi.ac.uk/pub/protists/release-${vers}/fasta/emiliania_huxleyi/cdna/Emiliania_huxleyi.Emiliana_huxleyi_CCMP1516_main_genome_assembly_v1.0.cdna.all.fa.gz"

## download host gff3
wget "https://ftp.ensemblgenomes.ebi.ac.uk/pub/protists/release-${vers}/gff3/emiliania_huxleyi/Emiliania_huxleyi.Emiliana_huxleyi_CCMP1516_main_genome_assembly_v1.0.${vers}.gff3.gz"

## download ncRNA
wget "http://ftp.ensemblgenomes.org/pub/protists/release-${vers}/fasta/emiliania_huxleyi/ncrna/Emiliania_huxleyi.Emiliana_huxleyi_CCMP1516_main_genome_assembly_v1.0.ncrna.fa.gz"

wget "https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-${vers}/protists/uniprot_report_EnsemblProtists.txt"
done

#check that nothing has changed
## diff <(grep Emiliania ~/scratch/CCMP1516/v60/uniprot_report_EnsemblProtists.txt) <(grep Emiliania ~/scratch/CCMP1516/v62/uniprot_report_EnsemblProtists.txt)
## zdiff Emiliania_huxleyi.Emiliana_huxleyi_CCMP1516_main_genome_assembly_v1.0.dna.toplevel.fa.gz ../v60/Emiliania_huxleyi.Emiliana_huxleyi_CCMP1516_main_genome_assembly_v1.0.dna.toplevel.fa.gz
## zdiff Emiliania_huxleyi.Emiliana_huxleyi_CCMP1516_main_genome_assembly_v1.0.cdna.all.fa.gz ../v62/Emiliania_huxleyi.Emiliana_huxleyi_CCMP1516_main_genome_assembly_v1.0.cdna.all.fa.gz > zdiff.out
#grep ">" zdiff.out|wc #12

mkdir -p ~/scratch/CCMP1516/ncbi/
cd ~/scratch/CCMP1516/ncbi/
#https://ngdc.cncb.ac.cn/p10k/sample/P10K-NCBI-001245

wget "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/372/725/GCA_000372725.1_Emiliana_huxleyi_CCMP1516_main_genome_assembly_v1.0/GCA_000372725.1_Emiliana_huxleyi_CCMP1516_main_genome_assembly_v1.0_genomic.gff.gz"

wget "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/372/725/GCA_000372725.1_Emiliana_huxleyi_CCMP1516_main_genome_assembly_v1.0/GCA_000372725.1_Emiliana_huxleyi_CCMP1516_main_genome_assembly_v1.0_protein.faa.gz"

wget "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/372/725/GCA_000372725.1_Emiliana_huxleyi_CCMP1516_main_genome_assembly_v1.0/GCA_000372725.1_Emiliana_huxleyi_CCMP1516_main_genome_assembly_v1.0_cds_from_genomic.fna.gz"

#check nr of seqs
#zgrep ">" GCA_000372725.1_Emiliana_huxleyi_CCMP1516_main_genome_assembly_v1.0_cds_from_genomic.fna.gz|wc -l# 38559
## EhV201 virus data

## install ncbi_datasets-cli: https://github.com/ncbi/datasets
conda create -n ncbi_datasets

conda activate ncbi_datasets

conda install -c conda-forge ncbi-datasets-cli

datasets download genome accession GCA_010974645.1 --include gff3,rna,cds,protein,genome,seq-report

mkdir -p ~/scratch/EhV201/

unzip ncbi_dataset.zip -d ~/scratch/EhV201/

rm ncbi_dataset.zip

conda deactivate
