#!/usr/bin/bash

# Ehux

##top link: https://protists.ensembl.org/Emiliania_huxleyi/Info/Index
##download host whole genome fasta

base_dir="/home/bbitarello" # change as needed

cd ${base_dir}

#download ensembl protists v60 (when project began) and v62 (when project ended)
for vers in 60 62;
do
mkdir -p ${base_dir}/scratch/CCMP1516/v${vers}/
cd ${base_dir}/scratch/CCMP1516/v${vers}

wget "http://ftp.ensemblgenomes.org/pub/protists/release-${vers}/fasta/emiliania_huxleyi/dna/Emiliania_huxleyi.Emiliana_huxleyi_CCMP1516_main_genome_assembly_v1.0.dna.toplevel.fa.gz"

##download host cdna fasta
wget "https://ftp.ensemblgenomes.ebi.ac.uk/pub/protists/release-${vers}/fasta/emiliania_huxleyi/cdna/Emiliania_huxleyi.Emiliana_huxleyi_CCMP1516_main_genome_assembly_v1.0.cdna.all.fa.gz"

##download host gff3
wget "https://ftp.ensemblgenomes.ebi.ac.uk/pub/protists/release-${vers}/gff3/emiliania_huxleyi/Emiliania_huxleyi.Emiliana_huxleyi_CCMP1516_main_genome_assembly_v1.0.${vers}.gff3.gz"

##download ncRNA
wget "http://ftp.ensemblgenomes.org/pub/protists/release-${vers}/fasta/emiliania_huxleyi/ncrna/Emiliania_huxleyi.Emiliana_huxleyi_CCMP1516_main_genome_assembly_v1.0.ncrna.fa.gz"

## download uniprot report
wget "https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-${vers}/protists/uniprot_report_EnsemblProtists.txt"
done

##check that nothing has changed
## diff <(grep Emiliania ~/scratch/CCMP1516/v60/uniprot_report_EnsemblProtists.txt) <(grep Emiliania ~/scratch/CCMP1516/v62/uniprot_report_EnsemblProtists.txt)
## zdiff Emiliania_huxleyi.Emiliana_huxleyi_CCMP1516_main_genome_assembly_v1.0.dna.toplevel.fa.gz ../v60/Emiliania_huxleyi.Emiliana_huxleyi_CCMP1516_main_genome_assembly_v1.0.dna.toplevel.fa.gz
## zdiff Emiliania_huxleyi.Emiliana_huxleyi_CCMP1516_main_genome_assembly_v1.0.cdna.all.fa.gz ../v62/Emiliania_huxleyi.Emiliana_huxleyi_CCMP1516_main_genome_assembly_v1.0.cdna.all.fa.gz > zdiff.out
##grep ">" zdiff.out|wc #12

conda activate ncbi_datasets

datasets download genome accession GCF_000372725.1 --include gff3,rna,cds,protein,genome,seq-report

unzip ncbi_dataset.zip -d ${base_dir}/scratch/CCMP1516/

rm ncbi_dataset.zip

cd ${base_dir}/scratch/CCMP1516/ncbi_dataset/data/GCF_000372725.1/

#check how many "genes"

grep ">" cds_from_genomic.fna |wc -l #38559

#check how many "hypothetical protein" annots:

grep "hypothetical protein" cds_from_genomic.fna |wc -l # 35788

#check nr of seqs
#zgrep ">" GCA_000372725.1_Emiliana_huxleyi_CCMP1516_main_genome_assembly_v1.0_cds_from_genomic.fna.gz|wc -l# 38559
## EhV201 virus data

# EhV201

mkdir -p ${base_dir}/scratch/EhV201/

conda activate ncbi_datasets

datasets download genome accession GCA_010974645.1 --include gff3,rna,cds,protein,genome,seq-report

unzip ncbi_dataset.zip -d ~/scratch/EhV201/

rm ncbi_dataset.zip

conda deactivate
