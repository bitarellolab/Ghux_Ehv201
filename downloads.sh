#!/usr/bin/bash

# Host data 

## top link: https://protists.ensembl.org/Emiliania_huxleyi/Info/Index
## download host whole genome fasta
mkdir -p ~/scratch/CCMP1516
cd ~/scratch/CCMP1516
wget http://ftp.ensemblgenomes.org/pub/protists/release-61/fasta/emiliania_huxleyi/dna/Emiliania_huxleyi.Emiliana_huxleyi_CCMP1516_main_genome_assembly_v1.0.dna.toplevel.fa.gz

## download host cdna fasta
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/protists/release-61/fasta/emiliania_huxleyi/cdna/Emiliania_huxleyi.Emiliana_huxleyi_CCMP1516_main_genome_assembly_v1.0.cdna.all.fa.gz

## download host gff3
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/protists/release-61/gff3/emiliania_huxleyi/Emiliania_huxleyi.Emiliana_huxleyi_CCMP1516_main_genome_assembly_v1.0.60.gff3.gz

## download ncRNA
wget http://ftp.ensemblgenomes.org/pub/protists/release-61/fasta/emiliania_huxleyi/ncrna/Emiliania_huxleyi.Emiliana_huxleyi_CCMP1516_main_genome_assembly_v1.0.ncrna.fa.gz

## ena
wget http://ftp.ensemblgenomes.org/pub/protists/release-61/tsv/emiliania_huxleyi/Emiliania_huxleyi.Emiliana_huxleyi_CCMP1516_main_genome_assembly_v1.0.61.ena.tsv.gz
wget http://ftp.ensemblgenomes.org/pub/protists/release-61/tsv/emiliania_huxleyi/README_ENA.tsv
wget http://ftp.ensemblgenomes.org/pub/protists/release-61/tsv/emiliania_huxleyi/Emiliania_huxleyi.Emiliana_huxleyi_CCMP1516_main_genome_assembly_v1.0.61.uniprot.tsv.gz
wget http://ftp.ensemblgenomes.org/pub/protists/release-61/tsv/emiliania_huxleyi/README_entrez.tsv
wget http://ftp.ensemblgenomes.org/pub/protists/release-61/tsv/emiliania_huxleyi/README_uniprot.tsv
wget http://ftp.ensemblgenomes.org/pub/protists/release-61/tsv/emiliania_huxleyi/README_refseq.tsv
## virus data

## install ncbi_datasets-cli: https://github.com/ncbi/datasets
cd ~/scratch/

conda create -n ncbi_datasets

conda activate ncbi_datasets

conda install -c conda-forge ncbi-datasets-cli

datasets download genome accession GCA_010974645.1 --include gff3,rna,cds,protein,genome,seq-report

unzip ncbi_dataset.zip -d ~/scratch/EhV201

rm ncbi_dataset.zip

mv ncbi_dataset/data/* .

rm -rf ncbi_dataset/

mv GCA_010974645.1/* .

rm -rf GCA_010974645.1

conda deactivate


# check for differences from previous ensembl protist release

zdiff Emiliania_huxleyi.Emiliana_huxleyi_CCMP1516_main_genome_assembly_v1.0.cdna.all.fa.gz ~/projects/rna-seq-host-virus/recreate_upenn/CCMP1516/Emiliania_huxleyi.Emiliana_huxleyi_CCMP1516_main_genome_assembly_v1.0.cdna.all.fa.gz
#no diff

zdiff Emiliania_huxleyi.Emiliana_huxleyi_CCMP1516_main_genome_assembly_v1.0.dna.toplevel.fa.gz ~/projects/rna-seq-host-virus/recreate_upenn/CCMP1516/Emiliania_huxleyi.Emiliana_huxleyi_CCMP1516_main_genome_assembly_v1.0.dna.toplevel.fa.gz
#no diff

conda create -n seqkit
conda activate seqkit
conda install -c bioconda seqkit
seqkit sort --quiet -s -i --two-pass ~/projects/rna-seq-host-virus/recreate_upenn/CCMP1516/Emiliania_huxleyi.Emiliana_huxleyi_CCMP1516_main_genome_assembly_v1.0.dna.toplevel.fa.gz > sorted.dna.v60.fa
seqkit sort --quiet -s -i --two-pass ~/projects/rna-seq-host-virus/recreate_upenn/CCMP1516/Emiliania_huxleyi.Emiliana_huxleyi_CCMP1516_main_genome_assembly_v1.0.ncrna.fa.gz |> sorted.ncrna.v60.fa
seqkit sort --quiet -s -i --two-pass ~/projects/rna-seq-host-virus/recreate_upenn/CCMP1516/Emiliania_huxleyi.Emiliana_huxleyi_CCMP1516_main_genome_assembly_v1.0.cdna.all.fa.gz > sorted.cdna.v60.fa
seqkit sort --quiet -s -i --two-pass ~/projects/rna-seq-host-virus/recreate_upenn/CCMP1516/Emiliania_huxleyi.Emiliana_huxleyi_CCMP1516_main_genome_assembly_v1.0.pep.all.fa.gz  > sorted.pep.v60.fa

seqkit sort --quiet -s -i --two-pass Emiliania_huxleyi.Emiliana_huxleyi_CCMP1516_main_genome_assembly_v1.0.ncrna.fa.gz > sorted.ncrna.v61.fa
seqkit sort --quiet -s -i --two-pass Emiliania_huxleyi.Emiliana_huxleyi_CCMP1516_main_genome_assembly_v1.0.dna.toplevel.fa.gz > sorted.dna.v61.fa
seqkit sort --quiet -s -i --two-pass Emiliania_huxleyi.Emiliana_huxleyi_CCMP1516_main_genome_assembly_v1.0.cdna.all.fa.gz > sorted.cdna.v61.fa
seqkit sort --quiet -s -i --two-pass Emiliania_huxleyi.Emiliana_huxleyi_CCMP1516_main_genome_assembly_v1.0.pep.all.fa.gz > sorted.pep.v61.fa

diff sorted.pep.v60.fa sorted.pep.v61.fa #nodiff
diff sorted.cdna.v60.fa sorted.cdna.v60.fa #nodiff
diff sorted.dna.v60.fa sorted.dna.v61.fa #nodiff
diff sorted.ncrna.v60.fa sorted.ncrna.v61.fa --suppress-common-lines> ncrnadiffs.txt
wc ncrnadiffs.txt #355
grep ">" ncrnadiffs.txt |wc #178



