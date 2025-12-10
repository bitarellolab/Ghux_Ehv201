#!/bin/bash

# Should be root directory of this repo
datapath='/home/bbitarello/scratch'
mypath=$(pwd)
mypath2=$(pwd)/salmon/salmonq_all

mkdir -p ${mypath2}

# Run j8: host for decoy and both transcriptomes combined (minus two transcripts manually removed from cds file: AET98196.1_308 and AET98250.1_362 (see seq. core notes).
# Save new file as cds_from_genomic_fixID.fna

mkdir -p  ${mypath2}/j_decoy8/

# transcriptome 1: host cDNA
zcat ${datapath}/CCMP1516/Emiliania_huxleyi.Emiliana_huxleyi_CCMP1516_main_genome_assembly_v1.0.cdna.all.fa.gz > ${mypath2}/j_decoy8/gentrome.fa
# transcriptome 2: host ncRNA
zcat ${datapath}/CCMP1516/Emiliania_huxleyi.Emiliana_huxleyi_CCMP1516_main_genome_assembly_v1.0.ncrna.fa.gz >> ${mypath2}/j_decoy8/gentrome.fa
# transcriptome 3: virus 
cat ${datapath}/EhV201/cds_from_genomic_fixID.fna >> ${mypath2}/j_decoy8/gentrome.fa
# genome 1: host genome
zcat ${datapath}/CCMP1516/Emiliania_huxleyi.Emiliana_huxleyi_CCMP1516_main_genome_assembly_v1.0.dna.toplevel.fa.gz >> ${mypath2}/j_decoy8/gentrome.fa

sed -i.bak -e 's/lcl|//g' ${mypath2}/j_decoy8/gentrome.fa


wc ${mypath2}/j_decoy8/gentrome.fa # 3609593 lines
gzip ${mypath2}/j_decoy8/gentrome.fa

#t: transcripts
#d: decoys
#i: index
#p: thead
#k: k-mer length

salmon index --gencode \
             -t ${mypath2}/j_decoy8/gentrome.fa.gz \
             -d ${mypath2}/h_decoy1/decoys.txt \
             -p 12 \
             -k 29 \
             -i ${mypath2}/j_decoy8/decoyIndex > ${mypath2}/j_decoy8/decoy8.out


