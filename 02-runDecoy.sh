#!/bin/bash

# Goals:

* prepare decoy files for salmon quantification
* use host for decoy and both transcriptomes combined
* run salmon index

# Should be root directory of this repo
datapath='/home/bbitarello/scratch'
#mypath=$(pwd)
mypath2=${datapath}/salmonI

mkdir -p ${mypath2}

#  Use host for decoy and both transcriptomes combined

## 1. Make one single transcriptome file + host genome as decoy.
 
#clean up
rm ${mypath2}/gentrome*

# transcriptome 1: host cDNA
zcat ${datapath}/CCMP1516/v62/Emiliania_huxleyi.Emiliana_huxleyi_CCMP1516_main_genome_assembly_v1.0.cdna.all.fa.gz > ${mypath2}/gentrome1.fa
wc -l ${mypath2}/gentrome1.fa #795311
# transcriptome 2: host ncRNA
zcat ${datapath}/CCMP1516/v62/Emiliania_huxleyi.Emiliana_huxleyi_CCMP1516_main_genome_assembly_v1.0.ncrna.fa.gz >> ${mypath2}/gentrome2.fa
wc -l ${mypath2}/gentrome2.fa #3179
# transcriptome 3: virus 
zcat ${datapath}/EhV201/cds_from_genomic.fna.gz >> ${mypath2}/gentrome3.fa
wc -l ${mypath2}/gentrome3.fa # 5217
# genome 1: host genome
zcat ${datapath}/CCMP1516/v62/Emiliania_huxleyi.Emiliana_huxleyi_CCMP1516_main_genome_assembly_v1.0.dna.toplevel.fa.gz >> ${mypath2}/gentrome4.fa
wc -l ${mypath2}/gentrome4.fa #2806184

cat ${mypath2}/gentrome1.fa ${mypath2}/gentrome2.fa ${mypath2}/gentrome3.fa ${mypath2}/gentrome4.fa > ${mypath2}/gentrome5.fa
wc -l ${mypath2}/gentrome5.fa
wc -l ${mypath2}/gentrome5.fa # 3609891=5217+3179+2806184+795311, all good

sed -i.bak -e 's/lcl|//g' ${mypath2}/gentrome5.fa

gzip ${mypath2}/gentrome5.fa

#2. Make decoy file
zgrep "^>" ${datapath}/CCMP1516/v62/Emiliania_huxleyi.Emiliana_huxleyi_CCMP1516_main_genome_assembly_v1.0.dna.toplevel.fa.gz \
   | cut -d " " -f 1 > ${mypath2}/decoys.txt
sed -i.bak -e 's/>//g' ${mypath2}/decoys.txt
wc  -l ${mypath2}/decoys.txt


#t: transcripts
#d: decoys
#i: index
#p: thead
#k: k-mer length

conda deactivate
conda activate py3.13


## salmon v1.11.4 (05/2026)
mkdir -p ${mypath2}/decoyIndex

salmon index --gencode \
             -t ${mypath2}/gentrome5.fa.gz \
             -d ${mypath2}/decoys.txt \
             -p 12 \
             -k 29 \
             -i ${mypath2}/decoyIndex > ${mypath2}/decoy.out


